#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (identical(Sys.getenv("CAPIVARA_DEV_LOAD", unset = "false"), "true") &&
      requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(Sys.getenv("CAPIVARA_DEV_ROOT", unset = getwd()), quiet = TRUE)
  } else {
    library(capivara)
  }
  library(FITSio)
  library(ggplot2)
  library(patchwork)
  library(reshape2)
  library(scales)
})

out_root <- Sys.getenv(
  "CAPIVARA_STRUCTURAL_OUT",
  unset = "/private/tmp/capivara2_structural_awareness_showcase"
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  list(
    id = "manga8482_9101_ring_like",
    label = "MaNGA 8482-9101 ring-like",
    path = Sys.getenv(
      "CAPIVARA_STRUCTURAL_RING_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/manga-8482-9101-LOGCUBE.fits"
    ),
    thin = as.integer(Sys.getenv("CAPIVARA_STRUCTURAL_THIN", unset = "20")),
    ncomp = as.integer(Sys.getenv("CAPIVARA_STRUCTURAL_NCOMP_RING", unset = "25"))
  ),
  list(
    id = "manga8083_12704_barred",
    label = "MaNGA 8083-12704 barred/asymmetric",
    path = Sys.getenv(
      "CAPIVARA_STRUCTURAL_BAR_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/manga-8083-12704-LOGCUBE.fits"
    ),
    thin = as.integer(Sys.getenv("CAPIVARA_STRUCTURAL_THIN", unset = "20")),
    ncomp = as.integer(Sys.getenv("CAPIVARA_STRUCTURAL_NCOMP_BAR", unset = "25"))
  )
)

vangogh <- c(
  "#17223B", "#244A7C", "#2F77A4", "#56A0A6", "#8AC5B5",
  "#F3D17C", "#F6A23A", "#C95F2E", "#8D2F1F", "#4B1E2B",
  "#6B559F", "#2A7F74", "#D8B365", "#5A3E2B"
)
as_df <- function(mat, value = "value") {
  df <- reshape2::melt(mat, varnames = c("x", "y"), value.name = value)
  names(df)[3] <- value
  df
}

plot_num <- function(mat, title = NULL, cols = vangogh, legend = FALSE) {
  ggplot(as_df(mat), aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(colours = cols, na.value = "white", oob = scales::squish) +
    theme_void(base_size = 12) +
    theme(
      legend.position = if (legend) "right" else "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.02),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    labs(title = title, fill = NULL)
}

plot_mask <- function(mask, title = NULL, fill = "#C95F2E") {
  plot_num(ifelse(mask, 1, NA_real_), title = title, cols = c(fill, fill))
}

plot_seg <- function(map, title = NULL) {
  df <- as_df(map, "bin")
  df$bin <- factor(df$bin)
  pal <- grDevices::colorRampPalette(vangogh)(max(2, length(levels(stats::na.omit(df$bin)))))
  ggplot(df, aes(x, y, fill = bin)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_manual(values = pal, na.value = "white") +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.02),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    labs(title = title)
}

plot_components <- function(catalogue) {
  plot_seg(catalogue$component_map, "neutral structure components")
}

read_manga_cube <- function(path, thin = 20L) {
  flux <- FITSio::readFITS(path, hdu = 1)$imDat
  wave <- as.numeric(FITSio::readFITS(path, hdu = 6)$imDat)
  keep <- wave >= 3700 & wave <= 9000
  flux <- flux[, , keep, drop = FALSE]
  wave <- wave[keep]
  flux[!is.finite(flux)] <- NA_real_

  if (thin > 1L) {
    groups <- split(seq_along(wave), ceiling(seq_along(wave) / thin))
    binned <- array(NA_real_, dim = c(dim(flux)[1], dim(flux)[2], length(groups)))
    wout <- numeric(length(groups))
    for (k in seq_along(groups)) {
      idx <- groups[[k]]
      binned[, , k] <- apply(
        flux[, , idx, drop = FALSE],
        c(1, 2),
        function(z) if (all(!is.finite(z))) NA_real_ else mean(z, na.rm = TRUE)
      )
      wout[k] <- mean(wave[idx], na.rm = TRUE)
    }
    flux <- binned
    wave <- wout
  }

  list(imDat = flux, hdr = list(source = basename(path)), axDat = list(wave = wave))
}

run_case <- function(case) {
  if (!file.exists(case$path)) {
    warning("Skipping missing cube: ", case$path)
    return(NULL)
  }

  message("Running structural-awareness showcase: ", case$id)
  out_dir <- file.path(out_root, case$id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  input <- read_manga_cube(case$path, thin = case$thin)
  support <- detect_support(
    input,
    method = "starlet",
    starlet_J = 5,
    starlet_scales = 2:5,
    include_coarse = FALSE,
    denoise_k = 1.2,
    mode = "soft",
    positive_only = TRUE
  )
  scores <- score_structures(
    input,
    support = support,
    structure_scales = 3:4,
    ridge_weight = 0.35,
    starlet_weight = 0.65
  )
  scores <- threshold_structures(scores, mode = "extended", low_quantile = 0.20)
  catalogue <- catalogue_structures(scores)
  segmentation <- segment_structures(
    input,
    scores,
    Ncomp = case$ncomp,
    feature_maps = c("combined_structure_score", "starlet_score", "ridge_score"),
    feature_repeats = 5L,
    knn_k = 50,
    auto_k = TRUE,
    max_k = 80,
    verbose = FALSE
  )

  utils::write.csv(
    catalogue$catalogue,
    file.path(out_dir, paste0(case$id, "_structural_catalogue.csv")),
    row.names = FALSE
  )

  panel <- (
    plot_num(scores$maps$collapsed, "whole-spectrum collapse") |
    plot_num(scores$maps$combined_structure_score, "combined structure score") |
      plot_num(scores$maps$ridge_score, "ridge score")
  ) / (
    plot_mask(scores$masks$structure_mask, "extended structure mask") |
      plot_components(catalogue) |
      plot_seg(segmentation$cluster_map, "structure-aware segmentation")
  ) +
    plot_annotation(title = case$label)

  ggsave(
    file.path(out_dir, paste0(case$id, "_structural_awareness_panel.png")),
    panel,
    width = 13,
    height = 7.2,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(case$id, "_structure_aware_segmentation.png")),
    plot_seg(segmentation$cluster_map, "structure-aware segmentation"),
    width = 5.2,
    height = 5,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(case$id, "_neutral_structure_components.png")),
    plot_components(catalogue),
    width = 5.2,
    height = 5,
    dpi = 320,
    bg = "white"
  )

  saveRDS(
    list(case = case, input = input, support = support, scores = scores,
         catalogue = catalogue, segmentation = segmentation),
    file.path(out_dir, paste0(case$id, "_structural_awareness_result.rds"))
  )

  data.frame(
    case = case$id,
    support_px = sum(scores$support_mask, na.rm = TRUE),
    structure_px = sum(scores$masks$structure_mask, na.rm = TRUE),
    n_catalogue = nrow(catalogue$catalogue),
    n_bins = length(unique(stats::na.omit(as.vector(segmentation$cluster_map)))),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, Filter(Negate(is.null), lapply(cases, run_case)))
utils::write.csv(summary, file.path(out_root, "structural_awareness_showcase_summary.csv"), row.names = FALSE)
print(summary)
message("Wrote structural-awareness showcase to: ", out_root)
