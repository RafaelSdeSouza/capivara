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
  "CAPIVARA_STRUCTURE_SUPPORT_OUT",
  unset = "/private/tmp/capivara2_structure_support_segmentation"
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

vangogh <- c(
  "#17223B", "#244A7C", "#2F77A4", "#56A0A6", "#8AC5B5",
  "#F3D17C", "#F6A23A", "#C95F2E", "#8D2F1F", "#4B1E2B",
  "#6B559F", "#2A7F74", "#D8B365", "#5A3E2B"
)

cases <- list(
  list(
    id = "bar_manga8135_12701",
    label = "modified-Ferrer support",
    title = "MaNGA 8135-12701",
    path = Sys.getenv(
      "CAPIVARA_BAR_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/sagui_capivara_MaNGA/manga-8135-12701-LOGCUBE.fits"
    ),
    mode = "bar",
    ncomp = as.integer(Sys.getenv("CAPIVARA_BAR_NCOMP", unset = "8"))
  ),
  list(
    id = "bulge_manga8443_6102",
    label = "central smooth support",
    title = "MaNGA 8443-6102",
    path = Sys.getenv(
      "CAPIVARA_BULGE_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/sagui_capivara_MaNGA/manga-8443-6102-LOGCUBE.fits"
    ),
    mode = "bulge",
    ncomp = as.integer(Sys.getenv("CAPIVARA_BULGE_NCOMP", unset = "8"))
  ),
  list(
    id = "ring_manga8482_9101",
    label = "annular support",
    title = "MaNGA 8482-9101",
    path = Sys.getenv(
      "CAPIVARA_RING_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/manga-8482-9101-LOGCUBE.fits"
    ),
    mode = "ring",
    ncomp = as.integer(Sys.getenv("CAPIVARA_RING_NCOMP", unset = "18"))
  )
)

as_df <- function(mat, value = "value") {
  df <- reshape2::melt(mat, varnames = c("row", "col"), value.name = value)
  names(df)[3] <- value
  df
}

plot_num <- function(mat, title = NULL, cols = vangogh) {
  ggplot(as_df(mat), aes(col, row, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_y_reverse() +
    scale_fill_gradientn(colours = cols, na.value = "white", oob = scales::squish) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
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
  pal <- grDevices::colorRampPalette(vangogh)(
    max(2, length(levels(stats::na.omit(df$bin))))
  )
  ggplot(df, aes(col, row, fill = bin)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_y_reverse() +
    scale_fill_manual(values = pal, na.value = "white") +
    theme_void(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.02),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    labs(title = title)
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

component_label <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  lab <- matrix(0L, nr, nc)
  sizes <- integer()
  cur <- 0L
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (!isTRUE(mask[i, j]) || lab[i, j] != 0L) next
      cur <- cur + 1L
      qi <- i
      qj <- j
      lab[i, j] <- cur
      head <- 1L
      size <- 0L
      while (head <= length(qi)) {
        a <- qi[head]
        b <- qj[head]
        head <- head + 1L
        size <- size + 1L
        for (d in list(c(1L, 0L), c(-1L, 0L), c(0L, 1L), c(0L, -1L))) {
          aa <- a + d[1]
          bb <- b + d[2]
          if (aa >= 1L && aa <= nr && bb >= 1L && bb <= nc &&
              isTRUE(mask[aa, bb]) && lab[aa, bb] == 0L) {
            lab[aa, bb] <- cur
            qi <- c(qi, aa)
            qj <- c(qj, bb)
          }
        }
      }
      sizes[cur] <- size
    }
  }
  list(label = lab, sizes = sizes)
}

drop_small_components <- function(mask, min_area = 12L) {
  cc <- component_label(mask)
  keep <- which(cc$sizes >= min_area)
  out <- cc$label %in% keep
  dim(out) <- dim(mask)
  out
}

largest_center_component <- function(mask, center, radius = 5L) {
  cc <- component_label(mask)
  if (!length(cc$sizes)) return(mask & FALSE)
  rr <- sqrt((row(mask) - center[1])^2 + (col(mask) - center[2])^2)
  labs <- unique(cc$label[mask & rr <= radius])
  labs <- labs[labs > 0]
  if (!length(labs)) labs <- which.max(cc$sizes)
  out <- cc$label == labs[which.max(cc$sizes[labs])]
  dim(out) <- dim(mask)
  out
}

make_support_mask <- function(input, mode) {
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

  rr <- sqrt((row(scores$support_mask) - scores$center[1])^2 +
               (col(scores$support_mask) - scores$center[2])^2)
  r95 <- stats::quantile(rr[scores$support_mask], 0.95, na.rm = TRUE)
  if (!is.finite(r95) || r95 <= 0) r95 <- max(rr[scores$support_mask], na.rm = TRUE)

  if (mode == "bar") {
    bar <- detect_bar(scores = scores, min_ellipticity = 0.35, solid_support = TRUE)
    mask <- bar$bar_mask
    score <- bar$bar_score
    extra <- list(bar = bar)
  } else if (mode == "bulge") {
    score <- scores$maps$central_smooth_score
    vals <- score[scores$support_mask & rr <= 0.45 * r95 & is.finite(score)]
    cut <- stats::quantile(vals, 0.45, na.rm = TRUE)
    mask <- scores$support_mask & rr <= 0.45 * r95 & score >= cut
    mask <- largest_center_component(mask, scores$center, radius = 6L)
    mask <- drop_small_components(mask, min_area = 20L)
    extra <- list()
  } else if (mode == "ring") {
    ring <- detect_ring(
      scores = scores,
      min_inner_radius_fraction = as.numeric(Sys.getenv("CAPIVARA_RING_INNER_FRAC", unset = "0.22")),
      threshold_quantile = as.numeric(Sys.getenv("CAPIVARA_RING_THRESHOLD_Q", unset = "0.05")),
      min_area = 24L
    )
    score <- ring$ring_score
    mask <- ring$ring_mask
    extra <- list(ring = ring)
  } else {
    stop("Unknown mode: ", mode)
  }

  scores$masks <- c(scores$masks, structure_support = list(mask))
  list(scores = scores, mask = mask, score = score, support = support, extra = extra)
}

run_case <- function(case) {
  if (!file.exists(case$path)) {
    warning("Skipping missing cube: ", case$path)
    return(NULL)
  }

  message("Running structure-support segmentation: ", case$id)
  out_dir <- file.path(out_root, case$id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  input <- read_manga_cube(
    case$path,
    thin = as.integer(Sys.getenv("CAPIVARA_STRUCTURE_THIN", unset = "20"))
  )
  support <- make_support_mask(input, case$mode)

  seg <- segment_structures(
    input,
    support$scores,
    Ncomp = min(case$ncomp, sum(support$mask, na.rm = TRUE)),
    mask_name = "structure_support",
    feature_maps = c("combined_structure_score", "ridge_score"),
    feature_repeats = 5L,
    knn_k = 40,
    auto_k = TRUE,
    max_k = 80,
    verbose = FALSE
  )

  panel <- (
    plot_num(support$scores$maps$collapsed, "whole-spectrum collapse") |
      plot_num(support$score, paste(case$label, "score")) |
      plot_mask(support$mask, paste(case$label, "mask"), "#C95F2E") |
      plot_seg(seg$cluster_map, paste(case$label, "segmentation"))
  ) +
    plot_annotation(title = case$title)

  ggsave(
    file.path(out_dir, paste0(case$id, "_support_segmentation_panel.png")),
    panel,
    width = 14,
    height = 4.2,
    dpi = 320,
    bg = "white"
  )
  ggsave(
    file.path(out_dir, paste0(case$id, "_segmentation.png")),
    plot_seg(seg$cluster_map, paste(case$label, "segmentation")),
    width = 5.2,
    height = 5,
    dpi = 320,
    bg = "white"
  )

  saveRDS(
    list(case = case, input = input, support = support, segmentation = seg),
    file.path(out_dir, paste0(case$id, "_result.rds"))
  )

  data.frame(
    case = case$id,
    mode = case$mode,
    support_px = sum(support$mask, na.rm = TRUE),
    requested_ncomp = case$ncomp,
    actual_ncomp = length(unique(stats::na.omit(as.vector(seg$cluster_map)))),
    stringsAsFactors = FALSE
  )
}

summary <- do.call(rbind, Filter(Negate(is.null), lapply(cases, run_case)))
utils::write.csv(summary, file.path(out_root, "structure_support_segmentation_summary.csv"), row.names = FALSE)
print(summary)
message("Wrote structure-support segmentation showcase to: ", out_root)
