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
  "CAPIVARA_BAR_OUT",
  unset = "/private/tmp/capivara2_bar_detection_showcase"
)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  list(
    id = "manga8135_12701_bar_candidate",
    label = "MaNGA 8135-12701 bar diagnostic",
    path = Sys.getenv(
      "CAPIVARA_BAR_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/sagui_capivara_MaNGA/manga-8135-12701-LOGCUBE.fits"
    ),
    thin = as.integer(Sys.getenv("CAPIVARA_BAR_THIN", unset = "20"))
  ),
  list(
    id = "manga8482_9101_ring_control",
    label = "MaNGA 8482-9101 ring-like control",
    path = Sys.getenv(
      "CAPIVARA_BAR_CONTROL_CUBE",
      unset = "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/manga-8482-9101-LOGCUBE.fits"
    ),
    thin = as.integer(Sys.getenv("CAPIVARA_BAR_THIN", unset = "20"))
  )
)

vangogh <- c(
  "#17223B", "#244A7C", "#2F77A4", "#56A0A6", "#8AC5B5",
  "#F3D17C", "#F6A23A", "#C95F2E", "#8D2F1F", "#4B1E2B"
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

  message("Running bar diagnostic: ", case$id)
  out_dir <- file.path(out_root, case$id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  input <- read_manga_cube(case$path, thin = case$thin)
  bar <- detect_bar(
    input = input,
    structure_scales = 3:4,
    threshold_quantile = 0.58,
    min_ellipticity = 0.35,
    max_pa_scatter_deg = 20,
    starlet_args = list(
      starlet_J = 5,
      starlet_scales = 2:5,
      include_coarse = FALSE,
      denoise_k = 1.2,
      mode = "soft",
      positive_only = TRUE
    )
  )

  utils::write.csv(
    bar$diagnostics,
    file.path(out_dir, paste0(case$id, "_bar_diagnostics.csv")),
    row.names = FALSE
  )
  utils::write.csv(
    bar$radial_profile,
    file.path(out_dir, paste0(case$id, "_bar_radial_profile.csv")),
    row.names = FALSE
  )

  panel <- (
    plot_num(bar$scores$maps$collapsed, "whole-spectrum collapse") |
      plot_num(bar$scores$maps$combined_structure_score, "combined structure score") |
      plot_num(bar$scores$maps$ridge_score, "ridge score")
  ) / (
    plot_num(bar$bar_score, "Ferrer bar evidence") |
      plot_mask(bar$candidate_mask, "thresholded candidate", "#C95F2E") |
      plot_mask(bar$bar_mask, "modified-Ferrer support", "#244A7C")
  ) +
    plot_annotation(
      title = case$label,
      subtitle = sprintf(
        "confidence %.2f; ellipticity %.2f; PA scatter %.1f deg",
        bar$diagnostics$confidence,
        bar$diagnostics$ellipticity,
        bar$diagnostics$pa_scatter_deg
      )
    )

  ggsave(
    file.path(out_dir, paste0(case$id, "_bar_detection_panel.png")),
    panel,
    width = 13,
    height = 7.2,
    dpi = 320,
    bg = "white"
  )
  saveRDS(bar, file.path(out_dir, paste0(case$id, "_bar_detection_result.rds")))

  cbind(case = case$id, bar$diagnostics)
}

summary <- do.call(rbind, Filter(Negate(is.null), lapply(cases, run_case)))
utils::write.csv(summary, file.path(out_root, "bar_detection_showcase_summary.csv"), row.names = FALSE)
print(summary)
message("Wrote bar-detection showcase to: ", out_root)
