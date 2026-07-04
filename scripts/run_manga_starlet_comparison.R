args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/run_manga_starlet_comparison.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

fits_path <- if (length(args) >= 1) args[[1]] else file.path(repo_root, "..", "sagui_capivara_MaNGA", "manga-8140-12703-LOGCUBE.fits")
png_path <- if (length(args) >= 2) args[[2]] else "/tmp/capivara_manga_starlet_comparison.png"
rds_path <- if (length(args) >= 3) args[[3]] else "/tmp/capivara_manga_starlet_comparison.rds"

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE)
} else {
  library(capivara)
}

palette_van_gogh_div <- function(n = 256) {
  stops <- c(
    "#80B7FF",
    "#547FFF",
    "#405CFF",
    "#263C8B",
    "#FFFAA3",
    "#FFDE38",
    "#BFA524"
  )
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}

matrix_to_df <- function(mat, panel) {
  df <- reshape2::melt(mat, varnames = c("Row", "Col"), value.name = "value")
  df$panel <- panel
  df
}

mat_to_df_sagui <- function(mat, panel) {
  mat2 <- mat
  mat2[mat2 <= 0] <- NA_real_
  matrix_to_df(mat2, panel)
}

normalize_panel <- function(v) {
  s <- asinh(v)
  rng <- range(s[is.finite(s)], na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(v)))
  }
  (s - rng[1]) / diff(rng)
}

x <- FITSio::readFITS(fits_path)
cube <- x$imDat

starlet_cfg <- list(
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  positive_only = TRUE
)

ncomp <- as.integer(Sys.getenv("CAPIVARA_STARLET_COMPARISON_NCOMP", unset = "25"))
knn_k <- as.integer(Sys.getenv("CAPIVARA_STARLET_COMPARISON_KNN", unset = "100"))
feature_window_min <- Sys.getenv("CAPIVARA_STARLET_FEATURE_MIN", unset = "")
feature_window_max <- Sys.getenv("CAPIVARA_STARLET_FEATURE_MAX", unset = "")
feature_window <- NULL
if (nzchar(feature_window_min) && nzchar(feature_window_max)) {
  feature_window <- c(as.numeric(feature_window_min), as.numeric(feature_window_max))
}

base_res <- segment_large(
  x,
  Ncomp = ncomp,
  feature_wavelength_range = feature_window,
  knn_k = knn_k,
  auto_k = FALSE,
  max_k = max(200, knn_k),
  verbose = FALSE
)
star_res <- do.call(
  segment_large,
  c(
    list(
      input = x,
      Ncomp = ncomp,
      use_starlet_mask = TRUE,
      mask_mode = "na",
      feature_wavelength_range = feature_window,
      knn_k = knn_k,
      auto_k = FALSE,
      max_k = max(200, knn_k),
      verbose = FALSE
    ),
    starlet_cfg
  )
)

collapsed <- collapse_white_light(cube)

white_df <- mat_to_df_sagui(collapsed, "White-light")
white_df$value_norm <- normalize_panel(white_df$value)

starlet_mask <- star_res$support_info$mask
mask_df <- matrix_to_df(starlet_mask * 1, "Starlet Mask")

white_plot <- ggplot2::ggplot(white_df, ggplot2::aes(x = Row, y = Col, fill = value_norm)) +
  ggplot2::geom_raster() +
  ggplot2::coord_fixed() +
  ggplot2::scale_fill_gradientn(
    colours = palette_van_gogh_div(256),
    limits = c(0, 1),
    na.value = "black"
  ) +
  ggplot2::labs(title = "White-light Full Frame", fill = NULL) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
  )

mask_plot <- ggplot2::ggplot(mask_df, ggplot2::aes(x = Row, y = Col, fill = factor(value))) +
  ggplot2::geom_raster() +
  ggplot2::coord_fixed() +
  ggplot2::scale_fill_manual(values = c("0" = "black", "1" = "#FFDE38"), na.value = "black") +
  ggplot2::labs(
    title = sprintf(
      "Starlet support mask (%d spaxels)",
      sum(starlet_mask, na.rm = TRUE)
    ),
    fill = NULL
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "none",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
  )

cluster_palette <- palette_van_gogh_div(ncomp)
base_plot <- plot_cluster(base_res, palette = cluster_palette) +
  ggplot2::ggtitle(sprintf("segment_large()\n%d valid spaxels", sum(!is.na(base_res$cluster_map)))) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15))

star_plot <- plot_cluster(star_res, palette = cluster_palette) +
  ggplot2::ggtitle(sprintf("segment_large(use_starlet_mask=TRUE)\n%d valid spaxels", sum(!is.na(star_res$cluster_map)))) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15))

grDevices::png(png_path, width = 1600, height = 1400, res = 170)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
print(white_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(mask_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
print(base_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
print(star_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
grid::popViewport()
grDevices::dev.off()

saveRDS(
  list(
    base = base_res,
    starlet = star_res,
    mask_source = "full_frame",
    starlet_cfg = starlet_cfg,
    ncomp = ncomp,
    knn_k = knn_k,
    feature_wavelength_range = feature_window,
    dims = dim(cube),
    fits_path = fits_path
  ),
  rds_path
)

cat("Saved comparison PNG to:", png_path, "\n")
cat("Saved result RDS to:", rds_path, "\n")
