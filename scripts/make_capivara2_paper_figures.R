#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

repo_dir <- normalizePath(
  Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara"),
  mustWork = TRUE
)
out_dir <- Sys.getenv(
  "CAPIVARA2_FIG_OUT",
  unset = file.path(repo_dir, "outputs", "capivara2_paper_figures", "individual")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paper_fig_dir <- Sys.getenv(
  "CAPIVARA2_PAPER_FIG_DIR",
  unset = "/Users/rd23aag/Documents/GitHub/Capivara2_0_paper/figures_capivara2"
)
copy_to_paper <- tolower(Sys.getenv("CAPIVARA2_COPY_TO_PAPER", unset = "false")) %in% c("1", "true", "yes")
make_starlet <- tolower(Sys.getenv("CAPIVARA2_MAKE_STARLET", unset = "true")) %in% c("1", "true", "yes")
make_kinematics <- tolower(Sys.getenv("CAPIVARA2_MAKE_KINEMATICS", unset = "true")) %in% c("1", "true", "yes")
make_ppxf <- tolower(Sys.getenv("CAPIVARA2_MAKE_PPXF", unset = "true")) %in% c("1", "true", "yes")
make_mosaics <- tolower(Sys.getenv("CAPIVARA2_MAKE_MOSAICS", unset = "false")) %in% c("1", "true", "yes")
include_oiii <- tolower(Sys.getenv("CAPIVARA2_INCLUDE_OIII", unset = "false")) %in% c("1", "true", "yes")

pkgload::load_all(repo_dir, quiet = TRUE)

vangogh_base <- c(
  "#0B132B", "#1C4E80", "#2F74B5", "#58A4B0", "#F2D06B",
  "#F2A541", "#C7771F", "#7A3E12", "#EFE6A4", "#244C89",
  "#1E8078", "#B84A28"
)

vangogh_discrete <- function(n) {
  if (n <= length(vangogh_base)) return(vangogh_base[seq_len(n)])
  grDevices::colorRampPalette(vangogh_base, space = "Lab")(n)
}

velocity_palette <- grDevices::colorRampPalette(
  c("#16345C", "#2C7FB8", "#D7E8E6", "#F2D06B", "#D46A2C", "#7F1D1D"),
  space = "Lab"
)(256)

ratio_palette <- grDevices::colorRampPalette(
  c("#16345C", "#68A6B8", "#F4E6A1", "#D46A2C", "#7F1D1D"),
  space = "Lab"
)(256)

bpt_palette <- c(
  sf = "#2F74B5",
  star_forming = "#2F74B5",
  composite = "#58A4B0",
  transition = "#58A4B0",
  agn_like = "#C7771F",
  seyfert = "#B84A28",
  liner = "#7A3E12",
  inactive = "#B8B8AA",
  unclassified = "#0B132B"
)

map_df <- function(mat, value_name = "value") {
  df <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  df[[value_name]] <- as.vector(mat)
  df
}

assignment_map <- function(assignments, value_col) {
  nx <- max(assignments$x, na.rm = TRUE)
  ny <- max(assignments$y, na.rm = TRUE)
  mat <- matrix(NA_real_, nrow = nx, ncol = ny)
  ok <- is.finite(assignments$x) & is.finite(assignments$y) & is.finite(assignments[[value_col]])
  mat[cbind(assignments$x[ok], assignments$y[ok])] <- assignments[[value_col]][ok]
  mat
}

panel_theme <- function(base_size = 8) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain", size = base_size + 1),
      plot.margin = margin(2, 2, 2, 2),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = unit(9, "pt"),
      legend.key.width = unit(16, "pt")
    )
}

save_clean_plot <- function(plot, stem, width = 2.4, height = 2.4, dpi = 450) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  invisible(c(png = png_file, pdf = pdf_file))
}

plot_discrete_map <- function(mat, label = NULL) {
  df <- map_df(mat)
  df$value <- factor(df$value)
  n <- length(stats::na.omit(unique(df$value)))
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(values = vangogh_discrete(max(n, 1)), na.value = "white", guide = "none", na.translate = FALSE) +
    labs(title = label) +
    panel_theme()
}

plot_continuous_map <- function(mat,
                                label = NULL,
                                palette = velocity_palette,
                                limits = NULL,
                                legend_title = NULL,
                                legend_labels = waiver()) {
  df <- map_df(mat)
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_gradientn(
      colours = palette,
      limits = limits,
      oob = scales::squish,
      na.value = "white",
      name = legend_title,
      labels = legend_labels
    ) +
    labs(title = label) +
    panel_theme()
}

plot_class_map <- function(mat, label = NULL) {
  df <- map_df(mat)
  df$value <- factor(df$value, exclude = NA)
  values <- bpt_palette[names(bpt_palette) %in% levels(df$value)]
  missing_levels <- setdiff(levels(df$value), names(values))
  if (length(missing_levels)) {
    values <- c(values, stats::setNames(vangogh_discrete(length(missing_levels)), missing_levels))
  }
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(
      values = values,
      na.value = "white",
      name = "BPT",
      na.translate = FALSE
    ) +
    labs(title = label) +
    panel_theme()
}

clean_case_label <- function(path) {
  sub("^manga-", "", basename(dirname(path)))
}

kinematic_case <- function(object_dir, line_dir, velocity_col = "np_centroid_median0") {
  rds <- file.path(repo_dir, "outputs", "capivara2_kinematics_prototype", object_dir, line_dir, "capivara2_kinematics_prototype.rds")
  x <- readRDS(rds)
  list(
    object = sub("^manga-", "", object_dir),
    line = line_dir,
    bins = x$segmentation$cluster_map,
    velocity = assignment_map(x$assignments, velocity_col)
  )
}

make_kinematic_panel <- function() {
  cases <- list(
    list(key = "7443_12703_halpha_nii", data = kinematic_case("manga-7443-12703-LOGCUBE", "ifu_halpha_nii_complex")),
    list(key = "10224_6104_halpha_nii", data = kinematic_case("manga-10224-6104-LOGCUBE", "ifu_halpha_nii_complex"))
  )
  if (include_oiii) {
    cases <- c(
      cases,
      list(
        list(key = "7443_12703_oiii5007", data = kinematic_case("manga-7443-12703-LOGCUBE", "ifu_oiii5007")),
        list(key = "10224_6104_oiii5007", data = kinematic_case("manga-10224-6104-LOGCUBE", "ifu_oiii5007"))
      )
    )
  }
  vlim <- c(-220, 220)

  mosaic_plots <- list()
  for (case in cases) {
    bin_plot <- plot_discrete_map(case$data$bins)
    vel_plot <- plot_continuous_map(
      case$data$velocity,
      legend_title = expression(km~s^{-1}),
      limits = vlim
    )
    save_clean_plot(bin_plot, paste0("capivara2_", case$key, "_bins_vangogh"))
    save_clean_plot(vel_plot, paste0("capivara2_", case$key, "_velocity_kms"))
    mosaic_plots <- c(mosaic_plots, list(bin_plot, vel_plot))
  }

  if (make_mosaics) {
    panel <- wrap_plots(mosaic_plots, ncol = 2, guides = "collect") &
      theme(legend.position = "right")
    ggsave(file.path(out_dir, "capivara2_kinematic_feature_panel.png"), panel, width = 4.8, height = 4.8, dpi = 450, bg = "white")
    ggsave(file.path(out_dir, "capivara2_kinematic_feature_panel.pdf"), panel, width = 4.8, height = 4.8, bg = "white")
  }
}

load_ppxf_maps <- function(path) {
  readRDS(file.path(repo_dir, path))$maps
}

make_ppxf_panel <- function() {
  m7443 <- load_ppxf_maps("outputs/capivara_ppxf_manga_smoke/maps/manga_7443_12703_ppxf_emission_laplace_maps.rds")
  m8135 <- load_ppxf_maps("outputs/capivara_ppxf_manga_8135_12701/maps/manga_8135_12701_ppxf_emission_laplace_maps.rds")

  plots <- list(
    "7443_12703_log_nii_ha" = plot_continuous_map(m7443$log_nii_ha_laplace, palette = ratio_palette, limits = c(-1.1, 0.4)),
    "7443_12703_log_oiii_hb" = plot_continuous_map(m7443$log_oiii_hb_laplace, palette = ratio_palette, limits = c(-0.8, 1.2)),
    "7443_12703_bpt_class" = plot_class_map(m7443$bpt_class_laplace),
    "8135_12701_log_nii_ha" = plot_continuous_map(m8135$log_nii_ha_laplace, palette = ratio_palette, limits = c(-1.1, 0.4)),
    "8135_12701_log_oiii_hb" = plot_continuous_map(m8135$log_oiii_hb_laplace, palette = ratio_palette, limits = c(-0.8, 1.2)),
    "8135_12701_bpt_class" = plot_class_map(m8135$bpt_class_laplace)
  )
  for (nm in names(plots)) {
    save_clean_plot(plots[[nm]], paste0("capivara2_ppxf_", nm))
  }

  if (make_mosaics) {
    panel <- wrap_plots(plots, ncol = 3, guides = "collect") &
      theme(legend.position = "right")
    ggsave(file.path(out_dir, "capivara2_ppxf_bpt_panel.png"), panel, width = 7.0, height = 4.1, dpi = 450, bg = "white")
    ggsave(file.path(out_dir, "capivara2_ppxf_bpt_panel.pdf"), panel, width = 7.0, height = 4.1, bg = "white")
  }
}

make_starlet_panel <- function() {
  fits_path <- Sys.getenv(
    "CAPIVARA2_STARLET_FITS",
    unset = "/Users/rd23aag/Documents/GitHub/sagui_capivara_MaNGA/manga-8140-12703-LOGCUBE.fits"
  )
  if (!file.exists(fits_path)) {
    warning("Skipping starlet paper panel; missing FITS file: ", fits_path)
    return(invisible(NULL))
  }
  x <- FITSio::readFITS(fits_path)
  collapse_sagui <- function(cube) collapse_white_light(cube, kclip = 1)
  starlet_cfg <- list(
    starlet_J = 5,
    starlet_scales = 2:5,
    include_coarse = FALSE,
    denoise_k = 2.5,
    positive_only = TRUE
  )
  base_res <- segment_large(x, Ncomp = 8, knn_k = 100, max_k = 200, verbose = FALSE)
  star_res <- do.call(
    segment_large,
    c(list(input = x, Ncomp = 8, use_starlet_mask = TRUE, collapse_fn = collapse_sagui, mask_mode = "na", knn_k = 100, max_k = 200, verbose = FALSE), starlet_cfg)
  )
  white <- collapse_sagui(x$imDat)
  white[white <= 0] <- NA_real_
  white_norm <- asinh(white)
  rng <- range(white_norm, finite = TRUE)
  white_norm <- (white_norm - rng[1]) / diff(rng)

  mask <- star_res$support_info$mask
  storage.mode(mask) <- "numeric"
  mask[mask == 0] <- NA_real_

  plots <- list(
    "white_light" = plot_continuous_map(white_norm, palette = ratio_palette, limits = c(0, 1), legend_title = NULL) +
      theme(legend.position = "none"),
    "starlet_support" = plot_continuous_map(mask, palette = c("#F2D06B", "#F2D06B"), limits = c(0, 1), legend_title = NULL) +
      theme(legend.position = "none"),
    "segment_large_bins" = plot_discrete_map(base_res$cluster_map),
    "starlet_segment_large_bins" = plot_discrete_map(star_res$cluster_map)
  )
  for (nm in names(plots)) {
    save_clean_plot(plots[[nm]], paste0("capivara2_starlet_", nm))
  }

  if (make_mosaics) {
    panel <- wrap_plots(plots, ncol = 4)
    ggsave(file.path(out_dir, "capivara2_starlet_support_panel.png"), panel, width = 7.4, height = 2.1, dpi = 450, bg = "white")
    ggsave(file.path(out_dir, "capivara2_starlet_support_panel.pdf"), panel, width = 7.4, height = 2.1, bg = "white")
  }
}

if (make_starlet) make_starlet_panel()
if (make_kinematics) make_kinematic_panel()
if (make_ppxf) make_ppxf_panel()

if (copy_to_paper && dir.exists(paper_fig_dir)) {
  files <- list.files(out_dir, pattern = "^capivara2_.*[.](png|pdf)$", full.names = TRUE)
  file.copy(files, paper_fig_dir, overwrite = TRUE)
}

message("Wrote paper figures to: ", out_dir)
if (copy_to_paper) {
  message("Copied paper figures to: ", paper_fig_dir)
}
