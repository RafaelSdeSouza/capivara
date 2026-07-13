#!/usr/bin/env Rscript

# Capivara full workflow tutorial
# Edit the block below, then click Source in RStudio.
# The script shows each step on screen and saves figures/tables beside the cube.

# ---- 0. User settings -------------------------------------------------------

repo_root <- "/Users/rd23aag/Documents/GitHub/capivara"
cube_path <- Sys.getenv(
  "CAPIVARA_TUTORIAL_CUBE",
  unset = "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/normal_bar/manga-8078-12703-LOGCUBE.fits"
)

object_id <- tools::file_path_sans_ext(basename(cube_path))
redshift <- NA_real_        # NA = try local MaNGA metadata/header where possible
emission_line <- "halpha"   # "halpha" or "oiii"

n_segments <- 25
knn_k <- 100
use_starlet_mask <- TRUE
starlet_scales <- 2:5
include_coarse_starlet <- FALSE

run_bar_model <- TRUE
segmentation_mode_for_bar <- "path_signature" # "kinematic", "path_signature", "spectral", or "all"

output_dir <- file.path(dirname(cube_path), "capivara_tutorial_outputs", object_id)

# ---- 1. Setup ---------------------------------------------------------------

if (!dir.exists(repo_root)) {
  stop("Edit `repo_root`; it does not exist: ", repo_root, call. = FALSE)
}
if (!file.exists(cube_path)) {
  stop("Edit `cube_path`; cube not found: ", cube_path, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
})

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE)
} else {
  stop("Install pkgload or install capivara before running this tutorial.", call. = FALSE)
}

source_extension <- function(path) {
  files <- list.files(path, pattern = "\\.R$", full.names = TRUE)
  for (ff in files) source(ff, local = .GlobalEnv)
}
source_extension(file.path(repo_root, "extensions", "capivaraKinematics", "R"))

line_rest <- switch(
  tolower(emission_line),
  halpha = 6562.80,
  ha = 6562.80,
  oiii = 5006.84,
  o3 = 5006.84,
  stop("Unknown emission_line. Use 'halpha' or 'oiii'.", call. = FALSE)
)

save_plot <- function(plot, filename, width = 7.0, height = 5.0) {
  path <- file.path(output_dir, filename)
  ggsave(path, plot, width = width, height = height, dpi = 320, bg = "white")
  path
}

message("Reading cube: ", cube_path)
cube <- FITSio::readFITS(cube_path)

# ---- 2. Segment the cube ----------------------------------------------------

message("Running Capivara segment_large()...")
seg <- segment_large(
  cube,
  Ncomp = n_segments,
  knn_k = knn_k,
  auto_k = FALSE,
  max_k = knn_k,
  feature_scale = "robust_col",
  use_starlet_mask = use_starlet_mask,
  starlet_scales = starlet_scales,
  include_coarse = include_coarse_starlet,
  valid_mode = "sagui",
  verbose = TRUE
)

p_segments <- plot_cluster(seg, palette = "starry_night") +
  ggtitle(sprintf("%s: Capivara segments (N = %d)", object_id, n_segments))
print(p_segments)
segment_png <- save_plot(p_segments, "01_capivara_segments.png", width = 6.2, height = 5.6)

saveRDS(seg, file.path(output_dir, "01_capivara_segments.rds"))
utils::write.csv(
  data.frame(x = row(seg$cluster_map), y = col(seg$cluster_map), segment = as.vector(seg$cluster_map)),
  file.path(output_dir, "01_capivara_segment_map.csv"),
  row.names = FALSE
)

# ---- 3. Extract region spectra for fitting ---------------------------------

message("Summarizing segment spectra...")
spectra <- summarize_cluster_spectra(seg)

sum_spectra <- as.data.frame(spectra$sum_spectra, check.names = FALSE)
sum_spectra$segment <- spectra$cluster_ids
sum_spectra <- tidyr::pivot_longer(sum_spectra, -segment, names_to = "wavelength", values_to = "flux")
sum_spectra$wavelength <- as.numeric(sum_spectra$wavelength)

p_spectra <- ggplot(sum_spectra, aes(wavelength, flux, group = segment, colour = factor(segment))) +
  geom_line(linewidth = 0.35, alpha = 0.75) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  labs(
    title = "Segment summed spectra",
    subtitle = "These are the spectra passed to pPXF, capivaraPPXF, or another fitter",
    x = "Observed wavelength",
    y = "Summed flux"
  )
print(p_spectra)
spectra_png <- save_plot(p_spectra, "02_segment_summed_spectra.png", width = 8.0, height = 4.8)
utils::write.csv(sum_spectra, file.path(output_dir, "02_segment_summed_spectra.csv"), row.names = FALSE)

# ---- 4. Quick spectral-fit visualization -----------------------------------

quick_line_fit <- function(wave, flux, lambda0, window = 900, continuum_inner = 1000, continuum_outer = 2200) {
  vel <- 299792.458 * (wave / lambda0 - 1)
  line_idx <- abs(vel) <= window
  cont_idx <- abs(vel) >= continuum_inner & abs(vel) <= continuum_outer
  if (sum(line_idx, na.rm = TRUE) < 5L || sum(cont_idx, na.rm = TRUE) < 5L) return(NULL)
  cont <- stats::median(flux[cont_idx], na.rm = TRUE)
  prof <- pmax(flux[line_idx] - cont, 0)
  v <- vel[line_idx]
  if (sum(prof, na.rm = TRUE) <= 0) return(NULL)
  mu <- sum(v * prof, na.rm = TRUE) / sum(prof, na.rm = TRUE)
  sig <- sqrt(sum(prof * (v - mu)^2, na.rm = TRUE) / sum(prof, na.rm = TRUE))
  amp <- max(prof, na.rm = TRUE)
  data.frame(
    velocity = v,
    observed = flux[line_idx],
    model = cont + amp * exp(-0.5 * ((v - mu) / sig)^2),
    continuum = cont,
    centroid_kms = mu,
    sigma_kms = sig
  )
}

lambda0 <- line_rest * ifelse(is.finite(redshift), 1 + redshift, 1)
if (!is.finite(redshift) && grepl("manga-", basename(cube_path), ignore.case = TRUE)) {
  z_info <- try(resolve_manga_redshift(cube_path, redshift = redshift), silent = TRUE)
  if (!inherits(z_info, "try-error") && is.finite(z_info$redshift)) {
    redshift <- z_info$redshift
    lambda0 <- line_rest * (1 + redshift)
  }
}

fit_rows <- list()
wide_sum <- spectra$sum_spectra
wave <- spectra$wavelength
for (i in seq_len(nrow(wide_sum))) {
  ff <- quick_line_fit(wave, as.numeric(wide_sum[i, ]), lambda0)
  if (is.null(ff)) next
  ff$segment <- rownames(wide_sum)[i]
  fit_rows[[length(fit_rows) + 1L]] <- ff
}
fit_df <- if (length(fit_rows)) do.call(rbind, fit_rows) else data.frame()

if (nrow(fit_df)) {
  keep_segments <- names(sort(table(fit_df$segment), decreasing = TRUE))[seq_len(min(9, length(unique(fit_df$segment))))]
  fit_show <- fit_df[fit_df$segment %in% keep_segments, , drop = FALSE]
  p_fit <- ggplot(fit_show, aes(velocity)) +
    geom_line(aes(y = observed), colour = "#071A3F", linewidth = 0.45) +
    geom_line(aes(y = model), colour = "#D34E24", linewidth = 0.65) +
    facet_wrap(~segment, scales = "free_y", ncol = 3) +
    theme_classic(base_size = 10) +
    labs(
      title = paste("Quick", emission_line, "fit visualization by segment"),
      subtitle = "Moment-based Gaussian quicklook; replace this block with pPXF/capivaraPPXF for science fits",
      x = "Velocity relative to line center (km/s)",
      y = "Flux"
    )
  print(p_fit)
  fit_png <- save_plot(p_fit, "03_quick_spectral_fit_examples.png", width = 8.0, height = 6.2)
  utils::write.csv(fit_df, file.path(output_dir, "03_quick_spectral_fit_examples.csv"), row.names = FALSE)
} else {
  warning("Quick line-fit visualization skipped: line not covered or redshift not resolved.")
}

# ---- 5. Kinematics and bar modelling ---------------------------------------

if (isTRUE(run_bar_model)) {
  message("Running Capivara kinematics + bar model...")
  bar <- run_manga_bar_model(
    cube_path = cube_path,
    redshift = redshift,
    emission_line = emission_line,
    segmentation_mode = segmentation_mode_for_bar,
    output_dir = file.path(output_dir, "04_kinematics_bar"),
    object_id = object_id,
    repo_root = repo_root,
    knn_k = knn_k,
    n_segments = n_segments,
    n_path_segments = max(35, n_segments),
    show_plots = FALSE
  )

  print(bar)
  plot(bar, which = "summary")
  plot(bar, which = "model")
  plot(bar, which = "components")
  saveRDS(bar, file.path(output_dir, "04_kinematics_bar_result.rds"))
}

message("Tutorial complete. Outputs saved in: ", output_dir)
message("Key figures:")
message("  ", segment_png)
message("  ", spectra_png)
if (exists("fit_png")) message("  ", fit_png)
