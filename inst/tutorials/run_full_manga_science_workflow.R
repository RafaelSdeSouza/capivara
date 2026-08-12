#!/usr/bin/env Rscript

# Capivara full MaNGA science workflow
#
# In RStudio: edit only the settings below, then click Source. Each figure is
# printed in the Plots pane and all products are saved below `output_dir`.

# ---- Settings ---------------------------------------------------------------

cube_path <- Sys.getenv(
  "CAPIVARA_CUBE_PATH",
  unset = "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/normal_bar/manga-8602-12705-LOGCUBE.fits"
)
output_dir <- Sys.getenv(
  "CAPIVARA_OUTPUT_DIR",
  unset = file.path(dirname(cube_path), "capivara_outputs", "manga8602_12705_full_analysis")
)

redshift <- 0.0318            # Use the measured MaNGA redshift for this target.
emission_line <- "oiii5007"  # This cube has broader usable [O III] support than H-alpha.

# Standard Capivara segmentation and pPXF use the same 25 starlet-supported bins.
n_segments <- 25
starlet_scales <- 2:5
include_coarse_starlet <- FALSE
clean_starlet_support <- TRUE  # Removes detached starlet islands, keeps the main galaxy.

# Use every pixel in the connected starlet-main footprint. This preserves the
# entire barred region while excluding detached starlet islands and blank sky.
bar_knn_k <- 50
bar_n_segments <- 25
bar_phi_deg <- NA_real_       # NA estimates a photometric prior from white light.
kinematic_support_mode <- "starlet"
line_flux_sigma <- 1.5        # Used only when `kinematic_support_mode = "line_flux"`.

# Set this to the local capivaraPPXF checkout. It is an optional companion
# package because pPXF itself has a separate licence.
ppxf_repo <- "/Users/rd23aag/Documents/GitHub/capivaraPPXF"
sps_file <- "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/outputs/capivara_ppxf_manga_smoke/sps_models/spectra_emiles_9.0.npz"

# ---- Setup ------------------------------------------------------------------

if (!file.exists(cube_path)) {
  stop("Cube not found: ", cube_path, call. = FALSE)
}
if (!dir.exists(ppxf_repo)) {
  stop("capivaraPPXF checkout not found: ", ppxf_repo, call. = FALSE)
}
if (!file.exists(sps_file)) {
  stop("Local E-MILES template grid not found: ", sps_file, call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
segmentation_dir <- file.path(output_dir, "01_segmentation")
ppxf_dir <- file.path(output_dir, "02_ppxf")
kinematics_dir <- file.path(output_dir, "03_kinematics_bar")
dir.create(segmentation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ppxf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(kinematics_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(patchwork)
  library(pkgload)
})

# Loading source checkouts keeps the workflow convenient while developing. For
# an installed workflow, `library(capivara)` and `library(capivaraPPXF)` work too.
capivara_repo <- Sys.getenv(
  "CAPIVARA_REPO",
  unset = "/Users/rd23aag/Documents/GitHub/capivara"
)
if (!dir.exists(capivara_repo) || !file.exists(file.path(capivara_repo, "DESCRIPTION"))) {
  stop("Capivara checkout not found: ", capivara_repo, call. = FALSE)
}
pkgload::load_all(capivara_repo, quiet = TRUE)
pkgload::load_all(ppxf_repo, quiet = TRUE)

# Keep the Python plotting cache in a writable temporary directory. This avoids
# a first-run Matplotlib cache failure on locked-down workstations.
ppxf_cache_dir <- Sys.getenv(
  "CAPIVARA_MPLCONFIGDIR",
  unset = file.path(tempdir(), "capivara-matplotlib")
)
dir.create(ppxf_cache_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(
  MPLCONFIGDIR = ppxf_cache_dir,
  XDG_CACHE_HOME = dirname(ppxf_cache_dir)
)

if (!capivaraPPXF::check_ppxf()) {
  stop(
    "Python pPXF is not available. Run capivaraPPXF::setup_ppxf() after accepting its licence.",
    call. = FALSE
  )
}

save_plot <- function(plot, path, width, height) {
  ggplot2::ggsave(path, plot, width = width, height = height, dpi = 320, bg = "white")
  plot
}

matrix_frame <- function(mat) {
  out <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  out$value <- as.vector(mat)
  out
}

path_diverging <- c("#12264A", "#28689E", "#83BBD0", "#F8F6F0", "#F2B275", "#D65B41", "#861F32")
path_sequential <- c("#10213F", "#1D5E84", "#2E9094", "#75C9AD", "#D4DD72", "#F2CA4E", "#E68635")

plot_map <- function(mat, title, label, diverging = FALSE, palette = NULL) {
  df <- matrix_frame(mat)
  finite <- df$value[is.finite(df$value)]
  if (length(finite) < 2L) {
    return(ggplot() + theme_void() + labs(title = title))
  }
  limits <- if (diverging) {
    q <- stats::quantile(abs(finite), 0.98, na.rm = TRUE)
    c(-q, q)
  } else {
    stats::quantile(finite, c(0.02, 0.98), na.rm = TRUE)
  }
  zero <- (0 - limits[1]) / diff(limits)
  zero <- max(0, min(1, zero))
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    labs(title = title, fill = label) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.key.width = grid::unit(0.32, "in"),
      legend.key.height = grid::unit(0.11, "in")
    )
  palette <- if (is.null(palette)) {
    if (diverging) path_diverging else path_sequential
  } else {
    palette
  }
  if (diverging) {
    p + scale_fill_gradientn(
      colours = palette,
      values = c(seq(0, zero, length.out = 4), seq(zero, 1, length.out = 4)[-1]),
      limits = limits,
      oob = scales::squish,
      na.value = "white"
    )
  } else {
    p + scale_fill_gradientn(
      colours = palette,
      limits = limits,
      oob = scales::squish,
      na.value = "white"
    )
  }
}

message("Resolving MaNGA metadata")
z_info <- capivara:::resolve_manga_redshift(cube_path, redshift = redshift)
redshift <- z_info$redshift
object_id <- z_info$plateifu

# ---- 1. Standard Capivara segmentation -------------------------------------

message("Reading cube and running standard starlet-supported segmentation")
cube <- FITSio::readFITS(cube_path)
starlet_support <- capivara::build_starlet_mask(
  cube,
  starlet_scales = starlet_scales,
  include_coarse = include_coarse_starlet
)
segmentation_mask <- starlet_support$mask
if (isTRUE(clean_starlet_support)) {
  segmentation_mask <- capivara:::clean_galaxy_support(
    segmentation_mask,
    preserve_input = FALSE
  )
}
segmentation_input <- cube
for (channel in seq_len(dim(segmentation_input$imDat)[3])) {
  image <- segmentation_input$imDat[, , channel]
  image[!segmentation_mask] <- NA_real_
  segmentation_input$imDat[, , channel] <- image
}

segmentation <- capivara::segment(
  segmentation_input,
  Ncomp = n_segments,
  redshift = redshift,
  use_starlet_mask = FALSE
)
# The clustering uses Capivara's normal MaNGA-sized backend on the connected
# starlet support only. pPXF then receives the untouched, full-resolution cube.
segmentation$original_cube <- cube
segmentation$starlet_support <- segmentation_mask
rm(segmentation_input)
gc(verbose = FALSE)

saveRDS(segmentation, file.path(segmentation_dir, "capivara_starlet_segments.rds"))
utils::write.csv(
  data.frame(
    row = row(segmentation$cluster_map),
    column = col(segmentation$cluster_map),
    segment = as.vector(segmentation$cluster_map)
  ),
  file.path(segmentation_dir, "capivara_starlet_segment_labels.csv"),
  row.names = FALSE
)

segments_plot <- capivara::plot_cluster(segmentation, palette = "starry_night") +
  ggtitle(sprintf("%s: Capivara starlet-supported segments", object_id))
print(segments_plot)
save_plot(
  segments_plot,
  file.path(segmentation_dir, "capivara_starlet_segments.png"),
  width = 6.0,
  height = 5.6
)

if (!is.null(starlet_support$collapsed) && !is.null(segmentation_mask)) {
  support_plot <- plot_map(starlet_support$collapsed, "White light", "flux") +
    geom_contour(
      data = transform(matrix_frame(ifelse(segmentation_mask, 1, 0)), z = value),
      aes(x = x, y = y, z = z),
      breaks = 0.5,
      colour = "#45D6D0",
      linewidth = 0.55,
      inherit.aes = FALSE
    ) +
    labs(subtitle = "Cyan: connected starlet support used for segmentation")
  print(support_plot)
  save_plot(
    support_plot,
    file.path(segmentation_dir, "starlet_support_on_white_light.png"),
    width = 6.0,
    height = 5.6
  )
}

# ---- 2. pPXF stellar population and kinematic maps -------------------------

message("Building Capivara-binned spectra and running pPXF")
ppxf_input <- capivaraPPXF::as_ppxf_input(
  segmentation,
  spectrum = "sum",
  redshift = redshift,
  metadata = list(source_cube = cube_path, object_id = object_id)
)
saveRDS(ppxf_input, file.path(ppxf_dir, "capivara_ppxf_input.rds"))

ppxf_population <- capivaraPPXF::fit_ppxf_population(
  ppxf_input,
  sps_file = sps_file,
  redshift = redshift,
  lam_range_rest = c(4800, 7400),
  fwhm_gal = 2.76,
  mdegree = 8,
  quiet = TRUE
)
ppxf_population$quality_flag <- ifelse(
  !ppxf_population$fit_ok,
  "fit_failed",
  ifelse(
    abs(ppxf_population$stellar_vel) >= 990 |
      ppxf_population$stellar_sigma <= 20.1 |
      ppxf_population$stellar_sigma >= 399.9,
    "parameter_at_bound",
    "interior_solution"
  )
)
utils::write.csv(
  ppxf_population,
  file.path(ppxf_dir, "ppxf_population_by_segment.csv"),
  row.names = FALSE
)
utils::write.csv(
  ppxf_population[, c("bin", "n_spaxels", "fit_ok", "quality_flag", "chi2", "message")],
  file.path(ppxf_dir, "ppxf_population_quality_by_segment.csv"),
  row.names = FALSE
)

ppxf_maps <- capivaraPPXF::map_ppxf_results(
  segmentation,
  ppxf_population[, c("bin", "stellar_vel_median0", "stellar_sigma", "mean_log_age", "mean_metal", "ml_r")]
)
saveRDS(ppxf_maps, file.path(ppxf_dir, "ppxf_population_maps.rds"))

ppxf_panel <- (
  plot_map(
    ppxf_maps$stellar_vel_median0,
    "pPXF stellar velocity",
    expression(Delta*v),
    TRUE,
    palette = path_diverging
  ) |
    plot_map(
      ppxf_maps$stellar_sigma,
      "pPXF stellar dispersion",
      expression(sigma),
      palette = c("#112447", "#1E6486", "#51AFA0", "#C7DD71", "#F1C44A")
    )
) / (
  plot_map(
    ppxf_maps$mean_log_age,
    "pPXF light-weighted log age",
    "log age",
    palette = path_sequential
  ) |
    plot_map(
      ppxf_maps$mean_metal,
      "pPXF light-weighted metallicity",
      "[M/H]",
      palette = path_sequential
    )
) + plot_annotation(
  title = sprintf("%s: pPXF on Capivara segments", object_id),
  subtitle = "Check ppxf_population_quality_by_segment.csv before interpreting bound-limited measurements."
)
print(ppxf_panel)
save_plot(
  ppxf_panel,
  file.path(ppxf_dir, "ppxf_population_maps.png"),
  width = 10.5,
  height = 9.0
)

# ---- 3. Native gas kinematics and bisymmetric bar hypothesis ----------------

message("Running native ", emission_line, " kinematics and the automatic-prior bar model")
bar_result <- capivara::run_manga_bar_model(
  cube_path = cube_path,
  redshift = redshift,
  emission_line = emission_line,
  segmentation_mode = "kinematic",
  bar_phi_deg = bar_phi_deg,
  output_dir = kinematics_dir,
  object_id = object_id,
  knn_k = bar_knn_k,
  n_segments = bar_n_segments,
  support_mode = kinematic_support_mode,
  line_flux_sigma = line_flux_sigma,
  model_control = list(
    # Cube arrays are stored as [y, x]; transpose maps them back to sky x/y
    # without rotating or flipping the cube.
    display_orientation = "transpose",
    use_bar_support_mask = TRUE,
    bar_support_width_deg = 25,
    robust_fit = TRUE,
    smooth_lambda = 10,
    second_order_lambda = 25
  ),
  show_plots = FALSE
)
saveRDS(bar_result, file.path(kinematics_dir, "capivara_bisymmetric_bar_result.rds"))
print(bar_result)
print(plot(bar_result, which = "all"))

# Individual ggplot panels make later paper-specific styling straightforward.
for (view in c("model", "components")) {
  panels <- capivara::kinematic_panels(bar_result, view = view)
  panel_dir <- file.path(kinematics_dir, "individual_panels", view)
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(panels)) {
    save_plot(panels[[name]], file.path(panel_dir, paste0(name, ".png")), 5.0, 4.0)
  }
}

manifest <- c(
  sprintf("Object: %s", object_id),
  sprintf("Cube: %s", normalizePath(cube_path)),
  sprintf("Redshift: %.7f (%s)", redshift, z_info$source),
  sprintf("Starlet support: %s (%d pixels)", if (clean_starlet_support) "main connected footprint" else "raw", sum(segmentation_mask)),
  "Segmentation backend: standard Capivara Ward clustering",
  sprintf("Starlet segmentation: %s", file.path(segmentation_dir, "capivara_starlet_segments.rds")),
  sprintf("pPXF table: %s", file.path(ppxf_dir, "ppxf_population_by_segment.csv")),
  sprintf("pPXF quality: %s", file.path(ppxf_dir, "ppxf_population_quality_by_segment.csv")),
  sprintf("Bar result: %s", file.path(kinematics_dir, "capivara_bisymmetric_bar_result.rds")),
  sprintf("Kinematic support mode: %s", kinematic_support_mode),
  sprintf("Display coordinate mapping: transpose ([y, x] cube storage to sky x/y)"),
  sprintf("Connected starlet footprint: %d spaxels", sum(bar_result$model_result$native$support_starlet, na.rm = TRUE)),
  sprintf("Kinematic support: %d spaxels", sum(bar_result$model_result$spaxels$valid, na.rm = TRUE)),
  sprintf("Bar angle source: %s", bar_result$model_result$bar_geometry$bar_status),
  sprintf("Bar angle prior (deg): %.3f", bar_result$model_result$bar_geometry$phi_b_deg),
  sprintf("Mean noncircular amplitude (km/s): %.4f", bar_result$model_result$fit$parameters$mean_V2)
)
writeLines(manifest, file.path(output_dir, "README.txt"))

object_prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(object_id))
output_index <- c(
  "Capivara output index",
  "",
  file.path("01_segmentation", "capivara_starlet_segments.png"),
  file.path("01_segmentation", "starlet_support_on_white_light.png"),
  file.path("02_ppxf", "ppxf_population_maps.png"),
  file.path("02_ppxf", "ppxf_population_by_segment.csv"),
  file.path("02_ppxf", "ppxf_population_quality_by_segment.csv"),
  file.path("03_kinematics_bar", paste0(object_prefix, "_", emission_line, "_capivara_kinematic_panel.png")),
  file.path("03_kinematics_bar", paste0(object_prefix, "_", emission_line, "_bisymmetric_bar_model.png")),
  file.path("03_kinematics_bar", paste0(object_prefix, "_", emission_line, "_bisymmetric_bar_components.png")),
  file.path("03_kinematics_bar", "individual_panels", "model"),
  file.path("03_kinematics_bar", "individual_panels", "components"),
  "README.txt",
  "TALK_NOTES.txt"
)
writeLines(output_index, file.path(output_dir, "OUTPUTS.txt"))

talk_notes <- c(
  sprintf("%s: Capivara kinematic summary", object_id),
  "",
  "1. Capivara finds a connected starlet support in the white-light cube.",
  "   Here that keeps the full central barred footprint while removing detached sky islands.",
  "2. For each supported spaxel, the gas velocity is the flux-weighted centroid",
  sprintf("   of the %s profile after local continuum subtraction.", emission_line),
  "3. Kinematic-aware segments group neighbouring spaxels with similar flux,",
  "   velocity, dispersion, and line-shape proxies. They are a diagnostic, not the fit itself.",
  "4. The axisymmetric model is a smooth rotating disk. The bisymmetric model",
  "   adds a two-sided, bar-like tangential and radial flow component.",
  "5. The fit is performed spaxel by spaxel. The starlet-derived bar direction",
  "   is a photometric prior used to evaluate the bar hypothesis, not a forced detection.",
  "",
  "How to read the result:",
  "- A coherent velocity gradient is the rotating gas disk.",
  "- Compare the disk and bisymmetric residuals: a meaningful bar term should",
  "  lower structured residuals and have non-negligible V2t/V2r amplitudes.",
  "- If the V2 terms remain small, report that this line map does not require",
  "  detectable bar streaming rather than claiming a bar-flow detection.",
  "",
  sprintf("This run: %d connected-starlet spaxels; %d gas-kinematic spaxels; mean |V2| = %.3f km/s.",
          sum(bar_result$model_result$native$support_starlet, na.rm = TRUE),
          sum(bar_result$model_result$spaxels$valid, na.rm = TRUE),
          bar_result$model_result$fit$parameters$mean_V2)
)
writeLines(talk_notes, file.path(output_dir, "TALK_NOTES.txt"))

message("\nCompleted. All products are in: ", output_dir)
