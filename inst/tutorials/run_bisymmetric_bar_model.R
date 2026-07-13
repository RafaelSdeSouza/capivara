#!/usr/bin/env Rscript

# Capivara bisymmetric bar model: one known barred galaxy
#
# Edit only this block, then click Source in RStudio. The bar angle must come
# from imaging. This script does not try to infer it from the disc PA.

cube_path <- "/path/to/manga-8078-12703-LOGCUBE.fits"
redshift <- NA_real_
emission_line <- "halpha"
bar_phi_deg <- 41 # in-plane angle relative to the disc major axis

output_dir <- file.path(
  dirname(cube_path),
  "capivara_outputs",
  tools::file_path_sans_ext(basename(cube_path)),
  "bisymmetric_bar"
)

# Nothing below this line needs editing.

if (!file.exists(cube_path)) {
  stop("Set `cube_path` to an existing FITS cube.", call. = FALSE)
}
if (!is.finite(bar_phi_deg)) {
  stop("Set `bar_phi_deg` from an imaging measurement.", call. = FALSE)
}

suppressPackageStartupMessages(library(capivara))

result <- run_manga_bar_model(
  cube_path = cube_path,
  redshift = redshift,
  emission_line = emission_line,
  segmentation_mode = "kinematic",
  bar_phi_deg = bar_phi_deg,
  output_dir = output_dir,
  knn_k = 50,
  n_segments = 25,
  show_plots = TRUE
)

print(result)
message("Saved figures and data in: ", result$output_dir)
