#!/usr/bin/env Rscript

# Capivara kinematics: one cube
#
# Edit only this block, then click Source in RStudio. The figures appear in
# the Plots pane and are saved next to the cube.

cube_path <- "/path/to/manga-8078-12703-LOGCUBE.fits"
redshift <- NA_real_       # Supply z for non-MaNGA data; NA resolves MaNGA locally.
emission_line <- "halpha" # e.g. "halpha", "oiii5007", "hbeta"

n_segments <- 25
n_path_segments <- 45
knn_k <- 100
segmentation_mode <- "all" # "kinematic", "path_signature", "spectral", or "all"

# Leave these as NA to estimate the geometry from the velocity field.
disc_pa_deg <- NA_real_
disc_inc_deg <- NA_real_
bar_phi_deg <- NA_real_

output_dir <- file.path(
  dirname(cube_path),
  "capivara_outputs",
  tools::file_path_sans_ext(basename(cube_path))
)

# Nothing below this line needs editing.

if (!file.exists(cube_path)) {
  stop("Set `cube_path` to an existing FITS cube.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(capivara)
})

result <- run_kinematic_analysis(
  cube_path = cube_path,
  redshift = redshift,
  emission_line = emission_line,
  segmentation_mode = segmentation_mode,
  output_dir = output_dir,
  knn_k = knn_k,
  n_segments = n_segments,
  n_path_segments = n_path_segments,
  disc_pa_deg = disc_pa_deg,
  disc_inc_deg = disc_inc_deg,
  bar_phi_deg = bar_phi_deg,
  show_plots = TRUE
)

print(result)
message("Saved figures and data in: ", result$output_dir)
