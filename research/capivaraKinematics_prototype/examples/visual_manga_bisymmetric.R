#!/usr/bin/env Rscript

# RStudio visual run: edit this block, click Source, see plots in the Plots pane.

repo_root <- "/Users/rd23aag/Documents/GitHub/capivara"
cube_path <- "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/bar_merger/manga-10218-12703-LOGCUBE.fits"

# Leave as NA to use the LOGCUBE header or the bundled local MaNGA DRPall table.
redshift <- 0.0461

emission_line <- "halpha"
segmentation_mode <- "kinematic" # "kinematic", "path_signature", "spectral", or "all"
knn_k <- 50
n_segments <- 25
starlet_scales <- "2:5"
include_coarse_starlet <- FALSE

display_orientation <- "rot90_cw"
disc_pa_deg <- NA_real_
disc_inc_deg <- NA_real_
bar_phi_deg <- NA_real_

# If you supply bar_phi_deg, this draws a two-sided bar-support prior mask.
use_bar_support_mask <- is.finite(bar_phi_deg)
bar_support_width_deg <- 25

robust_fit <- TRUE
smooth_lambda <- 10
second_order_lambda <- 25

# If left as NULL, products are saved beside the cube:
# dirname(cube_path)/capivara_outputs/<plateifu>
output_dir <- NULL

# Nothing below this line should need editing.

source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "capivara_kinematics_utils.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "manga_metadata.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "run_manga_bar_model.R"))

result <- run_manga_bar_model(
  cube_path = cube_path,
  redshift = redshift,
  emission_line = emission_line,
  segmentation_mode = segmentation_mode,
  output_dir = output_dir,
  repo_root = repo_root,
  knn_k = knn_k,
  n_segments = n_segments,
  starlet_scales = starlet_scales,
  include_coarse_starlet = include_coarse_starlet,
  display_orientation = display_orientation,
  disc_pa_deg = disc_pa_deg,
  disc_inc_deg = disc_inc_deg,
  bar_phi_deg = bar_phi_deg,
  use_bar_support_mask = use_bar_support_mask,
  bar_support_width_deg = bar_support_width_deg,
  robust_fit = robust_fit,
  smooth_lambda = smooth_lambda,
  second_order_lambda = second_order_lambda,
  show_plots = TRUE
)

print(result)
