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
knn_k <- 50
segmentation_mode <- "kinematic" # or "path_signature"
model <- "axisymmetric"          # or "bisymmetric_bar" for a known bar
support_mode <- "starlet"         # or "line_flux" for a clean gas-kinematics footprint
line_flux_sigma <- 3              # robust border-noise threshold for `line_flux`

# Keep this empty for the ordinary disc model. For a bar model, Capivara derives
# the prior from white light and never substitutes the disc PA.
model_control <- list()
# model_control <- list(bar_phi_deg = 41, disc_inc_deg = 38) # manual override

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
  model = model,
  output_dir = output_dir,
  knn_k = knn_k,
  n_segments = n_segments,
  n_path_segments = n_path_segments,
  support_mode = support_mode,
  line_flux_sigma = line_flux_sigma,
  model_control = model_control,
  show_plots = TRUE
)

print(result)

# Each entry is a ggplot object. Edit or save individual panels without
# rerunning the cube analysis.
panels <- kinematic_panels(result, view = "model")
print(panels$velocity)

panel_dir <- file.path(output_dir, "individual_panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
for (name in names(panels)) {
  ggplot2::ggsave(
    filename = file.path(panel_dir, paste0(name, ".png")),
    plot = panels[[name]],
    width = 5,
    height = 4,
    dpi = 320,
    bg = "white"
  )
}

# For `model = "bisymmetric_bar"`, use this instead to export the separate
# circular and bar-flow component maps:
# component_panels <- kinematic_panels(result, view = "components")

message("Saved figures and data in: ", result$output_dir)
