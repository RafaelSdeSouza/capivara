#!/usr/bin/env Rscript

# Capivara bisymmetric bar model: one known barred galaxy
#
# This explicit preview derives a white-light bar-angle prior and uses a
# placeholder inclination. Neither is a validated dynamical measurement.

cube_path <- Sys.getenv("CAPIVARA_CUBE_PATH", unset = "")
redshift <- NA_real_
emission_line <- "halpha"
phi_bar_disc_deg <- NA_real_ # Automatic from white light; set a measured in-plane angle to override.

output_dir <- file.path(
  dirname(cube_path),
  "capivara_outputs",
  tools::file_path_sans_ext(basename(cube_path)),
  "bisymmetric_bar"
)
if (!nzchar(cube_path)) stop("Set the explicit cube-path environment variable before running this tutorial.")

# Nothing below this line needs editing.

if (!file.exists(cube_path)) {
  stop("Set `cube_path` to an existing FITS cube.", call. = FALSE)
}
suppressPackageStartupMessages(library(capivara))

result <- run_manga_bar_model(
  cube_path = cube_path,
  redshift = redshift,
  emission_line = emission_line,
  segmentation_mode = "kinematic",
  phi_bar_disc_deg = phi_bar_disc_deg,
  output_dir = output_dir,
  knn_k = 50,
  n_segments = 25,
  model_control = list(analysis_mode = "preview"),
  show_plots = TRUE
)

print(result)

# The overview and components are ordinary ggplot objects, ready for your own
# paper layout.
model_panels <- kinematic_panels(result, view = "model")
component_panels <- kinematic_panels(result, view = "components")
print(component_panels$circular_component)

for (view in c("model", "components")) {
  panels <- if (identical(view, "model")) model_panels else component_panels
  panel_dir <- file.path(output_dir, "individual_panels", view)
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
}

message("Saved figures and data in: ", result$output_dir)
