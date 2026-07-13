#!/usr/bin/env Rscript

# Capivara bisymmetric bar model: one known barred galaxy
#
# Edit only this block, then click Source in RStudio. The bar angle must come
# from imaging. This script does not try to infer it from the disc PA.

cube_path <- "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/normal_bar/manga-8078-12703-LOGCUBE.fits"
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
