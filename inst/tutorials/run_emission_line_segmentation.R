#!/usr/bin/env Rscript

# Capivara emission-line segmentation
# Edit this block, click Source in RStudio, and inspect `seg`.

library(FITSio)
library(capivara)

cube_path <- "/path/to/your/cube.fits"
redshift <- 0.0
output_dir <- file.path(dirname(cube_path), "capivara_outputs", "emission_line_segments")

n_segments <- 50
knn_k <- 30
lines <- "agn"

# Nothing below this line should usually need editing.

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cube <- FITSio::readFITS(cube_path)

seg <- capivara::segment_emission_lines(
  cube,
  redshift = redshift,
  lines = lines,
  Ncomp = n_segments,
  knn_k = knn_k,
  feature_mode = "windows",
  line_window_kms = 500,
  spatial_weight = 0.05,
  verbose = TRUE
)

print(capivara::plot_cluster(seg))

saveRDS(seg, file.path(output_dir, "emission_line_segments.rds"))
FITSio::writeFITSim(
  seg$cluster_map * 1,
  file = file.path(output_dir, "emission_line_segments.fits"),
  type = "double"
)

png(file.path(output_dir, "emission_line_segments.png"), width = 1500, height = 1500, res = 220)
print(capivara::plot_cluster(seg))
dev.off()

message("Saved outputs in: ", output_dir)
