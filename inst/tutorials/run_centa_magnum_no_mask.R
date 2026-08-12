#!/usr/bin/env Rscript

# Cen A MAGNUM quick baseline:
# original cube flux, selected emission-region channels, no starlet mask,
# no supplied support mask. The only screening is Capivara's internal
# finite/signal validity check so empty spectra do not dominate the graph.

suppressPackageStartupMessages({
  library(FITSio)
  library(pkgload)
})

capivara_repo <- Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara")
pkgload::load_all(capivara_repo, quiet = TRUE)

cube_path <- Sys.getenv(
  "CAPIVARA_CUBE_PATH",
  unset = "/Users/rd23aag/Documents/GitHub/HUB_2026/Cecilia/CenA/CentA_magnum.fits"
)
output_dir <- Sys.getenv(
  "CAPIVARA_OUTPUT_DIR",
  unset = "/Users/rd23aag/Documents/GitHub/HUB_2026/Cecilia/CenA/capivara_outputs/CentA_magnum/starlet_vs_plain_emission_n80"
)

n_segments <- as.integer(Sys.getenv("CAPIVARA_NCOMP", unset = "50"))
knn_k <- as.integer(Sys.getenv("CAPIVARA_KNN", unset = "20"))
n_features <- as.integer(Sys.getenv("CAPIVARA_N_FEATURES", unset = "80"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cluster_palette <- function(n) {
  grDevices::hcl(
    h = seq(10, 370, length.out = n + 1L)[seq_len(n)],
    c = rep(c(95, 70, 110, 55), length.out = n),
    l = rep(c(35, 72, 48, 84, 58), length.out = n),
    fixup = TRUE
  )
}

feature_table <- file.path(output_dir, "selected_emission_features.csv")
if (!file.exists(feature_table)) {
  stop("Missing selected_emission_features.csv in: ", output_dir, call. = FALSE)
}

selected <- read.csv(feature_table)
selected <- selected[order(selected$score, decreasing = TRUE), ]
selected <- selected[seq_len(min(n_features, nrow(selected))), ]
selected_idx <- sort(selected$feature_index)
selected_wave <- selected$wavelength[match(selected_idx, selected$feature_index)]

message("Reading original cube: ", cube_path)
cube <- FITSio::readFITS(cube_path)

roi <- cube
roi$imDat <- cube$imDat[, , selected_idx, drop = FALSE]
roi$axDat <- list(wavelength = selected_wave)

message(
  "Running Capivara without starlet/support mask: Ncomp=", n_segments,
  ", k=", knn_k, ", features=", length(selected_idx)
)
seg <- capivara::segment_large(
  roi,
  Ncomp = n_segments,
  use_starlet_mask = FALSE,
  mask = NULL,
  valid_mode = "sagui",
  knn_k = knn_k,
  auto_k = FALSE,
  spatial_weight = 0.20,
  feature_scale = "robust_col",
  scale_fn = NULL,
  verbose = TRUE
)

seg$centa_magnum <- list(
  mode = "no_support_mask_original_flux_n50",
  source = "original CentA_magnum.fits flux values",
  support = "none supplied; valid_mode = sagui",
  selected_features = length(selected_idx),
  knn_k = knn_k
)

prefix <- file.path(output_dir, "no_support_mask_original_flux_n50")
saveRDS(seg, paste0(prefix, "_result.rds"))
FITSio::writeFITSim(seg$cluster_map * 1, file = paste0(prefix, "_cluster_map.fits"), type = "double")

png(paste0(prefix, "_segmentation.png"), width = 1500, height = 1500, res = 220)
print(capivara::plot_cluster(seg, palette = cluster_palette(seg$Ncomp)))
dev.off()

summary <- data.frame(
  mode = "no_support_mask_original_flux_n50",
  requested_segments = n_segments,
  actual_segments = seg$Ncomp,
  support_pixels = NA_integer_,
  segmented_pixels = sum(is.finite(seg$cluster_map)),
  knn_k = knn_k,
  selected_features = length(selected_idx),
  result_rds = paste0(prefix, "_result.rds"),
  cluster_map_fits = paste0(prefix, "_cluster_map.fits"),
  segmentation_png = paste0(prefix, "_segmentation.png")
)
utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
print(summary)

message("Saved no-mask outputs in: ", output_dir)
