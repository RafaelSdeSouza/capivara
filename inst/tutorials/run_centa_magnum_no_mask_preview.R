#!/usr/bin/env Rscript

# Cen A MAGNUM no-support-mask preview:
# learn Capivara sparse-Ward regions on a representative sample of valid pixels,
# then assign every valid pixel in the cube to the nearest learned region.
# This avoids imposing a sky/support mask while keeping the run interactive.

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
sample_pixels <- as.integer(Sys.getenv("CAPIVARA_SAMPLE_PIXELS", unset = "18000"))
set.seed(as.integer(Sys.getenv("CAPIVARA_SEED", unset = "20260729")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cluster_palette <- function(n) {
  grDevices::hcl(
    h = seq(10, 370, length.out = n + 1L)[seq_len(n)],
    c = rep(c(95, 70, 110, 55), length.out = n),
    l = rep(c(35, 72, 48, 84, 58), length.out = n),
    fixup = TRUE
  )
}

robust_col_scale <- function(x, center = NULL, scale = NULL) {
  if (is.null(center)) center <- apply(x, 2, stats::median, na.rm = TRUE)
  if (is.null(scale)) {
    scale <- apply(x, 2, stats::IQR, na.rm = TRUE)
    bad <- !is.finite(scale) | scale <= 0
    if (any(bad)) scale[bad] <- apply(x[, bad, drop = FALSE], 2, stats::mad, na.rm = TRUE)
    bad <- !is.finite(scale) | scale <= 0
    scale[bad] <- 1
  }
  list(x = sweep(sweep(x, 2, center, "-"), 2, scale, "/"), center = center, scale = scale)
}

assign_nearest <- function(features, centers, block = 5000L) {
  out <- integer(nrow(features))
  center_norm <- rowSums(centers^2)
  for (start in seq(1L, nrow(features), by = block)) {
    idx <- start:min(nrow(features), start + block - 1L)
    d <- sweep(-2 * tcrossprod(features[idx, , drop = FALSE], centers), 2, center_norm, "+")
    out[idx] <- max.col(-d, ties.method = "first")
  }
  out
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
nx <- dim(cube$imDat)[1]
ny <- dim(cube$imDat)[2]

roi <- cube
roi$imDat <- cube$imDat[, , selected_idx, drop = FALSE]
roi$axDat <- list(wavelength = selected_wave)

mat <- matrix(roi$imDat, nrow = nx * ny, ncol = length(selected_idx))
finite <- which(rowSums(is.finite(mat)) == ncol(mat))
message("Finite pixels with no support mask: ", length(finite))

sample_n <- min(sample_pixels, length(finite))
sample_valid_pos <- sort(sample(seq_along(finite), sample_n))
sample_idx <- finite[sample_valid_pos]

sample_cube <- roi
sample_cube$imDat[] <- NA_real_
sample_mat <- matrix(sample_cube$imDat, nrow = nx * ny, ncol = length(selected_idx))
sample_mat[sample_idx, ] <- mat[sample_idx, , drop = FALSE]
sample_cube$imDat <- array(sample_mat, dim = dim(roi$imDat))

message(
  "Learning sampled no-mask Capivara regions: Ncomp=", n_segments,
  ", k=", knn_k, ", sampled pixels=", sample_n
)
sample_seg <- capivara::segment_large(
  sample_cube,
  Ncomp = n_segments,
  use_starlet_mask = FALSE,
  mask = NULL,
  valid_mode = "finite",
  knn_k = knn_k,
  auto_k = FALSE,
  spatial_weight = 0.20,
  feature_scale = "robust_col",
  scale_fn = NULL,
  verbose = TRUE
)

labels_sample <- as.vector(sample_seg$cluster_map)[sample_idx]
ok_sample <- is.finite(labels_sample)
labels_sample <- as.integer(labels_sample[ok_sample])
sample_idx_ok <- sample_idx[ok_sample]

xy_all <- arrayInd(finite, .dim = c(nx, ny))
xy_sample <- arrayInd(sample_idx_ok, .dim = c(nx, ny))
x_sd <- stats::sd(xy_all[, 2]); if (!is.finite(x_sd) || x_sd <= 0) x_sd <- 1
y_sd <- stats::sd(xy_all[, 1]); if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
xy_features <- cbind(
  x = (xy_all[, 2] - mean(xy_all[, 2])) / x_sd,
  y = (xy_all[, 1] - mean(xy_all[, 1])) / y_sd
)
xy_sample_features <- cbind(
  x = (xy_sample[, 2] - mean(xy_all[, 2])) / x_sd,
  y = (xy_sample[, 1] - mean(xy_all[, 1])) / y_sd
)

scaled <- robust_col_scale(mat[finite, , drop = FALSE])
scaled_sample <- robust_col_scale(mat[sample_idx_ok, , drop = FALSE], scaled$center, scaled$scale)$x

full_features <- cbind(scaled$x, 0.20 * xy_features)
sample_features <- cbind(scaled_sample, 0.20 * xy_sample_features)

centers <- rowsum(sample_features, labels_sample) /
  as.vector(table(factor(labels_sample, levels = sort(unique(labels_sample)))))
centers <- centers[order(as.integer(rownames(centers))), , drop = FALSE]

message("Projecting sampled regions onto all finite no-mask pixels")
assigned <- assign_nearest(full_features, centers)

cluster_map <- matrix(NA_real_, nx, ny)
cluster_map[finite] <- assigned

seg <- sample_seg
seg$cluster_map <- cluster_map
seg$Ncomp <- length(unique(assigned))
seg$centa_magnum <- list(
  mode = "no_support_mask_original_flux_n50_preview",
  source = "original CentA_magnum.fits flux values",
  support = "none supplied; sampled Ward projected to all finite pixels",
  finite_pixels = length(finite),
  sampled_pixels = sample_n,
  selected_features = length(selected_idx),
  knn_k = knn_k
)

prefix <- file.path(output_dir, "no_support_mask_original_flux_n50_preview")
saveRDS(seg, paste0(prefix, "_result.rds"))
FITSio::writeFITSim(seg$cluster_map * 1, file = paste0(prefix, "_cluster_map.fits"), type = "double")

png(paste0(prefix, "_segmentation.png"), width = 1500, height = 1500, res = 220)
print(capivara::plot_cluster(seg, palette = cluster_palette(seg$Ncomp)))
dev.off()

summary <- data.frame(
  mode = "no_support_mask_original_flux_n50_preview",
  requested_segments = n_segments,
  actual_segments = seg$Ncomp,
  support_pixels = NA_integer_,
  finite_pixels = length(finite),
  sampled_pixels = sample_n,
  segmented_pixels = sum(is.finite(seg$cluster_map)),
  knn_k = knn_k,
  selected_features = length(selected_idx),
  result_rds = paste0(prefix, "_result.rds"),
  cluster_map_fits = paste0(prefix, "_cluster_map.fits"),
  segmentation_png = paste0(prefix, "_segmentation.png")
)
utils::write.csv(summary, paste0(prefix, "_summary.csv"), row.names = FALSE)
print(summary)

message("Saved no-mask preview outputs in: ", output_dir)
