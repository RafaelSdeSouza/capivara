#!/usr/bin/env Rscript

# Cen A MAGNUM: whole-spectral-range segmentation versus emission-window segmentation.
# The whole-spectrum map uses all wavelengths compressed into broad spectral bins
# so the no-mask comparison stays interactive on the full MAGNUM field.

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(patchwork)
  library(pkgload)
})

capivara_repo <- Sys.getenv("CAPIVARA_REPO", unset = "")
if (nzchar(capivara_repo)) pkgload::load_all(capivara_repo, quiet = TRUE) else library(capivara)

cube_path <- Sys.getenv(
  "CAPIVARA_CUBE_PATH",
  unset = ""
)
if (!nzchar(cube_path)) stop("Set the explicit cube-path environment variable before running this tutorial.")
moment_dir <- Sys.getenv(
  "CAPIVARA_MOMENT_DIR",
  unset = ""
)
if (!nzchar(moment_dir)) stop("Set CAPIVARA_MOMENT_DIR explicitly.")
emission_dir <- Sys.getenv(
  "CAPIVARA_EMISSION_DIR",
  unset = ""
)
if (!nzchar(emission_dir)) stop("Set CAPIVARA_EMISSION_DIR explicitly.")
output_dir <- Sys.getenv(
  "CAPIVARA_OUTPUT_DIR",
  unset = ""
)
if (!nzchar(output_dir)) stop("Set CAPIVARA_OUTPUT_DIR explicitly.")

redshift <- as.numeric(Sys.getenv("CAPIVARA_REDSHIFT", unset = "0.00183"))
n_segments <- as.integer(Sys.getenv("CAPIVARA_NCOMP", unset = "50"))
knn_k <- as.integer(Sys.getenv("CAPIVARA_KNN", unset = "20"))
sample_pixels <- as.integer(Sys.getenv("CAPIVARA_SAMPLE_PIXELS", unset = "18000"))
spectral_bins <- as.integer(Sys.getenv("CAPIVARA_SPECTRAL_BINS", unset = "180"))
set.seed(as.integer(Sys.getenv("CAPIVARA_SEED", unset = "20260730")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

segment_palette <- function(n) {
  grDevices::colorRampPalette(
    c("#0b1d3a", "#2048d8", "#46b5d1", "#f7f4d6", "#f2c230", "#8e6f12"),
    space = "Lab"
  )(n)
}

matrix_frame <- function(mat) {
  data.frame(
    x = rep(seq_len(nrow(mat)), times = ncol(mat)),
    y = rep(seq_len(ncol(mat)), each = nrow(mat)),
    value = as.vector(mat)
  )
}

plot_continuous <- function(mat, title) {
  vals <- mat[is.finite(mat)]
  lim <- stats::quantile(vals, c(0.01, 0.995), na.rm = TRUE)
  ggplot(matrix_frame(pmin(pmax(mat, lim[1]), lim[2])), aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(colors = c("#090417", "#25135f", "#6b1f8f", "#c43d59", "#ffb000"), na.value = "white") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
}

plot_segments <- function(mat, title, n = n_segments) {
  df <- matrix_frame(mat)
  df$value <- factor(df$value, levels = seq_len(n))
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_manual(values = segment_palette(n), na.value = "white", drop = FALSE) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
}

robust_col_scale <- function(x, center = NULL, scale = NULL) {
  if (is.null(center)) center <- apply(x, 2, stats::median, na.rm = TRUE)
  if (is.null(scale)) {
    scale <- vapply(seq_len(ncol(x)), function(j) {
      stats::mad(x[, j], center = center[j], constant = 1.4826, na.rm = TRUE)
    }, numeric(1))
    fallback <- stats::median(scale[is.finite(scale) & scale > 0], na.rm = TRUE)
    if (!is.finite(fallback) || fallback <= 0) fallback <- 1
    scale[!is.finite(scale) | scale <= 0] <- fallback
  }
  x <- sweep(x, 2, center, "-")
  x <- sweep(x, 2, scale, "/")
  x[!is.finite(x)] <- 0
  list(x = x, center = center, scale = scale)
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

bin_full_spectrum <- function(mat, n_bins) {
  groups <- cut(seq_len(ncol(mat)), breaks = n_bins, labels = FALSE)
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = max(groups))
  for (g in seq_len(ncol(out))) {
    idx <- which(groups == g)
    out[, g] <- rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  }
  out
}

run_whole_binned <- function(cube) {
  prefix <- file.path(output_dir, paste0("centa_whole_spectrum_binned_n", n_segments))
  fits_file <- paste0(prefix, "_cluster_map.fits")
  if (file.exists(fits_file)) {
    return(FITSio::readFITS(fits_file)$imDat)
  }

  nx <- dim(cube$imDat)[1]
  ny <- dim(cube$imDat)[2]
  nw <- dim(cube$imDat)[3]
  mat <- matrix(cube$imDat, nrow = nx * ny, ncol = nw)
  finite <- which(rowSums(is.finite(mat)) == ncol(mat))
  message("Whole-spectrum valid pixels: ", length(finite))
  message("Compressing whole spectrum to ", spectral_bins, " broad bins")
  binned <- bin_full_spectrum(mat[finite, , drop = FALSE], spectral_bins)
  scaled <- robust_col_scale(binned)$x

  ij <- arrayInd(finite, .dim = c(nx, ny))
  x_sd <- stats::sd(ij[, 2]); if (!is.finite(x_sd) || x_sd <= 0) x_sd <- 1
  y_sd <- stats::sd(ij[, 1]); if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
  xy <- cbind(
    x = (ij[, 2] - mean(ij[, 2])) / x_sd,
    y = (ij[, 1] - mean(ij[, 1])) / y_sd
  )
  full_features <- cbind(scaled, 0.05 * xy)

  sample_n <- min(sample_pixels, nrow(full_features))
  sample_pos <- sort(sample(seq_len(nrow(full_features)), sample_n))
  message("Learning whole-spectrum sampled Ward on ", sample_n, " spaxels")
  fit <- capivara:::.sparse_ward_cluster_matrix(
    full_features[sample_pos, , drop = FALSE],
    Ncomp = n_segments,
    knn_k = knn_k,
    auto_k = FALSE,
    verbose = TRUE
  )
  label_levels <- sort(unique(fit$labels))
  centers <- rowsum(full_features[sample_pos, , drop = FALSE], fit$labels) /
    as.vector(table(factor(fit$labels, levels = label_levels)))
  centers <- centers[match(label_levels, as.integer(rownames(centers))), , drop = FALSE]
  assigned <- label_levels[assign_nearest(full_features, centers)]

  cluster_map <- matrix(NA_real_, nx, ny)
  cluster_map[finite] <- assigned
  FITSio::writeFITSim(cluster_map * 1, file = fits_file, type = "double")
  saveRDS(
    list(
      cluster_map = cluster_map,
      Ncomp = length(unique(assigned)),
      backend_info = list(
        mode = "whole_spectrum_binned",
        spectral_bins = spectral_bins,
        valid_pixels = length(finite),
        sampled_pixels = sample_n,
        knn_k = fit$knn_k
      )
    ),
    paste0(prefix, "_result.rds")
  )
  cluster_map
}

message("Reading cube: ", cube_path)
cube <- FITSio::readFITS(cube_path)
moment <- readRDS(file.path(moment_dir, "continuum_subtracted_moment_info.rds"))$moment_map
whole_map <- run_whole_binned(cube)
emission_map <- FITSio::readFITS(file.path(emission_dir, "centa_emission_windows_n50_cluster_map.fits"))$imDat

panel <- plot_continuous(moment, "zero moment / emission reference") |
  plot_segments(whole_map, "whole spectral range") |
  plot_segments(emission_map, "AGN line windows")
panel <- panel + patchwork::plot_annotation(
  title = "Cen A MAGNUM: why emission windows beat whole-spectrum segmentation",
  subtitle = sprintf("z = %.5f, N = %d, same segment palette; whole spectrum binned into %d broad features", redshift, n_segments, spectral_bins)
)

ggsave(
  file.path(output_dir, "centa_whole_vs_emission_windows_n50.png"),
  panel,
  width = 15,
  height = 5.4,
  dpi = 320,
  bg = "white"
)

writeLines(
  c(
    "Cen A MAGNUM whole-spectrum versus emission-window segmentation",
    sprintf("Cube: %s", cube_path),
    sprintf("Redshift: %.5f", redshift),
    sprintf("Segments: %d", n_segments),
    "Whole-spectrum comparison uses all wavelengths compressed into broad bins for tractability.",
    "Emission-window comparison uses redshifted AGN line windows only.",
    "",
    "Main panel:",
    file.path(output_dir, "centa_whole_vs_emission_windows_n50.png")
  ),
  file.path(output_dir, "README.txt")
)

message("Saved comparison in: ", output_dir)
