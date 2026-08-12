#!/usr/bin/env Rscript

# Cen A MAGNUM Capivara comparison
#
# RStudio use: edit only the settings block, click Source, and inspect the
# output folder. The workflow compares plain moment-map support, starlet support,
# and a hybrid support that is usually safer for clumpy emission-line cubes.

# ---- Settings ---------------------------------------------------------------

cube_path <- Sys.getenv(
  "CAPIVARA_CUBE_PATH",
  unset = "/Users/rd23aag/Documents/GitHub/HUB_2026/Cecilia/CenA/CentA_magnum.fits"
)
output_dir <- Sys.getenv(
  "CAPIVARA_OUTPUT_DIR",
  unset = file.path(dirname(cube_path), "capivara_outputs", "CentA_magnum", "starlet_vs_plain_emission_n80")
)

n_segments <- as.integer(Sys.getenv("CAPIVARA_NCOMP", unset = "80"))
knn_k <- as.integer(Sys.getenv("CAPIVARA_KNN", unset = "50"))
continuum_k <- as.integer(Sys.getenv("CAPIVARA_CONTINUUM_K", unset = "151"))
n_features <- as.integer(Sys.getenv("CAPIVARA_EMISSION_FEATURES", unset = "300"))

moment_quantile <- as.numeric(Sys.getenv("CAPIVARA_MOMENT_QUANTILE", unset = "0.70"))
min_component_pixels <- as.integer(Sys.getenv("CAPIVARA_MIN_COMPONENT_PIXELS", unset = "4"))
spatial_weight <- as.numeric(Sys.getenv("CAPIVARA_SPATIAL_WEIGHT", unset = "0.20"))

# Small scales matter here: the moment map shows compact knots and filaments
# that broad white-light starlet support can suppress.
starlet_scales <- as.integer(strsplit(Sys.getenv("CAPIVARA_STARLET_SCALES", unset = "1,2,3,4"), ",")[[1]])
include_coarse_starlet <- FALSE

# Adaptive support is the recommended sky-cleaning gate for this cube: it asks
# whether emission-like signal persists across selected residual channels after
# robust border-sky normalization.
adaptive_z_threshold <- as.numeric(Sys.getenv("CAPIVARA_ADAPTIVE_Z", unset = "2.5"))
adaptive_min_persistence <- as.integer(Sys.getenv("CAPIVARA_ADAPTIVE_PERSISTENCE", unset = "2"))
adaptive_single_band_z <- as.numeric(Sys.getenv("CAPIVARA_ADAPTIVE_SINGLE_Z", unset = "5.5"))

# ---- Setup ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(patchwork)
  library(pkgload)
})

capivara_repo <- Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara")
pkgload::load_all(capivara_repo, quiet = TRUE)

if (!file.exists(cube_path)) {
  stop("Cube not found: ", cube_path, call. = FALSE)
}
if (continuum_k %% 2L == 0L) {
  continuum_k <- continuum_k + 1L
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

asinh_stretch <- function(x) {
  finite <- is.finite(x)
  if (!any(finite)) return(x)
  scale <- stats::quantile(abs(x[finite]), 0.99, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- max(abs(x[finite]), na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- 1
  asinh(x / scale)
}

impute_linear <- function(x) {
  good <- which(is.finite(x))
  if (length(good) < 2L) return(rep(NA_real_, length(x)))
  stats::approx(good, x[good], xout = seq_along(x), rule = 2)$y
}

positive_continuum_residual <- function(flux, k) {
  flux <- impute_linear(as.numeric(flux))
  if (!all(is.finite(flux))) return(rep(0, length(flux)))
  cont <- stats::runmed(flux, k = k, endrule = "median")
  resid <- flux - cont
  resid[!is.finite(resid) | resid < 0] <- 0
  resid
}

matrix_frame <- function(mat) {
  data.frame(
    x = rep(seq_len(nrow(mat)), times = ncol(mat)),
    y = rep(seq_len(ncol(mat)), each = nrow(mat)),
    value = as.vector(mat)
  )
}

plot_image <- function(mat, title, palette = c("#10021f", "#41217a", "#b12a90", "#ffb000", "#fff7bc")) {
  ggplot(matrix_frame(mat), aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(colours = palette, na.value = "black") +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
}

plot_mask <- function(mask, title) {
  mat <- ifelse(mask, 1, NA_real_)
  ggplot(matrix_frame(mat), aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradient(low = "white", high = "white", na.value = "black") +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
}

cluster_palette <- function(n) {
  grDevices::hcl(
    h = seq(10, 370, length.out = n + 1L)[seq_len(n)],
    c = rep(c(95, 70, 110, 55), length.out = n),
    l = rep(c(35, 72, 48, 84, 58), length.out = n),
    fixup = TRUE
  )
}

drop_tiny_components <- function(mask, min_pixels = 4L) {
  mask[is.na(mask)] <- FALSE
  nr <- nrow(mask)
  nc <- ncol(mask)
  seen <- matrix(FALSE, nr, nc)
  keep <- matrix(FALSE, nr, nc)
  dirs <- matrix(c(1, 0, -1, 0, 0, 1, 0, -1), ncol = 2, byrow = TRUE)

  for (idx in which(mask & !seen)) {
    start <- arrayInd(idx, .dim = dim(mask))
    qx <- integer(nr * nc)
    qy <- integer(nr * nc)
    members <- integer(nr * nc)
    head <- 1L
    tail <- 1L
    count <- 0L
    qx[tail] <- start[1]
    qy[tail] <- start[2]
    seen[start[1], start[2]] <- TRUE

    while (head <= tail) {
      x <- qx[head]
      y <- qy[head]
      head <- head + 1L
      count <- count + 1L
      members[count] <- x + (y - 1L) * nr
      for (d in seq_len(nrow(dirs))) {
        xx <- x + dirs[d, 1]
        yy <- y + dirs[d, 2]
        if (xx >= 1L && xx <= nr && yy >= 1L && yy <= nc && mask[xx, yy] && !seen[xx, yy]) {
          tail <- tail + 1L
          qx[tail] <- xx
          qy[tail] <- yy
          seen[xx, yy] <- TRUE
        }
      }
    }
    if (count >= min_pixels) {
      keep[members[seq_len(count)]] <- TRUE
    }
  }
  keep
}

run_segmentation <- function(name, support_mask, resid_cube, cube_header, wavelengths) {
  message("Running ", name, " segmentation with ", sum(support_mask), " support pixels")
  seg <- capivara::segment_large(
    list(imDat = resid_cube, hdr = cube_header, axDat = list(wavelength = wavelengths)),
    Ncomp = n_segments,
    use_starlet_mask = FALSE,
    mask = support_mask,
    valid_mode = "signal",
    knn_k = knn_k,
    auto_k = FALSE,
    spatial_weight = spatial_weight,
    feature_scale = "none",
    scale_fn = NULL,
    verbose = TRUE
  )
  seg$centa_magnum <- list(
    support_name = name,
    support_pixels = sum(support_mask),
    cube_path = cube_path,
    continuum_k = continuum_k,
    n_features = n_features,
    selected_wavelengths = wavelengths,
    spatial_weight = spatial_weight
  )

  prefix <- file.path(output_dir, name)
  saveRDS(seg, paste0(prefix, "_result.rds"))
  FITSio::writeFITSim(seg$cluster_map * 1, file = paste0(prefix, "_cluster_map.fits"), type = "double")

  png(paste0(prefix, "_segmentation.png"), width = 1500, height = 1500, res = 220)
  print(capivara::plot_cluster(seg, palette = cluster_palette(seg$Ncomp)))
  dev.off()

  data.frame(
    mode = name,
    requested_segments = n_segments,
    actual_segments = seg$Ncomp,
    support_pixels = sum(support_mask),
    segmented_pixels = sum(is.finite(seg$cluster_map)),
    knn_k = seg$backend_info$knn_k %||% knn_k,
    spatial_weight = spatial_weight,
    result_rds = paste0(prefix, "_result.rds"),
    cluster_map_fits = paste0(prefix, "_cluster_map.fits"),
    segmentation_png = paste0(prefix, "_segmentation.png")
  )
}

# ---- Build emission feature cube -------------------------------------------

message("Reading cube: ", cube_path)
cube <- FITSio::readFITS(cube_path)
nx <- dim(cube$imDat)[1]
ny <- dim(cube$imDat)[2]
nw <- dim(cube$imDat)[3]
wave <- tryCatch(FITSio::axVec(3, cube$axDat), error = function(e) seq_len(nw))

mat <- array(cube$imDat, dim = c(nx * ny, nw))
channel_score <- numeric(nw)
resid_white <- matrix(0, nx, ny)

message("First pass: continuum-subtracted moment map and feature ranking")
for (ii in seq_len(nrow(mat))) {
  resid <- positive_continuum_residual(mat[ii, ], continuum_k)
  channel_score <- channel_score + resid
  resid_white[ii] <- sum(resid, na.rm = TRUE)
  if (ii %% 10000L == 0L) message("  processed ", ii, " / ", nrow(mat), " spaxels")
}

selected <- sort(order(channel_score, decreasing = TRUE)[seq_len(min(n_features, nw))])
selected_wave <- wave[selected]
utils::write.csv(
  data.frame(feature_index = selected, wavelength = selected_wave, score = channel_score[selected]),
  file.path(output_dir, "selected_emission_features.csv"),
  row.names = FALSE
)
saveRDS(
  list(moment_map = resid_white, channel_score = channel_score, selected = selected, wavelength = wave),
  file.path(output_dir, "continuum_subtracted_moment_info.rds")
)

positive_moment <- resid_white[is.finite(resid_white) & resid_white > 0]
moment_threshold <- as.numeric(stats::quantile(positive_moment, moment_quantile, na.rm = TRUE))
plain_mask <- drop_tiny_components(resid_white > moment_threshold, min_component_pixels)

message("Second pass: selected emission-feature cube")
resid_selected <- matrix(0, nrow = nrow(mat), ncol = length(selected))
for (ii in seq_len(nrow(mat))) {
  resid_selected[ii, ] <- positive_continuum_residual(mat[ii, ], continuum_k)[selected]
  if (ii %% 10000L == 0L) message("  selected features for ", ii, " / ", nrow(mat), " spaxels")
}
resid_cube <- array(resid_selected, dim = c(nx, ny, length(selected)))
rm(resid_selected, mat)
gc(verbose = FALSE)

message("Building starlet support on continuum-subtracted emission cube")
starlet_info <- capivara::build_starlet_mask(
  list(imDat = resid_cube, hdr = cube$hdr, axDat = list(wavelength = selected_wave)),
  collapse_fn = function(x) apply(x, c(1, 2), sum, na.rm = TRUE),
  starlet_J = 5,
  starlet_scales = starlet_scales,
  include_coarse = include_coarse_starlet,
  denoise_k = 0,
  mode = "soft",
  positive_only = TRUE
)
starlet_mask <- drop_tiny_components(starlet_info$mask, min_component_pixels)

message("Building adaptive emission support")
adaptive_info <- capivara::build_adaptive_support(
  list(imDat = resid_cube, hdr = cube$hdr, axDat = list(wavelength = selected_wave)),
  bands = seq_len(dim(resid_cube)[3]),
  transform = "asinh",
  sky_method = "border",
  border_fraction = 0.10,
  z_threshold = adaptive_z_threshold,
  min_band_persistence = adaptive_min_persistence,
  single_band_z = adaptive_single_band_z,
  smooth_sigma = 0
)
adaptive_mask <- drop_tiny_components(adaptive_info$mask, min_component_pixels)

hybrid_mask <- drop_tiny_components(plain_mask | adaptive_mask, min_component_pixels)

mask_summary <- data.frame(
  mode = c("plain_moment_no_starlet", "starlet_emission", "adaptive_emission", "hybrid_adaptive_or_moment"),
  pixels = c(sum(plain_mask), sum(starlet_mask), sum(adaptive_mask), sum(hybrid_mask)),
  moment_quantile = moment_quantile,
  min_component_pixels = min_component_pixels,
  starlet_scales = paste(starlet_scales, collapse = ","),
  adaptive_z_threshold = adaptive_z_threshold,
  adaptive_min_persistence = adaptive_min_persistence,
  adaptive_single_band_z = adaptive_single_band_z
)
utils::write.csv(mask_summary, file.path(output_dir, "support_mask_summary.csv"), row.names = FALSE)
FITSio::writeFITSim(plain_mask * 1, file = file.path(output_dir, "plain_moment_no_starlet_support.fits"), type = "double")
FITSio::writeFITSim(starlet_mask * 1, file = file.path(output_dir, "starlet_emission_support.fits"), type = "double")
FITSio::writeFITSim(adaptive_mask * 1, file = file.path(output_dir, "adaptive_emission_support.fits"), type = "double")
FITSio::writeFITSim(hybrid_mask * 1, file = file.path(output_dir, "hybrid_adaptive_or_moment_support.fits"), type = "double")

moment_plot <- plot_image(asinh_stretch(resid_white), "continuum-subtracted moment")
ggsave(file.path(output_dir, "moment0_continuum_residual.png"), moment_plot, width = 5.8, height = 5.8, dpi = 320, bg = "white")
mask_panel <- (
  moment_plot |
    plot_mask(plain_mask, "plain moment support") |
    plot_mask(starlet_mask, "starlet emission support") |
    plot_mask(adaptive_mask, "adaptive emission support") |
    plot_mask(hybrid_mask, "hybrid support")
) + plot_annotation(
  title = "Cen A MAGNUM: support masks before Capivara",
  subtitle = "Hybrid = moment support OR adaptive support; useful when starlet misses compact emission knots."
)
print(mask_panel)
ggsave(file.path(output_dir, "support_mask_comparison.png"), mask_panel, width = 18.5, height = 4.8, dpi = 300, bg = "white")

# ---- Segmentations ----------------------------------------------------------

summaries <- list(
  run_segmentation("plain_moment_no_starlet", plain_mask, resid_cube, cube$hdr, selected_wave),
  run_segmentation("starlet_emission", starlet_mask, resid_cube, cube$hdr, selected_wave),
  run_segmentation("adaptive_emission", adaptive_mask, resid_cube, cube$hdr, selected_wave),
  run_segmentation("hybrid_adaptive_or_moment", hybrid_mask, resid_cube, cube$hdr, selected_wave)
)
summary <- do.call(rbind, summaries)
utils::write.csv(summary, file.path(output_dir, "segmentation_summary.csv"), row.names = FALSE)

overview <- (
  plot_image(asinh_stretch(resid_white), "moment map") |
    plot_mask(plain_mask, "plain support") |
    plot_mask(starlet_mask, "starlet support") |
    plot_mask(adaptive_mask, "adaptive support") |
    plot_mask(hybrid_mask, "hybrid support")
) + plot_annotation(title = "Cen A MAGNUM: support masks used for Capivara")
print(overview)
ggsave(file.path(output_dir, "capivara_support_overview.png"), overview, width = 18.5, height = 4.8, dpi = 260, bg = "white")

notes <- c(
  "Cen A MAGNUM Capivara support comparison",
  "",
  sprintf("Cube: %s", cube_path),
  sprintf("Continuum subtraction: per-spaxel rolling median, k = %d channels", continuum_k),
  sprintf("Emission features selected: %d strongest positive-residual channels", length(selected)),
  sprintf("Plain support: continuum-subtracted moment map > %.2f quantile", moment_quantile),
  sprintf("Starlet support: starlet scales %s on the continuum-subtracted emission cube", paste(starlet_scales, collapse = ",")),
  sprintf("Adaptive support: per-channel border-sky z > %.2f with persistence >= %d, single-channel rescue z > %.2f",
          adaptive_z_threshold, adaptive_min_persistence, adaptive_single_band_z),
  "Hybrid support: union of the plain moment support and adaptive support.",
  "",
  "Interpretation:",
  "Pure broad-band/white-light starlet can miss compact emission-line knots because the continuum dominates the collapsed image.",
  "For this cube, the safer scientific support is line/moment-aware: derive the support from positive continuum residuals, then use adaptive persistence to remove sky.",
  "The hybrid adaptive-or-moment support is the recommended first diagnostic when the moment map clearly contains real small structures."
)
writeLines(notes, file.path(output_dir, "README.txt"))

message("Completed. Outputs in: ", output_dir)
