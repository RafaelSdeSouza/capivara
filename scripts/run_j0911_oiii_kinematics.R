#!/usr/bin/env Rscript

# Minimal J0911 [OIII] kinematics + Capivara segmentation workflow.
# Edit only this block, then click Source in RStudio or run:
# Rscript run_j0911_oiii_kinematics.R

repo_root <- "/Users/rd23aag/Documents/GitHub/capivara"
cube_path <- "/Users/rd23aag/Documents/GitHub/HUB_2026/Emma/paper/J0911_30x30_cube.fits"
output_dir <- "/Users/rd23aag/Documents/GitHub/HUB_2026/Emma/paper"

object_id <- "J0911"
redshift <- 0.2622
line_rest <- 5006.84
line_label <- "[OIII] 5007"

knn_k <- 50
n_segments <- 25
n_path_segments <- 35

# Pure starlet is intentionally not used as the science mask here. The default
# science mask is built from the [OIII] line flux/SNR, with the SED mask saved
# only as a comparison product.
support_mode <- "line_snr" # "line_snr" or "sed"
mask_border_width <- 3
mask_sigma <- 1.0
mask_dilate <- 0
line_snr_min <- 3
line_flux_sigma <- 2
line_flux_quantile <- 0.25

line_window_kms <- 700
peak_search_kms <- 350
centroid_window_kms <- 220
continuum_inner_kms <- 900
continuum_outer_kms <- 2200
path_window_kms <- 700

run_path_signatures <- FALSE
run_disc_model <- FALSE
disc_pa_deg <- NA_real_
disc_inc_deg <- 60

# Nothing below should need editing.

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(gridExtra)
})
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Please install pkgload or install/load capivara before running this script.", call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "galaxy_mask.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "geometry.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "disc_model.R"))

`%||%` <- function(a, b) if (!is.null(a)) a else b

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- paste0(tolower(object_id), "_oiii5007")

hdr_get <- function(hdr, key, default = NA_character_) {
  i <- which(hdr == key)[1]
  if (!length(i) || is.na(i) || i >= length(hdr)) default else hdr[[i + 1L]]
}

read_cube <- function(path) {
  fits <- FITSio::readFITS(path, hdu = 1L, maxLines = 10000)
  stat <- tryCatch(FITSio::readFITS(path, hdu = 2L, maxLines = 10000)$imDat, error = function(e) NULL)
  cube <- fits$imDat
  if (length(dim(cube)) != 3L) {
    stop("Expected a 3-D DATA cube.", call. = FALSE)
  }
  crval <- as.numeric(hdr_get(fits$hdr, "CRVAL3"))
  crpix <- as.numeric(hdr_get(fits$hdr, "CRPIX3", "1"))
  cdelt <- as.numeric(hdr_get(fits$hdr, "CD3_3", hdr_get(fits$hdr, "CDELT3", "1")))
  wave <- crval + (seq_len(dim(cube)[3]) - crpix) * cdelt
  list(cube = cube, stat = stat, wave = wave, header = fits$hdr)
}

border_mask <- function(nx, ny, width = 3L) {
  x <- seq_len(nx)
  y <- seq_len(ny)
  outer(x <= width | x > nx - width, rep(TRUE, ny), `&`) |
    outer(rep(TRUE, nx), y <= width | y > ny - width, `&`)
}

dilate_mask <- function(mask, iterations = 1L) {
  out <- mask
  if (iterations <= 0L) return(out)
  for (it in seq_len(iterations)) {
    old <- out
    for (i in seq_len(nrow(out))) {
      for (j in seq_len(ncol(out))) {
        ii <- max(1, i - 1):min(nrow(out), i + 1)
        jj <- max(1, j - 1):min(ncol(out), j + 1)
        out[i, j] <- any(old[ii, jj], na.rm = TRUE)
      }
    }
  }
  out
}

build_sed_mask <- function(cube, border_width = 3L, sigma = 1.0, dilate = 1L) {
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  nw <- dim(cube)[3]
  border <- border_mask(nx, ny, border_width)
  mat <- matrix(cube, nrow = nx * ny, ncol = nw)
  sky_spec <- apply(mat[as.vector(border), , drop = FALSE], 2L, stats::median, na.rm = TRUE)
  clean <- sweep(mat, 2L, sky_spec, "-")
  white <- matrix(rowSums(pmax(clean, 0), na.rm = TRUE), nrow = nx, ncol = ny)

  b <- white[border & is.finite(white)]
  thr <- stats::median(b, na.rm = TRUE) + sigma * stats::mad(b, constant = 1.4826, na.rm = TRUE)
  support <- white > thr
  if (sum(support, na.rm = TRUE) < 60L) {
    thr <- stats::quantile(white[is.finite(white)], 0.65, na.rm = TRUE)
    support <- white > thr
  }
  support <- better_galaxy_mask(support, close_iterations = 2L, fill_holes = TRUE, preserve_input = TRUE)
  support <- dilate_mask(support, dilate)
  list(mask = support, white = white, sky_spec = sky_spec, threshold = thr)
}

matrix_df <- function(mat) {
  df <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  df$value <- as.vector(mat)
  df
}

robust_limits <- function(x, probs = c(0.02, 0.98), symmetric = FALSE) {
  vals <- x[is.finite(x)]
  if (!length(vals)) return(NULL)
  if (symmetric) {
    lim <- stats::quantile(abs(vals), probs[2], na.rm = TRUE)
    c(-lim, lim)
  } else {
    as.numeric(stats::quantile(vals, probs, na.rm = TRUE))
  }
}

save_map <- function(mat, file, title, palette = NULL, diverging = FALSE, discrete = FALSE, limits = NULL) {
  df <- matrix_df(mat)
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_y_reverse() +
    theme_void(base_size = 11) +
    labs(title = title, fill = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (discrete) {
    df$value <- factor(df$value)
    n <- max(3L, length(unique(stats::na.omit(df$value))))
    p <- ggplot(df, aes(x, y, fill = value)) +
      geom_raster() +
      coord_fixed(expand = FALSE) +
      scale_y_reverse() +
      scale_fill_manual(values = grDevices::colorRampPalette(c("#17335C", "#3F78A8", "#4DAF8E", "#F1D66B", "#C85B3C"), space = "Lab")(n), na.value = "grey93") +
      theme_void(base_size = 11) +
      labs(title = title, fill = "segment") +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  } else if (diverging) {
    p <- p + scale_fill_gradient2(
      low = "#17335C", mid = "#F7F4ED", high = "#8F1D1D",
      midpoint = 0, limits = limits, oob = scales::squish, na.value = "grey93"
    )
  } else {
    cols <- palette %||% grDevices::colorRampPalette(c("#101D3A", "#174A7C", "#1E7A9C", "#36A793", "#F0D343", "#C94E27", "#7F1D1D"), space = "Lab")(256)
    p <- p + scale_fill_gradientn(colours = cols, limits = limits, oob = scales::squish, na.value = "grey93")
  }
  ggsave(file, p, width = 5.2, height = 4.7, dpi = 320, bg = "white")
  p
}

compute_line_maps <- function(cube, wave, mask, stat = NULL) {
  lambda0 <- line_rest * (1 + redshift)
  vel <- 299792.458 * (wave / lambda0 - 1)
  line_idx <- which(abs(vel) <= line_window_kms)
  cont_idx <- which(abs(vel) > continuum_inner_kms & abs(vel) <= continuum_outer_kms)
  if (length(line_idx) < 5L || length(cont_idx) < 5L) {
    stop("Not enough wavelength channels around ", line_label, ".", call. = FALSE)
  }

  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  flux <- flux_err <- snr <- velocity <- sigma <- asymmetry <- h3_proxy <- h4_proxy <- matrix(NA_real_, nx, ny)

  for (idx in which(mask)) {
    ij <- arrayInd(idx, .dim = c(nx, ny))
    spec <- cube[ij[1], ij[2], ]
    cont <- stats::median(spec[cont_idx], na.rm = TRUE)
    if (!is.finite(cont)) next
    prof <- spec[line_idx] - cont
    prof[!is.finite(prof)] <- 0
    v <- vel[line_idx]
    search <- abs(v) <= peak_search_kms
    peak <- which(search)[which.max(prof[search])]
    if (!length(peak) || !is.finite(prof[peak]) || prof[peak] <= 0) next
    local <- abs(v - v[peak]) <= centroid_window_kms
    pos <- pmax(prof, 0)
    pos[!local] <- 0
    sum_pos <- sum(pos)
    if (!is.finite(sum_pos) || sum_pos <= 0) next
    mu <- sum(v * pos) / sum_pos
    sig <- sqrt(sum(pos * (v - mu)^2) / sum_pos)
    flux[ij[1], ij[2]] <- sum_pos
    if (!is.null(stat)) {
      var_prof <- stat[ij[1], ij[2], line_idx]
      ferr <- sqrt(sum(pmax(var_prof[local], 0), na.rm = TRUE))
      flux_err[ij[1], ij[2]] <- ferr
      snr[ij[1], ij[2]] <- sum_pos / ferr
    }
    velocity[ij[1], ij[2]] <- mu
    sigma[ij[1], ij[2]] <- sig
    asymmetry[ij[1], ij[2]] <- (sum(pos[v > 0]) - sum(pos[v < 0])) / (sum_pos + .Machine$double.eps)
    if (is.finite(sig) && sig > 0) {
      z <- (v - mu) / sig
      h3_proxy[ij[1], ij[2]] <- sum(pos * z^3) / sum_pos
      h4_proxy[ij[1], ij[2]] <- sum(pos * z^4) / sum_pos - 3
    }
  }
  valid <- mask & is.finite(flux) & flux > 0 & is.finite(velocity) & is.finite(sigma)
  list(lambda0 = lambda0, vel_grid = vel, line_idx = line_idx, flux = flux, flux_err = flux_err, snr = snr, velocity = velocity, sigma = sigma, asymmetry = asymmetry, h3_proxy = h3_proxy, h4_proxy = h4_proxy, valid = valid)
}

build_line_support <- function(kin, border_width = 3L, flux_sigma = 2, snr_min = 3,
                               flux_quantile = 0.25, dilate = 0L) {
  nx <- nrow(kin$flux)
  ny <- ncol(kin$flux)
  border <- border_mask(nx, ny, border_width)
  b <- kin$flux[border & is.finite(kin$flux)]
  border_thr <- stats::median(b, na.rm = TRUE) + flux_sigma * stats::mad(b, constant = 1.4826, na.rm = TRUE)
  quant_thr <- stats::quantile(kin$flux[is.finite(kin$flux)], flux_quantile, na.rm = TRUE)
  flux_thr <- max(border_thr, quant_thr, na.rm = TRUE)
  support <- is.finite(kin$flux) & kin$flux > flux_thr
  if (any(is.finite(kin$snr))) {
    support <- support & is.finite(kin$snr) & kin$snr >= snr_min
  }
  if (sum(support, na.rm = TRUE) < 25L) {
    flux_thr <- stats::quantile(kin$flux[is.finite(kin$flux)], 0.50, na.rm = TRUE)
    support <- is.finite(kin$flux) & kin$flux > flux_thr
    if (any(is.finite(kin$snr))) {
      support <- support & is.finite(kin$snr) & kin$snr >= snr_min
    }
  }
  support <- better_galaxy_mask(support, close_iterations = 1L, fill_holes = TRUE, preserve_input = TRUE)
  support <- dilate_mask(support, dilate)
  list(mask = support, flux_threshold = flux_thr, border_threshold = border_thr,
       quantile_threshold = quant_thr, flux_quantile = flux_quantile, snr_min = snr_min)
}

baseline_subtract <- function(v, y) {
  edge <- c(seq_len(min(3L, length(y))), seq.int(max(1L, length(y) - 2L), length(y)))
  y - stats::median(y[edge], na.rm = TRUE)
}

build_path_features <- function(cube, wave, observed_wave, mask) {
  if (!requireNamespace("spectropath", quietly = TRUE)) {
    stop("spectropath is needed for path-signature segmentation.", call. = FALSE)
  }
  features <- c("p2", "p3u", "p3F", "p4F", "p4T", "p_pm")
  vel <- 299792.458 * (wave / observed_wave - 1)
  keep <- abs(vel) <= path_window_kms
  vel <- vel[keep]
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  feature_cube <- array(NA_real_, dim = c(nx, ny, length(features)), dimnames = list(NULL, NULL, features))
  rows <- list()
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      if (!isTRUE(mask[i, j])) next
      prof <- as.numeric(cube[i, j, keep])
      if (sum(is.finite(prof)) < 8L) next
      prof[!is.finite(prof)] <- 0
      prof <- baseline_subtract(vel, prof)
      amp <- max(abs(prof), na.rm = TRUE)
      if (!is.finite(amp) || amp <= 0) next
      pf <- tryCatch(spectropath::path_features(cbind(vel, prof / amp), depth = 4, normalize = TRUE, notation = "paper"), error = function(e) NULL)
      if (is.null(pf)) next
      vals <- as.numeric(pf[1, features, drop = TRUE])
      feature_cube[i, j, ] <- vals
      rows[[length(rows) + 1L]] <- data.frame(x = i, y = j, t(vals), check.names = FALSE)
    }
  }
  table <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(feature_cube = feature_cube, table = table, features = features)
}

message("Reading cube...")
cube_obj <- read_cube(cube_path)
cube <- cube_obj$cube
wave <- cube_obj$wave
message(sprintf("Observed %s = %.2f A", line_label, line_rest * (1 + redshift)))

message("Building comparison SED/white-light mask...")
mask_obj <- build_sed_mask(cube, mask_border_width, mask_sigma, mask_dilate)

message("Computing OIII maps over the full field...")
full_field <- matrix(TRUE, dim(cube)[1], dim(cube)[2])
kin <- compute_line_maps(cube, wave, full_field, stat = cube_obj$stat)

message("Building science support mask...")
line_support <- build_line_support(kin, mask_border_width, line_flux_sigma, line_snr_min, line_flux_quantile, mask_dilate)
support <- if (identical(support_mode, "sed")) mask_obj$mask else line_support$mask
kin$valid <- support & is.finite(kin$flux) & kin$flux > 0 & is.finite(kin$velocity) & is.finite(kin$sigma)
kin$velocity[kin$valid] <- kin$velocity[kin$valid] - stats::median(kin$velocity[kin$valid], na.rm = TRUE)

message("Running traditional kinematic-aware Capivara segmentation...")
kin_cube <- array(NA_real_, dim = c(dim(cube)[1], dim(cube)[2], 6L))
kin_cube[, , 1] <- log10(pmax(kin$flux, 0) + 1e-6)
kin_cube[, , 2] <- kin$velocity
kin_cube[, , 3] <- kin$sigma
kin_cube[, , 4] <- kin$asymmetry
kin_cube[, , 5] <- kin$h3_proxy
kin_cube[, , 6] <- kin$h4_proxy
kin_seg <- capivara::segment_large(
  list(imDat = kin_cube, hdr = cube_obj$header, axDat = NULL),
  Ncomp = n_segments,
  knn_k = knn_k,
  auto_k = FALSE,
  max_k = knn_k,
  feature_scale = "robust_col",
  spatial_weight = 0.15,
  mask = kin$valid,
  valid_mode = "finite",
  verbose = TRUE
)

path_seg <- NULL
path_features <- NULL
if (isTRUE(run_path_signatures)) {
  message("Running path-signature Capivara segmentation...")
  path_features <- build_path_features(cube, wave, kin$lambda0, kin$valid)
  path_seg <- capivara::segment_large(
    list(imDat = path_features$feature_cube, hdr = cube_obj$header, axDat = NULL),
    Ncomp = n_path_segments,
    knn_k = knn_k,
    auto_k = FALSE,
    max_k = knn_k,
    feature_scale = "robust_col",
    spatial_weight = 0.10,
    mask = kin$valid,
    valid_mode = "finite",
    verbose = TRUE
  )
}

message("Saving maps...")
p_white <- save_map(log10(mask_obj$white + 1), file.path(output_dir, paste0(prefix, "_sed_white_light.png")), "SED sky-subtracted white light")
p_sedmask <- save_map(ifelse(mask_obj$mask, 1, NA_real_), file.path(output_dir, paste0(prefix, "_sed_mask_comparison.png")), sprintf("SED mask comparison: %d spaxels", sum(mask_obj$mask)), discrete = TRUE)
p_mask <- save_map(ifelse(support, 1, NA_real_), file.path(output_dir, paste0(prefix, "_science_support_mask.png")), sprintf("%s support: %d spaxels", support_mode, sum(support)), discrete = TRUE)
p_flux <- save_map(log10(pmax(kin$flux, 0) + 1e-6), file.path(output_dir, paste0(prefix, "_flux_log.png")), paste(line_label, "log flux"), limits = robust_limits(log10(pmax(kin$flux, 0) + 1e-6)))
p_snr <- save_map(kin$snr, file.path(output_dir, paste0(prefix, "_snr.png")), paste(line_label, "line S/N"), limits = robust_limits(kin$snr))
p_vel <- save_map(kin$velocity, file.path(output_dir, paste0(prefix, "_velocity.png")), paste(line_label, "velocity, centered"), diverging = TRUE, limits = robust_limits(kin$velocity, symmetric = TRUE))
p_sig <- save_map(kin$sigma, file.path(output_dir, paste0(prefix, "_sigma.png")), paste(line_label, "sigma"), limits = robust_limits(kin$sigma))
p_asym <- save_map(kin$asymmetry, file.path(output_dir, paste0(prefix, "_asymmetry.png")), paste(line_label, "red-blue asymmetry"), diverging = TRUE, limits = c(-1, 1))
p_kinseg <- save_map(kin_seg$cluster_map, file.path(output_dir, sprintf("%s_kinematic_segments_n%d.png", prefix, n_segments)), "Kinematic-aware segments", discrete = TRUE)
p_pathseg <- NULL
if (!is.null(path_seg)) {
  p_pathseg <- save_map(path_seg$cluster_map, file.path(output_dir, sprintf("%s_path_signature_segments_n%d.png", prefix, n_path_segments)), "Path-signature segments", discrete = TRUE)
}

disc_result <- NULL
disc_plots <- list()
if (isTRUE(run_disc_model) && sum(kin$valid, na.rm = TRUE) >= 40L) {
  message("Trying exploratory rotating-disc model...")
  spaxels <- expand.grid(x = seq_len(dim(cube)[1]), y = seq_len(dim(cube)[2]))
  idx <- cbind(spaxels$x, spaxels$y)
  spaxels$velocity <- kin$velocity[idx]
  spaxels$flux <- kin$flux[idx]
  spaxels$valid <- kin$valid[idx]
  spaxels$seg_class <- "disc"
  geometry <- estimate_disc_geometry(
    spaxels,
    geometry = list(pa_deg = disc_pa_deg, inc_deg = disc_inc_deg, coordinate_convention = "nirvana"),
    allow_placeholder_inclination = TRUE,
    placeholder_inc_deg = disc_inc_deg
  )
  spaxels <- cbind(spaxels, deproject_coordinates(spaxels$x, spaxels$y, geometry))
  fit <- tryCatch(
    fit_axisymmetric_piecewise_model(spaxels, geometry, n_rings = 8, smooth_lambda = 5, fixed_vsys = geometry$vsys, robust = TRUE),
    error = function(e) e
  )
  if (!inherits(fit, "error")) {
    disc_result <- list(geometry = geometry, fit = fit)
    model_mat <- resid_mat <- matrix(NA_real_, dim(cube)[1], dim(cube)[2])
    model_mat[idx] <- fit$spaxels$v_axisym_model
    resid_mat[idx] <- fit$spaxels$v_axisym_resid
    disc_plots$model <- save_map(model_mat, file.path(output_dir, paste0(prefix, "_axisym_disc_model.png")), "Exploratory rotating-disc model", diverging = TRUE, limits = robust_limits(c(kin$velocity, model_mat), symmetric = TRUE))
    disc_plots$resid <- save_map(resid_mat, file.path(output_dir, paste0(prefix, "_axisym_disc_residual.png")), "OIII velocity - disc model", diverging = TRUE, limits = robust_limits(resid_mat, symmetric = TRUE))
    prof <- fit$profile
    curve <- ggplot(prof, aes(R, Vt)) +
      geom_line(linewidth = 0.8, colour = "#17335C") +
      geom_point(size = 1.8, colour = "#C85B3C") +
      theme_classic(base_size = 11) +
      labs(title = "Exploratory OIII rotation curve", x = "Radius (spaxels)", y = "Vt (km/s)")
    ggsave(file.path(output_dir, paste0(prefix, "_rotation_curve.png")), curve, width = 5.8, height = 4.2, dpi = 320, bg = "white")
    utils::write.csv(prof, file.path(output_dir, paste0(prefix, "_rotation_curve.csv")), row.names = FALSE)
  } else {
    message("Disc model skipped: ", conditionMessage(fit))
  }
}

message("Saving tables and panel...")
tab <- expand.grid(x = seq_len(dim(cube)[1]), y = seq_len(dim(cube)[2]))
idx <- cbind(tab$x, tab$y)
tab$sed_mask <- mask_obj$mask[idx]
tab$science_support <- support[idx]
tab$oiii_valid <- kin$valid[idx]
tab$oiii_flux <- kin$flux[idx]
tab$oiii_snr <- kin$snr[idx]
tab$oiii_velocity <- kin$velocity[idx]
tab$oiii_sigma <- kin$sigma[idx]
tab$oiii_asymmetry <- kin$asymmetry[idx]
tab$oiii_h3_proxy <- kin$h3_proxy[idx]
tab$oiii_h4_proxy <- kin$h4_proxy[idx]
tab$kinematic_segment <- kin_seg$cluster_map[idx]
if (!is.null(path_seg)) tab$path_signature_segment <- path_seg$cluster_map[idx]
utils::write.csv(tab, file.path(output_dir, paste0(prefix, "_spaxel_table.csv")), row.names = FALSE)
if (!is.null(path_features) && nrow(path_features$table)) {
  utils::write.csv(path_features$table, file.path(output_dir, paste0(prefix, "_path_signature_features.csv")), row.names = FALSE)
}

panel_items <- list(p_white, p_sedmask, p_mask, p_flux, p_snr, p_vel, p_sig, p_kinseg)
if (!is.null(p_pathseg)) panel_items <- c(panel_items, list(p_pathseg))
panel <- do.call(gridExtra::arrangeGrob, c(panel_items, ncol = 3))
ggsave(file.path(output_dir, paste0(prefix, "_summary_panel.png")), panel, width = 11.5, height = 8.0, dpi = 320, bg = "white")

result <- list(
  object_id = object_id,
  cube_path = cube_path,
  redshift = redshift,
  line_rest = line_rest,
  observed_wave = kin$lambda0,
  support = support,
  support_mode = support_mode,
  line_support = line_support,
  sed_mask = mask_obj$mask,
  white_light = mask_obj$white,
  kinematics = kin,
  kinematic_segmentation = kin_seg,
  path_features = path_features,
  path_segmentation = path_seg,
  disc_model = disc_result
)
saveRDS(result, file.path(output_dir, paste0(prefix, "_capivara_oiii_results.rds")))

message("Done. Outputs saved in: ", output_dir)
