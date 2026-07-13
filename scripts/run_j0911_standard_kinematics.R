#!/usr/bin/env Rscript

# Minimal standard kinematics for J0911.
# Edit this block, then Source in RStudio or run:
# Rscript run_j0911_standard_kinematics.R
#
# Orientation: maps are plotted in native cube coordinates only:
# x = FITS/R array dimension 1, y = FITS/R array dimension 2.
# No transpose, rot90, y reversal, or image-style flip is applied.

repo_root <- "/Users/rd23aag/Documents/GitHub/capivara"
cube_path <- "/Users/rd23aag/Documents/GitHub/HUB_2026/Emma/paper/J0911_30x30_cube.fits"
output_dir <- "/Users/rd23aag/Documents/GitHub/HUB_2026/Emma/paper"

object_id <- "J0911"
redshift <- 0.2622
lines <- list(
  list(slug = "oiii", label = "[OIII] 5007+4959", rest = c(5006.84, 4958.91)),
  list(slug = "halpha", label = "Halpha", rest = 6562.80)
)

n_segments <- 16
n_path_segments <- 18
knn_k <- 50
support_flux_quantile <- 0.70
support_snr_min <- 15
run_support_disk_model <- TRUE
disk_inc_deg <- 60
disk_n_rings <- 6
disk_smooth_lambda <- 5
segmentation_laplacian_steps <- 2
segmentation_laplacian_alpha <- 0.25
outflow_candidate_quantile <- 0.55
outflow_candidate_flux_quantile <- 0.35
outflow_candidate_snr_quantile <- 0.35
outflow_candidate_wing_quantile <- 0.50
outflow_candidate_gaussian_excess_quantile <- 0.50
outflow_candidate_disk_residual_quantile <- 0.50
outflow_disk_residual_weight <- 1.00
outflow_disk_residual_probability_weight <- 0.20
outflow_min_candidate_pixels <- 20
outflow_probability_threshold <- 0.60
outflow_probability_smoothing_steps <- 4

line_window_kms <- 700
peak_search_kms <- 350
centroid_window_kms <- 220
continuum_inner_kms <- 900
continuum_outer_kms <- 2200
profile_grid_kms <- seq(-500, 500, by = 25)

# Nothing below should need editing.

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(gridExtra)
})
pkgload::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "galaxy_mask.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "geometry.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "disc_model.R"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

hdr_get <- function(hdr, key, default = NA_character_) {
  i <- which(hdr == key)[1]
  if (!length(i) || is.na(i) || i >= length(hdr)) default else hdr[[i + 1L]]
}

read_cube <- function(path) {
  fits <- FITSio::readFITS(path, hdu = 1L, maxLines = 10000)
  stat <- tryCatch(FITSio::readFITS(path, hdu = 2L, maxLines = 10000)$imDat, error = function(e) NULL)
  cube <- fits$imDat
  crval <- as.numeric(hdr_get(fits$hdr, "CRVAL3"))
  crpix <- as.numeric(hdr_get(fits$hdr, "CRPIX3", "1"))
  cdelt <- as.numeric(hdr_get(fits$hdr, "CD3_3", hdr_get(fits$hdr, "CDELT3", "1")))
  wave <- crval + (seq_len(dim(cube)[3]) - crpix) * cdelt
  list(cube = cube, stat = stat, wave = wave, header = fits$hdr)
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

positive_robust_z_map <- function(x, mask) {
  vals <- x[mask & is.finite(x)]
  out <- matrix(NA_real_, nrow(x), ncol(x))
  if (length(vals) < 5L) return(out)
  center <- stats::median(vals, na.rm = TRUE)
  scale <- stats::mad(vals, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- stats::IQR(vals, na.rm = TRUE) / 1.349
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(vals, na.rm = TRUE)
  out[] <- pmax((x - center) / (scale + .Machine$double.eps), 0)
  out[!mask] <- NA_real_
  out
}

smooth_map_laplacian <- function(mat, mask, steps = segmentation_laplacian_steps, alpha = segmentation_laplacian_alpha) {
  out <- mat
  if (!is.finite(alpha) || alpha <= 0 || steps <= 0L) return(out)
  for (step in seq_len(as.integer(steps))) {
    old <- out
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        if (!isTRUE(mask[i, j]) || !is.finite(old[i, j])) next
        vals <- numeric()
        if (i > 1L && isTRUE(mask[i - 1L, j]) && is.finite(old[i - 1L, j])) vals <- c(vals, old[i - 1L, j])
        if (i < nrow(mat) && isTRUE(mask[i + 1L, j]) && is.finite(old[i + 1L, j])) vals <- c(vals, old[i + 1L, j])
        if (j > 1L && isTRUE(mask[i, j - 1L]) && is.finite(old[i, j - 1L])) vals <- c(vals, old[i, j - 1L])
        if (j < ncol(mat) && isTRUE(mask[i, j + 1L]) && is.finite(old[i, j + 1L])) vals <- c(vals, old[i, j + 1L])
        if (length(vals)) out[i, j] <- (1 - alpha) * old[i, j] + alpha * mean(vals)
      }
    }
  }
  out
}

smooth_feature_cube <- function(x, mask) {
  y <- x
  for (k in seq_len(dim(x)[3])) {
    y[, , k] <- smooth_map_laplacian(x[, , k], mask)
  }
  y
}

erode_mask <- function(mask, iterations = 1L) {
  out <- matrix(mask %in% TRUE, nrow = nrow(mask), ncol = ncol(mask))
  for (step in seq_len(as.integer(iterations))) {
    old <- out
    for (i in seq_len(nrow(mask))) {
      for (j in seq_len(ncol(mask))) {
        if (!old[i, j]) {
          out[i, j] <- FALSE
          next
        }
        neighbours <- c(
          i > 1L && old[i - 1L, j],
          i < nrow(mask) && old[i + 1L, j],
          j > 1L && old[i, j - 1L],
          j < ncol(mask) && old[i, j + 1L]
        )
        out[i, j] <- all(neighbours)
      }
    }
  }
  out
}

rank_probability_map <- function(x, mask) {
  out <- matrix(NA_real_, nrow(x), ncol(x))
  vals <- x[mask & is.finite(x)]
  if (length(vals) < 2L) return(out)
  ranks <- rank(vals, ties.method = "average")
  out[mask & is.finite(x)] <- (ranks - 1) / max(1, length(vals) - 1)
  out
}

save_map <- function(mat, file, title, diverging = FALSE, discrete = FALSE, limits = NULL) {
  df <- matrix_df(mat)
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    theme_void(base_size = 11) +
    labs(title = title, fill = NULL) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (discrete) {
    df$value <- factor(df$value)
    n <- max(3L, length(unique(stats::na.omit(df$value))))
    p <- ggplot(df, aes(x, y, fill = value)) +
      geom_raster() +
      coord_fixed(expand = FALSE) +
      scale_fill_manual(
        values = grDevices::colorRampPalette(c("#071A3F", "#004C99", "#00A6A6", "#F6D743", "#F07C24", "#B11226", "#5A0B23"), space = "Lab")(n),
        na.value = "grey93"
      ) +
      theme_void(base_size = 11) +
      labs(title = title, fill = "segment") +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
  } else if (diverging) {
    p <- p + scale_fill_gradient2(
      low = "#062B70", mid = "#FFF8E7", high = "#B11226",
      midpoint = 0, limits = limits, oob = scales::squish, na.value = "grey93"
    )
  } else {
    p <- p + scale_fill_gradientn(
      colours = grDevices::colorRampPalette(c("#050B2E", "#003C7A", "#007E9E", "#00A878", "#F4D03F", "#F97316", "#B11226", "#4A061D"), space = "Lab")(256),
      limits = limits, oob = scales::squish, na.value = "grey93"
    )
  }
  ggsave(file, p, width = 5.2, height = 4.7, dpi = 320, bg = "white")
  p
}

component_indices <- function(wave, rest_wave) {
  lambda0 <- rest_wave * (1 + redshift)
  vel_grid <- 299792.458 * (wave / lambda0 - 1)
  line_idx <- which(abs(vel_grid) <= line_window_kms)
  cont_idx <- which(abs(vel_grid) > continuum_inner_kms & abs(vel_grid) <= continuum_outer_kms)
  if (length(line_idx) < 5L || length(cont_idx) < 5L) {
    stop("Insufficient wavelength coverage near ", round(lambda0, 2), " A", call. = FALSE)
  }
  list(lambda0 = lambda0, vel_grid = vel_grid, line_idx = line_idx, cont_idx = cont_idx)
}

compute_line_maps <- function(cube, stat, wave, line) {
  parts <- lapply(line$rest, function(rest_wave) component_indices(wave, rest_wave))
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  flux <- flux_err <- snr <- velocity <- sigma <- asymmetry <- wing_fraction <- path_score <- matrix(NA_real_, nx, ny)
  blue_wing <- red_wing <- signed_blue_excess <- signed_red_excess <- signed_path_response <- wing_asymmetry <- path_skewness <- path_kurtosis <- h3 <- h4 <- peak_offset <- gaussian_wing_excess <- outflow_score <- matrix(NA_real_, nx, ny)
  path_profiles <- array(NA_real_, dim = c(nx, ny, length(profile_grid_kms)))

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      spec <- cube[i, j, ]
      prof_sum <- rep(0, length(profile_grid_kms))
      var_sum <- 0
      for (part in parts) {
        cont <- stats::median(spec[part$cont_idx], na.rm = TRUE)
        if (!is.finite(cont)) next
        prof <- spec[part$line_idx] - cont
        prof[!is.finite(prof)] <- 0
        pos <- pmax(prof, 0)
        v <- part$vel_grid[part$line_idx]
        prof_sum <- prof_sum + approx(v, pos, xout = profile_grid_kms, yleft = 0, yright = 0, ties = "ordered")$y
        if (!is.null(stat)) {
          var_prof <- stat[i, j, part$line_idx]
          var_sum <- var_sum + sum(pmax(var_prof[abs(v) <= centroid_window_kms], 0), na.rm = TRUE)
        }
      }

      prof_sum[!is.finite(prof_sum)] <- 0
      search <- abs(profile_grid_kms) <= peak_search_kms
      peak <- which(search)[which.max(prof_sum[search])]
      if (!length(peak) || !is.finite(prof_sum[peak]) || prof_sum[peak] <= 0) next

      local <- abs(profile_grid_kms - profile_grid_kms[peak]) <= centroid_window_kms
      pos <- prof_sum
      pos[!local] <- 0
      sum_pos <- sum(pos)
      if (!is.finite(sum_pos) || sum_pos <= 0) next

      mu <- sum(profile_grid_kms * pos) / sum_pos
      sig <- sqrt(sum(pos * (profile_grid_kms - mu)^2) / sum_pos)
      norm_prof <- prof_sum / (sum(prof_sum) + .Machine$double.eps)
      centered <- profile_grid_kms - mu
      wing <- abs(centered) > max(80, 1.5 * sig) & abs(centered) <= max(profile_grid_kms)
      blue <- centered < -max(80, 1.5 * sig)
      red <- centered > max(80, 1.5 * sig)
      wing_frac <- sum(norm_prof[wing], na.rm = TRUE)
      blue_frac <- sum(norm_prof[blue], na.rm = TRUE)
      red_frac <- sum(norm_prof[red], na.rm = TRUE)
      wing_asym <- (red_frac - blue_frac) / (red_frac + blue_frac + .Machine$double.eps)
      centered_scale <- centered / (sig + .Machine$double.eps)
      H3 <- (2 * sqrt(2) * centered_scale^3 - 3 * sqrt(2) * centered_scale) / sqrt(6)
      H4 <- (4 * centered_scale^4 - 12 * centered_scale^2 + 3) / sqrt(24)
      gaussian_core <- stats::dnorm(profile_grid_kms, mean = mu, sd = sig)
      gaussian_core <- gaussian_core / (sum(gaussian_core) + .Machine$double.eps)
      gaussian_residual <- norm_prof - gaussian_core
      profile_excess <- pmax(norm_prof - gaussian_core, 0)
      gaussian_excess <- sum(profile_excess[wing], na.rm = TRUE)
      blue_excess_signed <- sum(gaussian_residual[blue], na.rm = TRUE)
      red_excess_signed <- sum(gaussian_residual[red], na.rm = TRUE)

      flux[i, j] <- sum_pos
      velocity[i, j] <- mu
      sigma[i, j] <- sig
      asymmetry[i, j] <- (sum(pos[profile_grid_kms > mu]) - sum(pos[profile_grid_kms < mu])) / (sum_pos + .Machine$double.eps)
      wing_fraction[i, j] <- wing_frac
      blue_wing[i, j] <- blue_frac
      red_wing[i, j] <- red_frac
      signed_blue_excess[i, j] <- blue_excess_signed
      signed_red_excess[i, j] <- red_excess_signed
      signed_path_response[i, j] <- red_excess_signed - blue_excess_signed
      wing_asymmetry[i, j] <- wing_asym
      path_skewness[i, j] <- sum(norm_prof * centered_scale^3, na.rm = TRUE)
      path_kurtosis[i, j] <- sum(norm_prof * centered_scale^4, na.rm = TRUE) - 3
      h3[i, j] <- sum(gaussian_residual * H3, na.rm = TRUE) / (sum(gaussian_core * H3^2, na.rm = TRUE) + .Machine$double.eps)
      h4[i, j] <- sum(gaussian_residual * H4, na.rm = TRUE) / (sum(gaussian_core * H4^2, na.rm = TRUE) + .Machine$double.eps)
      peak_offset[i, j] <- profile_grid_kms[peak] - mu
      gaussian_wing_excess[i, j] <- gaussian_excess
      path_profiles[i, j, ] <- norm_prof
      if (!is.null(stat)) {
        ferr <- sqrt(var_sum)
        flux_err[i, j] <- ferr
        snr[i, j] <- sum_pos / ferr
      }
      path_score[i, j] <- wing_frac * log1p(ifelse(is.finite(snr[i, j]), snr[i, j], sum_pos))
    }
  }

  flux_thr <- stats::quantile(flux[is.finite(flux)], support_flux_quantile, na.rm = TRUE)
  support <- is.finite(flux) & flux > flux_thr
  if (any(is.finite(snr))) {
    support <- support & is.finite(snr) & snr >= support_snr_min
  }
  support <- better_galaxy_mask(support, close_iterations = 1L, fill_holes = TRUE, preserve_input = FALSE)
  valid <- support & is.finite(flux) & flux > 0 & is.finite(velocity) & is.finite(sigma)
  velocity[valid] <- velocity[valid] - stats::median(velocity[valid], na.rm = TRUE)

  positive_robust_z <- function(x) {
    vals <- x[valid & is.finite(x)]
    if (length(vals) < 5L) return(matrix(NA_real_, nrow(x), ncol(x)))
    center <- stats::median(vals, na.rm = TRUE)
    scale <- stats::mad(vals, center = center, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) {
      scale <- stats::IQR(vals, na.rm = TRUE) / 1.349
    }
    if (!is.finite(scale) || scale <= 0) scale <- stats::sd(vals, na.rm = TRUE)
    z <- (x - center) / (scale + .Machine$double.eps)
    pmax(z, 0)
  }

  outflow_score <- positive_robust_z(wing_fraction) +
    1.40 * positive_robust_z(gaussian_wing_excess) +
    0.80 * positive_robust_z(abs(wing_asymmetry)) +
    0.65 * positive_robust_z(abs(path_skewness)) +
    0.60 * positive_robust_z(path_kurtosis) +
    0.55 * positive_robust_z(sigma)
  outflow_score[!valid] <- NA_real_
  outflow_score_smooth <- smooth_map_laplacian(outflow_score, valid)

  list(
    lambda0 = vapply(parts, `[[`, numeric(1), "lambda0"),
    flux = flux,
    flux_err = flux_err,
    snr = snr,
    velocity = velocity,
    sigma = sigma,
    asymmetry = asymmetry,
    wing_fraction = wing_fraction,
    blue_wing = blue_wing,
    red_wing = red_wing,
    signed_blue_excess = signed_blue_excess,
    signed_red_excess = signed_red_excess,
    signed_path_response = signed_path_response,
    wing_asymmetry = wing_asymmetry,
    path_skewness = path_skewness,
    path_kurtosis = path_kurtosis,
    h3 = h3,
    h4 = h4,
    peak_offset = peak_offset,
    gaussian_wing_excess = gaussian_wing_excess,
    path_score = path_score,
    outflow_score = outflow_score,
    outflow_score_smooth = outflow_score_smooth,
    path_profiles = path_profiles,
    support = support,
    valid = valid
  )
}

fit_support_disk_model <- function(kin, dims, slug, line_label) {
  if (!isTRUE(run_support_disk_model) || sum(kin$valid, na.rm = TRUE) < 40L) {
    return(NULL)
  }

  spaxels <- expand.grid(x = seq_len(dims[1]), y = seq_len(dims[2]))
  idx <- cbind(spaxels$x, spaxels$y)
  spaxels$velocity <- kin$velocity[idx]
  spaxels$flux <- kin$flux[idx]
  spaxels$snr <- kin$snr[idx]
  spaxels$sigma <- kin$sigma[idx]
  spaxels$valid <- kin$valid[idx]
  spaxels$seg_class <- "disc"
  spaxels$velocity_error <- pmax(spaxels$sigma / pmax(spaxels$snr, 1), 5)

  geometry <- estimate_disc_geometry(
    spaxels,
    geometry = list(inc_deg = disk_inc_deg, coordinate_convention = "nirvana"),
    allow_placeholder_inclination = TRUE,
    placeholder_inc_deg = disk_inc_deg
  )
  spaxels <- cbind(spaxels, deproject_coordinates(spaxels$x, spaxels$y, geometry))
  fit <- tryCatch(
    fit_axisymmetric_piecewise_model(
      spaxels,
      geometry,
      n_rings = disk_n_rings,
      smooth_lambda = disk_smooth_lambda,
      fixed_vsys = geometry$vsys,
      robust = TRUE
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    message("Support-mask disk model skipped for ", line_label, ": ", conditionMessage(fit))
    return(NULL)
  }

  model_mat <- resid_mat <- matrix(NA_real_, dims[1], dims[2])
  model_mat[idx] <- fit$spaxels$v_axisym_model
  resid_mat[idx] <- fit$spaxels$v_axisym_resid
  velocity_map <- ifelse(kin$valid, kin$velocity, NA_real_)
  model_mat <- ifelse(kin$valid, model_mat, NA_real_)
  resid_mat <- ifelse(kin$valid, resid_mat, NA_real_)
  vel_lim <- robust_limits(c(velocity_map, model_mat), symmetric = TRUE)
  res_lim <- robust_limits(resid_mat, symmetric = TRUE)

  p_model <- save_map(model_mat, file.path(output_dir, paste0(slug, "_support_disk_model.png")), paste(line_label, "support-mask disk model"), diverging = TRUE, limits = vel_lim)
  p_resid <- save_map(resid_mat, file.path(output_dir, paste0(slug, "_support_disk_residual.png")), paste(line_label, "velocity - disk model"), diverging = TRUE, limits = res_lim)
  p_data <- save_map(velocity_map, file.path(output_dir, paste0(slug, "_support_disk_data_velocity.png")), paste(line_label, "data velocity"), diverging = TRUE, limits = vel_lim)

  prof <- fit$profile
  prof$Vt_signed <- prof$Vt
  prof$Vt_speed <- abs(prof$Vt)
  curve <- ggplot(prof, aes(R, Vt_speed)) +
    geom_line(linewidth = 0.9, colour = "#071A3F") +
    geom_point(size = 2.0, colour = "#B11226") +
    theme_classic(base_size = 11) +
    labs(title = paste(line_label, "support-mask rotation speed"), x = "Radius (spaxels)", y = "|Vt| (km/s)")
  ggsave(file.path(output_dir, paste0(slug, "_support_disk_rotation_curve.png")), curve, width = 5.8, height = 4.2, dpi = 320, bg = "white")
  utils::write.csv(prof, file.path(output_dir, paste0(slug, "_support_disk_rotation_curve.csv")), row.names = FALSE)

  panel <- gridExtra::arrangeGrob(p_data, p_model, p_resid, curve, ncol = 4)
  ggsave(file.path(output_dir, paste0(slug, "_support_disk_model_panel.png")), panel, width = 13.5, height = 3.7, dpi = 320, bg = "white")

  list(geometry = geometry, fit = fit, model = model_mat, residual = resid_mat, panel = panel)
}

build_spectropath_feature_cube <- function(kin, feature_names = c("p2", "p3u", "p3F", "p4F", "p4T", "p_pm")) {
  nx <- nrow(kin$valid)
  ny <- ncol(kin$valid)
  feature_cube <- array(NA_real_, dim = c(nx, ny, length(feature_names)), dimnames = list(NULL, NULL, feature_names))
  rows <- list()

  if (!requireNamespace("spectropath", quietly = TRUE)) {
    message("spectropath is not installed; signed native path maps will still be produced.")
    return(list(feature_cube = feature_cube, table = data.frame(), features = feature_names, available = FALSE))
  }

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      if (!isTRUE(kin$valid[i, j])) next
      profile <- as.numeric(kin$path_profiles[i, j, ])
      if (sum(is.finite(profile)) < 8L) next
      profile[!is.finite(profile)] <- 0
      if (sum(profile) <= 0) next
      profile <- profile / (sum(profile) + .Machine$double.eps)
      mu <- sum(profile_grid_kms * profile) / (sum(profile) + .Machine$double.eps)
      signed_path <- cbind(profile_grid_kms - mu, profile)
      pf <- tryCatch(
        spectropath::path_features(signed_path, depth = 4, normalize = TRUE, notation = "paper"),
        error = function(e) NULL
      )
      if (is.null(pf)) next
      present <- intersect(feature_names, colnames(pf))
      if (!length(present)) next
      vals <- rep(NA_real_, length(feature_names))
      names(vals) <- feature_names
      vals[present] <- as.numeric(pf[1, present, drop = TRUE])
      feature_cube[i, j, ] <- vals
      rows[[length(rows) + 1L]] <- data.frame(x = i, y = j, t(vals), check.names = FALSE)
    }
  }

  list(
    feature_cube = feature_cube,
    table = if (length(rows)) do.call(rbind, rows) else data.frame(),
    features = feature_names,
    available = TRUE
  )
}

probability_gradient_map <- function(probability, mask) {
  grad <- matrix(NA_real_, nrow(probability), ncol(probability))
  for (i in seq_len(nrow(probability))) {
    for (j in seq_len(ncol(probability))) {
      if (!isTRUE(mask[i, j]) || !is.finite(probability[i, j])) next
      vals <- numeric()
      if (i > 1L && isTRUE(mask[i - 1L, j]) && is.finite(probability[i - 1L, j])) vals <- c(vals, probability[i - 1L, j])
      if (i < nrow(probability) && isTRUE(mask[i + 1L, j]) && is.finite(probability[i + 1L, j])) vals <- c(vals, probability[i + 1L, j])
      if (j > 1L && isTRUE(mask[i, j - 1L]) && is.finite(probability[i, j - 1L])) vals <- c(vals, probability[i, j - 1L])
      if (j < ncol(probability) && isTRUE(mask[i, j + 1L]) && is.finite(probability[i, j + 1L])) vals <- c(vals, probability[i, j + 1L])
      if (length(vals)) grad[i, j] <- max(abs(probability[i, j] - vals), na.rm = TRUE)
    }
  }
  grad
}

select_atlas_spaxels <- function(score, mask, n = 9L, min_sep = 3) {
  idx <- which(mask & is.finite(score), arr.ind = TRUE)
  if (!nrow(idx)) return(idx)
  ord <- order(score[idx], decreasing = TRUE, na.last = NA)
  idx <- idx[ord, , drop = FALSE]
  keep <- matrix(numeric(0), ncol = 2)
  for (r in seq_len(nrow(idx))) {
    xy <- idx[r, , drop = FALSE]
    if (!nrow(keep)) {
      keep <- rbind(keep, xy)
    } else {
      dist <- sqrt((keep[, 1] - xy[1, 1])^2 + (keep[, 2] - xy[1, 2])^2)
      if (all(dist >= min_sep)) keep <- rbind(keep, xy)
    }
    if (nrow(keep) >= n) break
  }
  keep
}

save_profile_atlas <- function(kin, outflow_probability, candidate, slug, line_label) {
  probability_gradient <- probability_gradient_map(outflow_probability, kin$valid)
  atlas_score <- outflow_probability * probability_gradient
  atlas_mask <- kin$valid & is.finite(atlas_score) & outflow_probability >= stats::quantile(outflow_probability[kin$valid], 0.55, na.rm = TRUE)
  atlas_points <- select_atlas_spaxels(atlas_score, atlas_mask, n = 9L, min_sep = 3)
  if (nrow(atlas_points) < 4L) {
    atlas_points <- select_atlas_spaxels(outflow_probability, kin$valid, n = 9L, min_sep = 3)
  }
  if (!nrow(atlas_points)) return(NULL)

  atlas_rows <- list()
  summary_rows <- vector("list", nrow(atlas_points))
  for (k in seq_len(nrow(atlas_points))) {
    i <- atlas_points[k, 1]
    j <- atlas_points[k, 2]
    profile <- kin$path_profiles[i, j, ]
    profile[!is.finite(profile)] <- 0
    profile <- profile / (sum(profile) + .Machine$double.eps)
    mu <- sum(profile_grid_kms * profile) / (sum(profile) + .Machine$double.eps)
    sig <- sqrt(sum(profile * (profile_grid_kms - mu)^2) / (sum(profile) + .Machine$double.eps))
    gaussian <- stats::dnorm(profile_grid_kms, mean = mu, sd = sig)
    gaussian <- gaussian / (sum(gaussian) + .Machine$double.eps)
    label <- sprintf("x=%d y=%d | p=%.2f | h3=%.2f h4=%.2f", i, j, outflow_probability[i, j], kin$h3[i, j], kin$h4[i, j])
    atlas_rows[[length(atlas_rows) + 1L]] <- data.frame(
      spaxel = label,
      x = i,
      y = j,
      velocity = profile_grid_kms - mu,
      profile = profile,
      gaussian = gaussian,
      stringsAsFactors = FALSE
    )
    summary_rows[[k]] <- data.frame(
      x = i,
      y = j,
      outflow_probability = outflow_probability[i, j],
      probability_gradient = probability_gradient[i, j],
      h3 = kin$h3[i, j],
      h4 = kin$h4[i, j],
      gaussian_wing_excess = kin$gaussian_wing_excess[i, j],
      wing_fraction = kin$wing_fraction[i, j],
      candidate = isTRUE(candidate[i, j])
    )
  }
  atlas_df <- do.call(rbind, atlas_rows)
  summary_df <- do.call(rbind, summary_rows)
  utils::write.csv(summary_df, file.path(output_dir, paste0(slug, "_profile_atlas_spaxels.csv")), row.names = FALSE)

  atlas_long <- rbind(
    data.frame(spaxel = atlas_df$spaxel, velocity = atlas_df$velocity, flux = atlas_df$profile, curve = "profile"),
    data.frame(spaxel = atlas_df$spaxel, velocity = atlas_df$velocity, flux = atlas_df$gaussian, curve = "Gaussian reference")
  )
  atlas_plot <- ggplot(atlas_long, aes(velocity, flux, colour = curve, linewidth = curve)) +
    geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.25) +
    geom_line() +
    facet_wrap(~spaxel, ncol = 3, scales = "free_y") +
    scale_colour_manual(values = c("profile" = "#071A3F", "Gaussian reference" = "#D34E24")) +
    scale_linewidth_manual(values = c("profile" = 0.8, "Gaussian reference" = 0.55)) +
    coord_cartesian(xlim = c(-360, 360)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(size = 8)) +
    labs(
      title = paste(line_label, "profiles at high outflow-probability change"),
      subtitle = "Native cube orientation; profiles are recentered on their local centroid",
      x = "Velocity relative to local centroid (km/s)",
      y = "Normalized line profile",
      colour = NULL,
      linewidth = NULL
  )
  ggsave(file.path(output_dir, paste0(slug, "_profile_atlas_high_change.png")), atlas_plot, width = 10.5, height = 7.2, dpi = 320, bg = "white")

  footprint_rows <- list()
  for (k in seq_len(nrow(atlas_points))) {
    i <- atlas_points[k, 1]
    j <- atlas_points[k, 2]
    profile <- kin$path_profiles[i, j, ]
    profile[!is.finite(profile)] <- 0
    profile <- profile / (sum(profile) + .Machine$double.eps)
    mu <- sum(profile_grid_kms * profile) / (sum(profile) + .Machine$double.eps)
    sig <- sqrt(sum(profile * (profile_grid_kms - mu)^2) / (sum(profile) + .Machine$double.eps))
    gaussian <- stats::dnorm(profile_grid_kms, mean = mu, sd = sig)
    gaussian <- gaussian / (sum(gaussian) + .Machine$double.eps)
    velocity <- profile_grid_kms - mu
    local_max <- max(c(profile, gaussian), na.rm = TRUE)
    if (!is.finite(local_max) || local_max <= 0) next
    scale_x <- 1.45
    scale_y <- 1.65
    footprint_rows[[length(footprint_rows) + 1L]] <- data.frame(
      x = i + (velocity / 360) * scale_x,
      y = j + ((profile / local_max) - 0.5) * scale_y,
      curve = "profile",
      spaxel = sprintf("%d,%d", i, j)
    )
    footprint_rows[[length(footprint_rows) + 1L]] <- data.frame(
      x = i + (velocity / 360) * scale_x,
      y = j + ((gaussian / local_max) - 0.5) * scale_y,
      curve = "Gaussian reference",
      spaxel = sprintf("%d,%d", i, j)
    )
  }
  footprint_curves <- if (length(footprint_rows)) do.call(rbind, footprint_rows) else data.frame()
  footprint_df <- matrix_df(outflow_probability)
  point_df <- data.frame(
    x = atlas_points[, 1],
    y = atlas_points[, 2],
    label = seq_len(nrow(atlas_points))
  )
  footprint_plot <- ggplot(footprint_df, aes(x, y)) +
    geom_raster(aes(fill = value), alpha = 0.72) +
    geom_contour(
      data = matrix_df(ifelse(kin$valid, 1, NA_real_)),
      aes(x = x, y = y, z = value),
      breaks = 0.5,
      colour = "#071A3F",
      linewidth = 0.25,
      inherit.aes = FALSE
    ) +
    geom_path(
      data = footprint_curves,
      aes(x, y, colour = curve, group = interaction(spaxel, curve)),
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    geom_point(data = point_df, aes(x, y), inherit.aes = FALSE, size = 1.1, colour = "#071A3F") +
    geom_text(data = point_df, aes(x, y + 1.45, label = label), inherit.aes = FALSE, size = 2.5, colour = "#071A3F") +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours = grDevices::colorRampPalette(c("#050B2E", "#003C7A", "#007E9E", "#00A878", "#F4D03F", "#F97316", "#B11226", "#4A061D"), space = "Lab")(256),
      limits = c(0, 1),
      oob = scales::squish,
      na.value = "grey94"
    ) +
    scale_colour_manual(values = c("profile" = "#071A3F", "Gaussian reference" = "#D34E24")) +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(
      title = paste(line_label, "profile atlas on native cube footprint"),
      subtitle = "Mini-profiles are drawn at their selected spaxel positions; no rotation or y reversal",
      fill = "outflow p",
      colour = NULL
    )
  ggsave(file.path(output_dir, paste0(slug, "_profile_atlas_footprint.png")), footprint_plot, width = 8.0, height = 7.5, dpi = 320, bg = "white")

  list(plot = atlas_plot, footprint_plot = footprint_plot, spaxels = summary_df, gradient = probability_gradient)
}

run_one_line <- function(cube_obj, line) {
  slug <- paste0(tolower(object_id), "_", line$slug)
  obs <- paste(round(line$rest * (1 + redshift), 2), collapse = ", ")
  message(sprintf("Processing %s at observed %s A", line$label, obs))
  kin <- compute_line_maps(cube_obj$cube, cube_obj$stat, cube_obj$wave, line)
  spectropath_features <- build_spectropath_feature_cube(kin)

  feature_cube <- array(NA_real_, dim = c(dim(cube_obj$cube)[1], dim(cube_obj$cube)[2], 4L))
  feature_cube[, , 1] <- log10(pmax(kin$flux, 0) + 1e-6)
  feature_cube[, , 2] <- kin$velocity
  feature_cube[, , 3] <- kin$sigma
  feature_cube[, , 4] <- kin$asymmetry
  feature_cube <- smooth_feature_cube(feature_cube, kin$valid)
  seg <- capivara::segment_large(
    list(imDat = feature_cube, hdr = cube_obj$header, axDat = NULL),
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

  n_native_path_features <- length(profile_grid_kms) + 15L
  n_spectropath_features <- length(spectropath_features$features)
  path_feature_cube <- array(NA_real_, dim = c(dim(cube_obj$cube)[1], dim(cube_obj$cube)[2], n_native_path_features + n_spectropath_features))
  path_feature_cube[, , seq_along(profile_grid_kms)] <- kin$path_profiles
  path_feature_cube[, , length(profile_grid_kms) + 1L] <- kin$velocity
  path_feature_cube[, , length(profile_grid_kms) + 2L] <- kin$sigma
  path_feature_cube[, , length(profile_grid_kms) + 3L] <- kin$asymmetry
  path_feature_cube[, , length(profile_grid_kms) + 4L] <- kin$wing_fraction
  path_feature_cube[, , length(profile_grid_kms) + 5L] <- kin$blue_wing
  path_feature_cube[, , length(profile_grid_kms) + 6L] <- kin$red_wing
  path_feature_cube[, , length(profile_grid_kms) + 7L] <- kin$wing_asymmetry
  path_feature_cube[, , length(profile_grid_kms) + 8L] <- kin$path_skewness
  path_feature_cube[, , length(profile_grid_kms) + 9L] <- kin$path_kurtosis
  path_feature_cube[, , length(profile_grid_kms) + 10L] <- kin$gaussian_wing_excess
  path_feature_cube[, , length(profile_grid_kms) + 11L] <- kin$h3
  path_feature_cube[, , length(profile_grid_kms) + 12L] <- kin$h4
  path_feature_cube[, , length(profile_grid_kms) + 13L] <- kin$signed_blue_excess
  path_feature_cube[, , length(profile_grid_kms) + 14L] <- kin$signed_red_excess
  path_feature_cube[, , length(profile_grid_kms) + 15L] <- kin$signed_path_response
  if (n_spectropath_features) {
    path_feature_cube[, , (n_native_path_features + 1L):(n_native_path_features + n_spectropath_features)] <- spectropath_features$feature_cube
  }
  path_feature_cube <- smooth_feature_cube(path_feature_cube, kin$valid)
  path_seg <- capivara::segment_large(
    list(imDat = path_feature_cube, hdr = cube_obj$header, axDat = NULL),
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

  p_flux <- save_map(log10(pmax(kin$flux, 0) + 1e-6), file.path(output_dir, paste0(slug, "_flux_log.png")), paste(line$label, "log flux"), limits = robust_limits(log10(pmax(kin$flux, 0) + 1e-6)))
  p_snr <- save_map(kin$snr, file.path(output_dir, paste0(slug, "_snr.png")), paste(line$label, "line S/N"), limits = robust_limits(kin$snr))
  p_support <- save_map(ifelse(kin$support, 1, NA_real_), file.path(output_dir, paste0(slug, "_support_mask.png")), sprintf("%s support: %d spaxels", line$label, sum(kin$support)), discrete = TRUE)
  velocity_map <- ifelse(kin$valid, kin$velocity, NA_real_)
  sigma_map <- ifelse(kin$valid, kin$sigma, NA_real_)
  path_score_map <- ifelse(kin$valid, kin$path_score, NA_real_)
  outflow_score_map <- ifelse(kin$valid, kin$outflow_score_smooth, NA_real_)
  wing_asymmetry_map <- ifelse(kin$valid, kin$wing_asymmetry, NA_real_)
  red_wing_map <- ifelse(kin$valid, kin$red_wing, NA_real_)
  blue_wing_map <- ifelse(kin$valid, kin$blue_wing, NA_real_)
  signed_blue_map <- ifelse(kin$valid, kin$signed_blue_excess, NA_real_)
  signed_red_map <- ifelse(kin$valid, kin$signed_red_excess, NA_real_)
  signed_response_map <- ifelse(kin$valid, kin$signed_path_response, NA_real_)
  gaussian_excess_map <- ifelse(kin$valid, kin$gaussian_wing_excess, NA_real_)
  h3_map <- ifelse(kin$valid, kin$h3, NA_real_)
  h4_map <- ifelse(kin$valid, kin$h4, NA_real_)
  p_vel <- save_map(velocity_map, file.path(output_dir, paste0(slug, "_velocity.png")), paste(line$label, "velocity, centered"), diverging = TRUE, limits = robust_limits(velocity_map, symmetric = TRUE))
  p_sig <- save_map(sigma_map, file.path(output_dir, paste0(slug, "_sigma.png")), paste(line$label, "sigma"), limits = robust_limits(sigma_map))
  p_seg <- save_map(seg$cluster_map, file.path(output_dir, sprintf("%s_kinematic_segments_n%d.png", slug, n_segments)), paste(line$label, "kinematic-aware segments"), discrete = TRUE)
  p_path_score <- save_map(path_score_map, file.path(output_dir, paste0(slug, "_path_wing_score.png")), paste(line$label, "path wing score"), limits = robust_limits(path_score_map))
  p_path_seg <- save_map(path_seg$cluster_map, file.path(output_dir, sprintf("%s_path_segments_n%d.png", slug, n_path_segments)), paste(line$label, "path/profile segments"), discrete = TRUE)
  p_outflow_score <- save_map(outflow_score_map, file.path(output_dir, paste0(slug, "_path_outflow_score.png")), paste(line$label, "path outflow score"), limits = robust_limits(outflow_score_map, probs = c(0.02, 0.96)))
  p_wing_asym <- save_map(wing_asymmetry_map, file.path(output_dir, paste0(slug, "_path_wing_asymmetry.png")), paste(line$label, "red-blue wing balance"), diverging = TRUE, limits = robust_limits(wing_asymmetry_map, symmetric = TRUE))
  p_red_wing <- save_map(red_wing_map, file.path(output_dir, paste0(slug, "_path_red_wing.png")), paste(line$label, "red-wing coefficient"), limits = robust_limits(red_wing_map, probs = c(0.02, 0.96)))
  p_blue_wing <- save_map(blue_wing_map, file.path(output_dir, paste0(slug, "_path_blue_wing.png")), paste(line$label, "blue-wing coefficient"), limits = robust_limits(blue_wing_map, probs = c(0.02, 0.96)))
  p_signed_blue <- save_map(signed_blue_map, file.path(output_dir, paste0(slug, "_signed_blue_residual_path.png")), paste(line$label, "signed blue residual path"), diverging = TRUE, limits = robust_limits(signed_blue_map, symmetric = TRUE))
  p_signed_red <- save_map(signed_red_map, file.path(output_dir, paste0(slug, "_signed_red_residual_path.png")), paste(line$label, "signed red residual path"), diverging = TRUE, limits = robust_limits(signed_red_map, symmetric = TRUE))
  p_signed_response <- save_map(signed_response_map, file.path(output_dir, paste0(slug, "_signed_red_minus_blue_path.png")), paste(line$label, "signed red-blue path"), diverging = TRUE, limits = robust_limits(signed_response_map, symmetric = TRUE))
  p_gaussian_excess <- save_map(gaussian_excess_map, file.path(output_dir, paste0(slug, "_gaussian_wing_excess.png")), paste(line$label, "Gaussian wing excess"), limits = robust_limits(gaussian_excess_map, probs = c(0.02, 0.96)))
  p_h3 <- save_map(h3_map, file.path(output_dir, paste0(slug, "_h3_map.png")), paste(line$label, "h3 profile shape"), diverging = TRUE, limits = robust_limits(h3_map, symmetric = TRUE))
  p_h4 <- save_map(h4_map, file.path(output_dir, paste0(slug, "_h4_map.png")), paste(line$label, "h4 profile shape"), diverging = TRUE, limits = robust_limits(h4_map, symmetric = TRUE))
  spectropath_plots <- list()
  if (isTRUE(spectropath_features$available) && n_spectropath_features) {
    for (ff in spectropath_features$features) {
      fmap <- spectropath_features$feature_cube[, , ff]
      spectropath_plots[[ff]] <- save_map(
        fmap,
        file.path(output_dir, paste0(slug, "_spectropath_", ff, ".png")),
        paste(line$label, "SpectroPath", ff),
        diverging = TRUE,
        limits = robust_limits(fmap, symmetric = TRUE)
      )
    }
    if (nrow(spectropath_features$table)) {
      utils::write.csv(spectropath_features$table, file.path(output_dir, paste0(slug, "_spectropath_features.csv")), row.names = FALSE)
    }
    spectropath_panel <- gridExtra::arrangeGrob(grobs = spectropath_plots, ncol = 3)
    ggsave(file.path(output_dir, paste0(slug, "_spectropath_feature_panel.png")), spectropath_panel, width = 11.5, height = 4.9, dpi = 320, bg = "white")
  }
  disk_model <- fit_support_disk_model(kin, dim(cube_obj$cube)[1:2], slug, line$label)
  disk_residual_map <- if (!is.null(disk_model)) disk_model$residual else matrix(NA_real_, nrow(kin$valid), ncol(kin$valid))
  disk_residual_abs <- ifelse(kin$valid, abs(disk_residual_map), NA_real_)
  disk_residual_z <- positive_robust_z_map(disk_residual_abs, kin$valid)
  non_gaussian_wing_velocity <- 500 * (kin$red_wing - kin$blue_wing)
  non_gaussian_wing_velocity[!kin$valid] <- NA_real_

  flux_confidence <- rank_probability_map(log10(pmax(kin$flux, 0) + 1e-6), kin$valid)
  snr_confidence <- rank_probability_map(log1p(pmax(kin$snr, 0)), kin$valid)
  support_confidence <- sqrt(pmax(flux_confidence, 0) * pmax(snr_confidence, 0))
  support_confidence <- smooth_map_laplacian(support_confidence, kin$valid, steps = 2L, alpha = segmentation_laplacian_alpha)

  profile_outflow_score <-
    1.50 * positive_robust_z_map(kin$gaussian_wing_excess, kin$valid) +
    1.00 * positive_robust_z_map(kin$wing_fraction, kin$valid) +
    0.80 * positive_robust_z_map(abs(kin$wing_asymmetry), kin$valid) +
    0.55 * positive_robust_z_map(kin$sigma, kin$valid) +
    0.35 * positive_robust_z_map(abs(non_gaussian_wing_velocity), kin$valid)
  profile_outflow_score <- smooth_map_laplacian(profile_outflow_score, kin$valid, steps = 2L, alpha = segmentation_laplacian_alpha)

  disk_subtracted_outflow_score <- (profile_outflow_score + outflow_disk_residual_probability_weight * disk_residual_z) * support_confidence
  disk_subtracted_outflow_score[!kin$valid] <- NA_real_
  outflow_probability <- rank_probability_map(disk_subtracted_outflow_score, kin$valid)
  outflow_probability <- smooth_map_laplacian(
    outflow_probability,
    kin$valid,
    steps = outflow_probability_smoothing_steps,
    alpha = segmentation_laplacian_alpha
  )
  outflow_probability <- pmin(pmax(outflow_probability, 0), 1)
  p_disk_resid_abs <- save_map(disk_residual_abs, file.path(output_dir, paste0(slug, "_support_disk_residual_abs.png")), paste(line$label, "|velocity - disk|"), limits = robust_limits(disk_residual_abs, probs = c(0.02, 0.96)))
  p_disk_sub_score <- save_map(disk_subtracted_outflow_score, file.path(output_dir, paste0(slug, "_disk_subtracted_outflow_score.png")), paste(line$label, "disk-subtracted outflow score"), limits = robust_limits(disk_subtracted_outflow_score, probs = c(0.02, 0.96)))
  p_outflow_prob <- save_map(outflow_probability, file.path(output_dir, paste0(slug, "_outflow_probability.png")), paste(line$label, "relative outflow probability"), limits = c(0, 1))
  p_non_gauss_vel <- save_map(non_gaussian_wing_velocity, file.path(output_dir, paste0(slug, "_non_gaussian_wing_velocity.png")), paste(line$label, "non-Gaussian wing velocity"), diverging = TRUE, limits = robust_limits(non_gaussian_wing_velocity, symmetric = TRUE))

  interior <- erode_mask(kin$valid, iterations = 1L)
  flux_thr <- stats::quantile(kin$flux[kin$valid], outflow_candidate_flux_quantile, na.rm = TRUE)
  snr_thr <- stats::quantile(kin$snr[kin$valid], outflow_candidate_snr_quantile, na.rm = TRUE)
  wing_thr <- stats::quantile(kin$wing_fraction[kin$valid], outflow_candidate_wing_quantile, na.rm = TRUE)
  gaussian_thr <- stats::quantile(kin$gaussian_wing_excess[kin$valid], outflow_candidate_gaussian_excess_quantile, na.rm = TRUE)
  disk_resid_thr <- stats::quantile(disk_residual_abs[kin$valid], outflow_candidate_disk_residual_quantile, na.rm = TRUE)

  make_candidate_pool <- function(mask_gate) {
    profile_evidence <- kin$wing_fraction >= wing_thr | kin$gaussian_wing_excess >= gaussian_thr
    disk_evidence <- disk_residual_abs >= disk_resid_thr
    mask_gate &
      kin$flux >= flux_thr &
      kin$snr >= snr_thr &
      (profile_evidence | disk_evidence) &
      is.finite(disk_subtracted_outflow_score)
  }

  candidate_pool <- make_candidate_pool(interior)
  if (sum(candidate_pool, na.rm = TRUE) < outflow_min_candidate_pixels) {
    candidate_pool <- make_candidate_pool(kin$valid)
  }
  if (sum(candidate_pool, na.rm = TRUE) >= 5L) {
    candidate_thr <- max(outflow_probability_threshold, stats::quantile(outflow_probability[candidate_pool], outflow_candidate_quantile, na.rm = TRUE))
    candidate <- candidate_pool & outflow_probability >= candidate_thr
    if (sum(candidate, na.rm = TRUE) < outflow_min_candidate_pixels) {
      candidate_thr <- stats::quantile(outflow_probability[candidate_pool], 1 - outflow_min_candidate_pixels / sum(candidate_pool, na.rm = TRUE), na.rm = TRUE)
      candidate <- candidate_pool & outflow_probability >= candidate_thr
    }
  } else {
    candidate <- matrix(FALSE, nrow = nrow(kin$valid), ncol = ncol(kin$valid))
  }
  p_candidate <- save_map(ifelse(candidate, outflow_probability, NA_real_), file.path(output_dir, paste0(slug, "_path_outflow_candidate_mask.png")), paste(line$label, "likely outflow zone"), limits = c(0, 1))
  profile_atlas <- save_profile_atlas(kin, outflow_probability, candidate, slug, line$label)

  panel <- gridExtra::arrangeGrob(p_flux, p_snr, p_support, p_vel, p_sig, p_seg, p_outflow_prob, p_path_seg, p_candidate, ncol = 3)
  ggsave(file.path(output_dir, paste0(slug, "_standard_kinematic_panel.png")), panel, width = 11.5, height = 7.2, dpi = 320, bg = "white")
  signed_panel <- gridExtra::arrangeGrob(p_signed_blue, p_signed_red, p_signed_response, p_red_wing, p_blue_wing, p_wing_asym, ncol = 3)
  ggsave(file.path(output_dir, paste0(slug, "_signed_path_response_panel.png")), signed_panel, width = 11.5, height = 4.9, dpi = 320, bg = "white")
  outflow_panel <- gridExtra::arrangeGrob(p_outflow_prob, p_disk_sub_score, p_candidate, p_gaussian_excess, p_signed_response, p_non_gauss_vel, p_h3, p_h4, p_wing_asym, ncol = 3)
  ggsave(file.path(output_dir, paste0(slug, "_path_outflow_diagnostic_panel.png")), outflow_panel, width = 11.5, height = 7.2, dpi = 320, bg = "white")

  tab <- expand.grid(x = seq_len(dim(cube_obj$cube)[1]), y = seq_len(dim(cube_obj$cube)[2]))
  idx <- cbind(tab$x, tab$y)
  tab$support <- kin$support[idx]
  tab$valid <- kin$valid[idx]
  tab$flux <- kin$flux[idx]
  tab$snr <- kin$snr[idx]
  tab$velocity <- kin$velocity[idx]
  tab$sigma <- kin$sigma[idx]
  tab$asymmetry <- kin$asymmetry[idx]
  tab$wing_fraction <- kin$wing_fraction[idx]
  tab$blue_wing <- kin$blue_wing[idx]
  tab$red_wing <- kin$red_wing[idx]
  tab$signed_blue_excess <- kin$signed_blue_excess[idx]
  tab$signed_red_excess <- kin$signed_red_excess[idx]
  tab$signed_path_response <- kin$signed_path_response[idx]
  tab$wing_asymmetry <- kin$wing_asymmetry[idx]
  tab$path_skewness <- kin$path_skewness[idx]
  tab$path_kurtosis <- kin$path_kurtosis[idx]
  tab$h3 <- kin$h3[idx]
  tab$h4 <- kin$h4[idx]
  tab$peak_offset <- kin$peak_offset[idx]
  tab$gaussian_wing_excess <- kin$gaussian_wing_excess[idx]
  tab$path_score <- kin$path_score[idx]
  tab$path_outflow_score <- kin$outflow_score[idx]
  tab$path_outflow_score_smooth <- kin$outflow_score_smooth[idx]
  tab$disk_residual_abs <- disk_residual_abs[idx]
  tab$disk_residual_z <- disk_residual_z[idx]
  tab$disk_subtracted_outflow_score <- disk_subtracted_outflow_score[idx]
  tab$outflow_probability <- outflow_probability[idx]
  tab$non_gaussian_wing_velocity <- non_gaussian_wing_velocity[idx]
  tab$path_outflow_candidate_pool <- candidate_pool[idx]
  tab$path_outflow_candidate <- candidate[idx]
  tab$kinematic_segment <- seg$cluster_map[idx]
  tab$path_segment <- path_seg$cluster_map[idx]
  if (isTRUE(spectropath_features$available) && n_spectropath_features) {
    for (ff in spectropath_features$features) {
      tab[[paste0("spectropath_", ff)]] <- spectropath_features$feature_cube[, , ff][idx]
    }
  }
  if (!is.null(disk_model)) {
    tab$support_disk_model <- disk_model$model[idx]
    tab$support_disk_residual <- disk_model$residual[idx]
  }
  utils::write.csv(tab, file.path(output_dir, paste0(slug, "_standard_kinematic_spaxels.csv")), row.names = FALSE)

  list(line = line, kinematics = kin, segmentation = seg, path_segmentation = path_seg, spectropath_features = spectropath_features, disk_model = disk_model, profile_atlas = profile_atlas, panel = panel)
}

message("Reading cube...")
cube_obj <- read_cube(cube_path)
results <- lapply(lines, function(line) run_one_line(cube_obj, line))
names(results) <- vapply(lines, `[[`, character(1), "slug")
saveRDS(results, file.path(output_dir, paste0(tolower(object_id), "_standard_kinematics_results.rds")))
message("Done. Outputs saved in: ", output_dir)
