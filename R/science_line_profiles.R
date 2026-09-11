# Systemic optical Doppler coordinate, evaluated at native samples.
# The line and axis must use the same air/vacuum convention.
.systemic_line_coordinate <- function(wave, rest_wave, redshift, wavelength_frame,
                                        line_name, systemic_redshift_source,
                                        velocity_window_kms, profile_centering_mode,
                                        wavelength_medium, input_wavelength_medium = wavelength_medium) {
  wavelength_frame <- match.arg(wavelength_frame, c("observed", "rest"))
  profile_centering_mode <- match.arg(profile_centering_mode, c("systemic", "local_centroid"))
  wavelength_medium <- match.arg(wavelength_medium, c("air", "vacuum"))
  input_wavelength_medium <- match.arg(input_wavelength_medium, c("air", "vacuum"))
  if (input_wavelength_medium != wavelength_medium) {
    stop("Laboratory and input wavelength medium differ; convert explicitly before redshift/velocity calculation.")
  }
  registry <- .capivara_match_emission_lines(line_name, wavelength_medium)
  if (nrow(registry) != 1L || abs(rest_wave - registry$rest_wavelength) > 1e-6) {
    stop("Laboratory wavelength does not match the line registry in the declared medium (possible air/vacuum mismatch).")
  }
  if (!.valid_systemic_redshift(redshift)) stop("A finite systemic `redshift` > -1 is required.")
  if (!is.numeric(rest_wave) || length(rest_wave) != 1L || !is.finite(rest_wave) || rest_wave <= 0) {
    stop("A positive finite laboratory `rest_wave` is required.")
  }
  if (!is.numeric(wave) || length(wave) < 2L || any(!is.finite(wave)) || any(diff(wave) <= 0)) {
    stop("`wave` must contain increasing finite native wavelengths.")
  }
  if (length(velocity_window_kms) != 1L || !is.finite(velocity_window_kms) || velocity_window_kms <= 0) {
    stop("`velocity_window_kms` must be a positive finite half-width.")
  }
  if (length(systemic_redshift_source) != 1L || is.na(systemic_redshift_source) || !nzchar(systemic_redshift_source)) {
    stop("A nonempty `systemic_redshift_source` is required.")
  }
  lab_vacuum <- if (wavelength_medium == "air") convert_wavelength_medium(rest_wave, "air", "vacuum") else rest_wave
  wave_vacuum <- if (wavelength_medium == "air") convert_wavelength_medium(wave, "air", "vacuum") else wave
  observed_centre <- lab_vacuum * (1 + redshift)
  if (wavelength_medium == "air") observed_centre <- convert_wavelength_medium(observed_centre, "vacuum", "air")
  native_centre <- if (wavelength_frame == "observed") observed_centre else rest_wave
  native_vacuum_centre <- if (wavelength_frame == "observed") lab_vacuum*(1+redshift) else lab_vacuum
  velocity <- 299792.458 * (wave_vacuum / native_vacuum_centre - 1)
  ix <- which(abs(velocity) <= velocity_window_kms)
  rest <- if (wavelength_frame == "observed") wave_vacuum / (1 + redshift) else wave_vacuum
  if (wavelength_medium == "air") rest <- convert_wavelength_medium(rest,"vacuum","air")
  native_bounds <- if (length(ix)) range(wave[ix]) else c(NA_real_, NA_real_)
  rest_bounds <- if (length(ix)) range(rest[ix]) else c(NA_real_, NA_real_)
  list(velocity = velocity,
    wavelength_provenance = list(input_wavelength_frame = wavelength_frame,
      requested_feature_wavelength_range = rest_wave * (1 + c(-1, 1) * velocity_window_kms / 299792.458),
      requested_feature_wavelength_frame = "rest", systemic_redshift = redshift,
      selected_native_wavelength_min = native_bounds[1], selected_native_wavelength_max = native_bounds[2],
      selected_rest_wavelength_min = rest_bounds[1], selected_rest_wavelength_max = rest_bounds[2],
      selected_channel_indices = ix, number_of_selected_channels = length(ix),
      wavelength_medium = wavelength_medium, resampled = FALSE),
    provenance = list(line_name = line_name, line_rest_wavelength = rest_wave,
      systemic_redshift = redshift, systemic_redshift_source = systemic_redshift_source,
      input_wavelength_frame = wavelength_frame, wavelength_medium = wavelength_medium,
      input_wavelength_medium = input_wavelength_medium, velocity_wavelength_medium = "vacuum",
      medium_conversion = if (wavelength_medium == "air") "explicit air-to-vacuum before redshift and optical velocity" else "none",
      line_identifier = registry$name,
      line_reference = registry$reference, line_registry_version = registry$registry_version,
      observed_line_centre = observed_centre, native_line_centre = native_centre,
      velocity_definition = "optical: c * (lambda_obs / (lambda0 * (1 + z_sys)) - 1); c = 299792.458 km/s",
      velocity_window_kms = c(-velocity_window_kms, velocity_window_kms),
      profile_centering_mode = profile_centering_mode, resampled = FALSE))
}

.centre_velocity_profile <- function(velocity, profile, mode = c("systemic", "local_centroid")) {
  mode <- match.arg(mode)
  path <- spectropath::as_spectral_path(cbind(velocity, profile), remove_nonfinite = FALSE)
  offset <- if (mode == "local_centroid") spectropath::line_moments(path)$centroid else 0
  if (!is.finite(offset)) stop("Local centring requires a finite positive-profile centroid.")
  path[, 1] <- path[, 1] - offset
  list(path = path, centroid_offset_kms = offset)
}

compute_line_maps <- function(cube,
                              wave,
                              mask,
                              z,
                              rest_wave,
                              line_name,
                              line_window_kms = 600,
                              peak_search_kms = 350,
                              centroid_window_kms = 260,
                              cont_inner_kms = 800,
                              cont_outer_kms = 1400,
                              wavelength_frame = "observed",
                              systemic_redshift_source = "explicit argument",
                              wavelength_medium = "vacuum") {
  coordinate <- .systemic_line_coordinate(wave, rest_wave, z, wavelength_frame,
    line_name, systemic_redshift_source, line_window_kms, "systemic", wavelength_medium)
  lambda0 <- coordinate$provenance$observed_line_centre
  vel <- coordinate$velocity
  line_idx <- which(abs(vel) <= line_window_kms)
  cont_idx <- which(abs(vel) > cont_inner_kms & abs(vel) <= cont_outer_kms)
  if (length(line_idx) < 5L || length(cont_idx) < 5L) {
    stop("Insufficient wavelength channels for ", line_name, " line/continuum windows.")
  }

  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  flux <- matrix(NA_real_, nx, ny)
  velocity <- matrix(NA_real_, nx, ny)
  sigma <- matrix(NA_real_, nx, ny)
  asymmetry <- matrix(NA_real_, nx, ny)
  h3_proxy <- matrix(NA_real_, nx, ny)
  h4_proxy <- matrix(NA_real_, nx, ny)

  for (idx in which(mask)) {
    ij <- arrayInd(idx, .dim = c(nx, ny))
    i <- ij[1]
    j <- ij[2]
    spec <- cube[i, j, ]
    if (!all(is.finite(spec[line_idx]))) next

    cont <- stats::median(spec[cont_idx], na.rm = TRUE)
    if (!is.finite(cont)) next
    line <- spec[line_idx] - cont
    v <- vel[line_idx]
    search <- abs(v) <= peak_search_kms
    if (sum(search) < 3L) {
      search <- rep(TRUE, length(v))
    }
    search_idx <- which(search)
    peak <- search_idx[which.max(line[search])]
    if (!length(peak) || !is.finite(line[peak]) || line[peak] <= 0) next

    local <- abs(v - v[peak]) <= centroid_window_kms
    if (sum(local) < 3L) {
      local <- rep(TRUE, length(v))
    }
    pos <- pmax(line, 0)
    pos[!local] <- 0
    sum_pos <- sum(pos)
    if (!is.finite(sum_pos) || sum_pos <= 0) next

    mu <- sum(v * pos) / sum_pos
    sig <- sqrt(sum(pos * (v - mu)^2) / sum_pos)
    blue <- sum(pos[v < 0])
    red <- sum(pos[v > 0])

    flux[i, j] <- sum_pos
    velocity[i, j] <- mu
    sigma[i, j] <- sig
    asymmetry[i, j] <- (red - blue) / (red + blue + .Machine$double.eps)
    if (is.finite(sig) && sig > 0) {
      zvel <- (v - mu) / sig
      h3_proxy[i, j] <- sum(pos * zvel^3) / sum_pos
      h4_proxy[i, j] <- sum(pos * zvel^4) / sum_pos - 3
    }
  }

  ok <- mask & is.finite(flux) & flux > 0 & is.finite(velocity) & is.finite(sigma)
  median_offset <- if (any(ok)) stats::median(velocity[ok], na.rm = TRUE) else NA_real_
  velocity_median_centered <- velocity - median_offset

  list(
    frame_provenance = coordinate$provenance,
    wavelength_provenance = coordinate$wavelength_provenance,
    lambda0 = lambda0,
    line_name = line_name,
    rest_wave = rest_wave,
    line_window_kms = line_window_kms,
    cont_inner_kms = cont_inner_kms,
    cont_outer_kms = cont_outer_kms,
    velocity_grid = vel,
    line_idx = line_idx,
    cont_idx = cont_idx,
    flux = flux,
    velocity = velocity,
    velocity_systemic = velocity,
    velocity_median_centered = velocity_median_centered,
    median_velocity_offset_kms = median_offset,
    sigma = sigma,
    asymmetry = asymmetry,
    h3_proxy = h3_proxy,
    h4_proxy = h4_proxy,
    valid = ok
  )
}

baseline_subtract <- function(v, y, n_edge = 2L) {
  edge <- c(seq_len(min(n_edge, length(y))), seq.int(max(1L, length(y) - n_edge + 1L), length(y)))
  base <- stats::median(y[edge], na.rm = TRUE)
  y - base
}

nonparam_velocity <- function(v, y) {
  line <- baseline_subtract(v, y)
  pos <- pmax(line, 0)
  flux <- sum(pos, na.rm = TRUE)
  if (!is.finite(flux) || flux <= 0) {
    return(data.frame(
      np_flux = NA_real_,
      np_peak_v = NA_real_,
      np_centroid = NA_real_,
      np_sigma = NA_real_,
      np_w80 = NA_real_
    ))
  }

  centroid <- sum(v * pos, na.rm = TRUE) / flux
  sigma <- sqrt(sum((v - centroid)^2 * pos, na.rm = TRUE) / flux)
  ord <- order(v)
  cum <- cumsum(pos[ord]) / flux
  qv <- stats::approx(cum, v[ord], xout = c(0.1, 0.5, 0.9), ties = "ordered", rule = 2)$y

  data.frame(
    np_flux = flux,
    np_peak_v = v[which.max(line)],
    np_centroid = centroid,
    np_sigma = sigma,
    np_w80 = qv[3] - qv[1]
  )
}

build_path_feature_cube <- function(cube,
                                    wave,
                                    observed_wave,
                                    mask,
                                    max_abs_velocity = 600,
                                    feature_names = c("p2", "p3u", "p3F", "p4F", "p4T", "p_pm"),
                                    rest_wave, redshift, line_name,
                                    wavelength_frame = "observed",
                                    systemic_redshift_source = "explicit argument",
                                    profile_centering_mode = c("systemic", "local_centroid"),
                                    wavelength_medium = "vacuum") {
  profile_centering_mode <- match.arg(profile_centering_mode)
  coordinate <- .systemic_line_coordinate(wave, rest_wave, redshift, wavelength_frame,
    line_name, systemic_redshift_source, max_abs_velocity, profile_centering_mode, wavelength_medium)
  if (!isTRUE(all.equal(observed_wave, coordinate$provenance$observed_line_centre, tolerance = 1e-12))) {
    stop("`observed_wave` must equal rest_wave * (1 + redshift); possible double application of redshift.")
  }
  vel <- coordinate$velocity
  keep <- abs(vel) <= max_abs_velocity
  selected_channel_indices <- which(keep)
  vel <- vel[keep]
  line_cube <- cube[, , keep, drop = FALSE]
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  feature_cube <- array(NA_real_, dim = c(nx, ny, length(feature_names)), dimnames = list(NULL, NULL, feature_names))
  rows <- vector("list", sum(mask, na.rm = TRUE))
  n <- 0L

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      if (!isTRUE(mask[i, j])) next
      flux <- as.numeric(line_cube[i, j, ])
      if (length(flux) < 8L || !all(is.finite(flux))) next
      line <- baseline_subtract(vel, flux)
      amp <- max(abs(line), na.rm = TRUE)
      if (!is.finite(amp) || amp <= 0) next
      profile <- line / amp
      represented <- .centre_velocity_profile(vel, profile, profile_centering_mode)
      pf <- tryCatch(
        spectropath::path_features(represented$path, depth = 4, normalize = TRUE, notation = "paper"),
        error = function(e) NULL
      )
      if (is.null(pf)) next
      vals <- as.numeric(pf[1, feature_names, drop = TRUE])
      feature_cube[i, j, ] <- vals

      np <- nonparam_velocity(vel, line)
      n <- n + 1L
      rows[[n]] <- data.frame(
        x = i,
        y = j,
        line_flux = np$np_flux,
        centroid_kms = np$np_centroid,
        profile_coordinate_offset_kms = represented$centroid_offset_kms,
        represented_centroid_kms = spectropath::line_moments(represented$path)$centroid,
        sigma_kms = np$np_sigma,
        w80_kms = np$np_w80,
        pf[, feature_names, drop = FALSE],
        check.names = FALSE
      )
    }
  }

  table <- if (n) do.call(rbind, rows[seq_len(n)]) else data.frame()

  list(
    frame_provenance = coordinate$provenance,
    wavelength_provenance = coordinate$wavelength_provenance,
    selected_channel_indices = selected_channel_indices,
    profile_centering_mode = profile_centering_mode,
    feature_cube = feature_cube,
    table = table,
    velocity = vel,
    line_cube = line_cube,
    features = feature_names
  )
}

