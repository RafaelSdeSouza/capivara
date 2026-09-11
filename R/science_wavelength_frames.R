# Coordinates stay on the native channel grid. No flux interpolation occurs here.
.valid_systemic_redshift <- function(z) {
  is.numeric(z) && length(z) == 1L && is.finite(z) && z > -1
}

.input_wavelength_frame <- function(cubedat, wavelength_frame = NULL, required = FALSE) {
  stored <- cubedat$wavelength_frame
  if (!is.null(wavelength_frame)) {
    wavelength_frame <- match.arg(wavelength_frame, c("observed", "rest"))
    if (!is.null(stored) && !identical(stored, wavelength_frame)) {
      stop("`wavelength_frame` conflicts with input$wavelength_frame; correct the metadata explicitly.", call. = FALSE)
    }
    stored <- wavelength_frame
  }
  if (!is.null(stored)) return(match.arg(stored, c("observed", "rest")))
  if (required) stop("Specify `wavelength_frame` ('observed' or 'rest') or input$wavelength_frame.", call. = FALSE)
  "unknown"
}

.subset_cubedat_wavelength_range <- function(cubedat, feature_wavelength_range = NULL,
                                             wavelength_frame = NULL,
                                             feature_wavelength_frame = NULL,
                                             redshift = NA_real_) {
  cubedat <- .as_cubedat(cubedat)
  cube <- cubedat$imDat
  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("`input$imDat` must be a 3D array with dimensions (n_row, n_col, n_wave).")
  }
  if (!is.null(feature_wavelength_range) && is.null(feature_wavelength_frame)) {
    stop("Ambiguous historical wavelength selection: explicitly set `feature_wavelength_frame` to 'observed' or 'rest'. Historical bounds selected native channels; redshift was unused.", call. = FALSE)
  }
  if (!is.null(feature_wavelength_frame)) {
    feature_wavelength_frame <- match.arg(feature_wavelength_frame, c("observed", "rest"))
  }
  input_frame <- .input_wavelength_frame(cubedat, wavelength_frame,
                                          required = !is.null(feature_wavelength_frame))
  valid_z <- .valid_systemic_redshift(redshift)
  if (identical(feature_wavelength_frame, "rest") && !valid_z) {
    stop("Rest-frame feature selection requires an explicit finite scalar `redshift` > -1 (including z = 0).", call. = FALSE)
  }
  if (!(is.null(redshift) || (is.numeric(redshift) && length(redshift) == 1L && is.na(redshift))) && !valid_z) {
    stop("`redshift` must be a finite scalar > -1, or NA when unused.", call. = FALSE)
  }
  n_wave <- dim(cube)[3]
  physical_axis <- !is.null(cubedat$wavelength) || !is.null(cubedat$axDat)
  if (!physical_axis && (!is.null(feature_wavelength_range) || !is.null(feature_wavelength_frame))) {
    stop("Physical wavelength selection requires input$wavelength or valid spectral axis metadata.", call. = FALSE)
  }
  wavelengths <- .wavelength_axis(cubedat, n_wave)
  native_range <- feature_wavelength_range
  if (!is.null(native_range)) {
    if (!is.numeric(native_range) || length(native_range) != 2L ||
        any(!is.finite(native_range)) || any(native_range <= 0) || native_range[1] > native_range[2]) {
      stop("`feature_wavelength_range` must be two increasing positive finite wavelengths.", call. = FALSE)
    }
    if (!identical(input_frame, feature_wavelength_frame)) {
      if (!valid_z) stop("Converting between wavelength frames requires a valid `redshift`.", call. = FALSE)
      native_range <- if (input_frame == "observed") native_range * (1 + redshift) else native_range / (1 + redshift)
    }
  }
  wave_idx <- .wavelength_range_index(cubedat, n_wave, native_range, "feature_wavelength_range")
  native_wave <- wavelengths[wave_idx]
  rest_wave <- if (!physical_axis) rep(NA_real_, length(wave_idx)) else if (input_frame == "rest") {
    native_wave
  } else if (input_frame == "observed" && valid_z) native_wave / (1 + redshift) else rep(NA_real_, length(wave_idx))
  bounds <- function(x) if (all(is.na(x))) c(NA_real_, NA_real_) else range(x)
  nb <- if (physical_axis) bounds(native_wave) else c(NA_real_, NA_real_)
  rb <- bounds(rest_wave)
  provenance <- list(
    input_wavelength_frame = input_frame,
    requested_feature_wavelength_range = feature_wavelength_range,
    requested_feature_wavelength_frame = if (is.null(feature_wavelength_frame)) input_frame else feature_wavelength_frame,
    systemic_redshift = if (valid_z) redshift else NA_real_,
    selected_native_wavelength_min = nb[1], selected_native_wavelength_max = nb[2],
    selected_rest_wavelength_min = rb[1], selected_rest_wavelength_max = rb[2],
    selected_channel_indices = wave_idx, number_of_selected_channels = length(wave_idx),
    native_selection_range = native_range,
    wavelength_medium = if (is.null(cubedat$wavelength_medium)) "unspecified" else cubedat$wavelength_medium,
    coordinate_kind = if (physical_axis) "wavelength" else "channel_or_feature_index",
    resampled = FALSE, channel_index_convention = "one-based input spectral axis"
  )
  out <- cubedat
  if (input_frame != "unknown") out$wavelength_frame <- input_frame
  if (!is.null(feature_wavelength_range)) {
    out$imDat <- cube[, , wave_idx, drop = FALSE]
    out$wavelength <- native_wave
    # The explicit vector is authoritative, including for logarithmic grids.
    out$axDat <- NULL
  }
  list(cubedat = out, wave_idx = wave_idx, selected_wavelengths = native_wave,
       provenance = provenance)
}

.subset_variance_channels <- function(var_cube, full_input, wave_idx) {
  out <- .as_cubedat(var_cube)
  if (!identical(dim(out$imDat), dim(full_input$imDat))) {
    stop("`var_cube` must have the same dimensions as the full input cube.")
  }
  out$imDat <- out$imDat[, , wave_idx, drop = FALSE]
  out$wavelength <- .wavelength_axis(full_input, dim(full_input$imDat)[3])[wave_idx]
  out$axDat <- NULL
  out
}
