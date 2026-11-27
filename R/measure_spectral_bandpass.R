#' Measure a native-grid spectral bandpass
#'
#' Fit a linear continuum to two sidebands and integrate the excess flux in a
#' declared line window. The result describes the observed profile within that
#' window; blends, stellar absorption and flux outside the window remain in the
#' measurement. No kinematic or intrinsic-line inference is performed.
#'
#' @param wavelength Increasing native wavelengths in Angstrom.
#' @param flux Flux density at each native sample.
#' @param variance Per-sample variance in flux-density units squared. This must
#'   include the spatial covariance contribution if the spectrum was summed
#'   over spatial pixels and that contribution is known.
#' @param line_window Two wavelengths bounding the integration window.
#' @param continuum_windows A two-row matrix with lower and upper wavelength
#'   bounds for the blue and red continuum sidebands.
#' @param flux_unit Declared input flux-density unit.
#' @param spectral_covariance Optional full channel-by-channel covariance
#'   matrix, with diagonal equal to `variance`. If omitted, distinct spectral
#'   channels are treated as uncorrelated for uncertainty propagation.
#' @param lsf_sigma_angstrom Optional scalar or native-channel vector of
#'   instrumental Gaussian sigma in Angstrom. It is reported, not deconvolved.
#' @param wavelength_frame Either `"observed"` or `"rest"`; the window bounds
#'   must use this same frame. The equivalent width follows this frame.
#' @param wavelength_medium Either `"vacuum"` or `"air"`; the windows must use
#'   the same convention. No wavelength conversion is performed.
#' @return A `capivara_bandpass_measurement` list. A positive integrated excess
#'   denotes emission above the fitted continuum; a negative value denotes
#'   absorption or a deficit. The equivalent width has the same sign and is
#'   unavailable when the fitted continuum is nonpositive. Missing samples in
#'   either sideband or the line window return `incomplete_support` and missing
#'   measurements. The reported 2-by-2 covariance covers integrated excess and
#'   continuum at the window centre. Equivalent-width variance uses a local
#'   linear propagation. The result is descriptive, not a physical line fit.
#' @export
measure_spectral_bandpass <- function(wavelength, flux, variance,
                                      line_window, continuum_windows,
                                      flux_unit, spectral_covariance = NULL,
                                      lsf_sigma_angstrom = NULL,
                                      wavelength_frame = "observed",
                                      wavelength_medium = "vacuum") {
  n <- length(wavelength)
  if (!is.numeric(wavelength) || n < 6L || any(!is.finite(wavelength)) ||
      any(diff(wavelength) <= 0))
    stop("`wavelength` must be a finite increasing native grid.", call. = FALSE)
  if (!is.numeric(flux) || length(flux) != n ||
      !is.numeric(variance) || length(variance) != n)
    stop("`flux` and `variance` must match the wavelength grid.", call. = FALSE)
  if (!is.numeric(line_window) || length(line_window) != 2L ||
      any(!is.finite(line_window)) || line_window[1L] >= line_window[2L])
    stop("`line_window` must have increasing finite bounds.", call. = FALSE)
  if (!is.matrix(continuum_windows) || !is.numeric(continuum_windows) ||
      !identical(dim(continuum_windows), c(2L, 2L)) ||
      any(!is.finite(continuum_windows)) ||
      any(continuum_windows[, 1L] >= continuum_windows[, 2L]) ||
      continuum_windows[1L, 2L] >= line_window[1L] ||
      line_window[2L] >= continuum_windows[2L, 1L])
    stop("Supply two ordered, disjoint continuum sidebands around the line window.", call. = FALSE)
  if (line_window[1L] <= wavelength[1L] ||
      line_window[2L] >= wavelength[n])
    stop("The line window must lie inside the native wavelength grid.", call. = FALSE)
  if (!is.character(flux_unit) || length(flux_unit) != 1L ||
      is.na(flux_unit) || !nzchar(flux_unit))
    stop("Declare a nonempty `flux_unit`.", call. = FALSE)
  wavelength_frame <- match.arg(wavelength_frame, c("observed", "rest"))
  wavelength_medium <- match.arg(wavelength_medium, c("vacuum", "air"))
  edges <- c(wavelength[1L] - diff(wavelength)[1L] / 2,
             (wavelength[-1L] + wavelength[-n]) / 2,
             wavelength[n] + diff(wavelength)[n - 1L] / 2)
  line <- which(edges[-1L] > line_window[1L] &
                  edges[-(n + 1L)] < line_window[2L])
  blue <- which(wavelength >= continuum_windows[1L, 1L] &
                  wavelength <= continuum_windows[1L, 2L])
  red <- which(wavelength >= continuum_windows[2L, 1L] &
                 wavelength <= continuum_windows[2L, 2L])
  if (length(line) < 2L || length(blue) < 2L || length(red) < 2L ||
      length(intersect(line, c(blue, red))))
    stop("Windows must contain at least two distinct native samples each.", call. = FALSE)
  selected <- sort(unique(c(line, blue, red)))
  lsf_line <- if (is.null(lsf_sigma_angstrom)) NULL else
    if (length(lsf_sigma_angstrom) == 1L) lsf_sigma_angstrom else
      lsf_sigma_angstrom[line]
  if (!is.null(lsf_sigma_angstrom)) {
    if (!is.numeric(lsf_sigma_angstrom) ||
        !length(lsf_sigma_angstrom) %in% c(1L, n) ||
        any(!is.finite(lsf_line)) || any(lsf_line <= 0))
      stop("`lsf_sigma_angstrom` must be positive at every selected line sample.", call. = FALSE)
  }
  lsf <- if (is.null(lsf_sigma_angstrom)) NA_real_ else
    stats::median(lsf_line)
  covariance_model <- if (is.null(spectral_covariance))
    "diagonal spectral covariance assumed" else "supplied spectral covariance"
  result <- list(status = "ok", wavelength_frame = wavelength_frame,
                 wavelength_medium = wavelength_medium,
                 flux_density_unit = flux_unit,
                 integrated_unit = paste0("(", flux_unit, ") * Angstrom"),
                 equivalent_width_unit = "Angstrom",
                 line_window = line_window,
                 continuum_windows = continuum_windows,
                 line_channel_indices = line,
                 continuum_channel_indices = c(blue, red),
                 covariance_model = covariance_model,
                 instrumental_sigma_angstrom = lsf,
                 window_width_over_lsf_fwhm = if (is.na(lsf)) NA_real_ else
                   diff(line_window) / (2 * sqrt(2 * log(2)) * lsf),
                 resolution_policy = "Observed native profile; LSF not deconvolved; window may omit line wings",
                 interpretation = "Window-integrated excess above a linear sideband continuum; blends and absorption retained",
                 line_flux = NA_real_, line_flux_se = NA_real_,
                 continuum_density = NA_real_, continuum_density_se = NA_real_,
                 equivalent_width = NA_real_, equivalent_width_se = NA_real_,
                 measurement_covariance = matrix(NA_real_, 2L, 2L,
                   dimnames = list(c("line_flux", "continuum_density"),
                                   c("line_flux", "continuum_density"))))
  if (any(!is.finite(flux[selected])) ||
      any(!is.finite(variance[selected]) | variance[selected] <= 0)) {
    result$status <- "incomplete_support"
    return(structure(result, class = "capivara_bandpass_measurement"))
  }
  if (is.null(spectral_covariance)) {
    C <- diag(variance[selected], nrow = length(selected))
  } else {
    if (!is.matrix(spectral_covariance) || !is.numeric(spectral_covariance) ||
        !identical(dim(spectral_covariance), c(n, n)))
      stop("`spectral_covariance` must be a native-channel square matrix.", call. = FALSE)
    C <- spectral_covariance[selected, selected, drop = FALSE]
    if (any(!is.finite(C)) ||
        !isTRUE(all.equal(C, t(C), tolerance = 1e-9)) ||
        !isTRUE(all.equal(diag(C), variance[selected], tolerance = 1e-7)))
      stop("Spectral covariance must be finite, symmetric and match `variance`.", call. = FALSE)
    if (min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) <
        -1e-8 * max(diag(C)))
      stop("Spectral covariance must be positive semidefinite.", call. = FALSE)
  }
  centre <- mean(line_window)
  continuum <- c(blue, red)
  ci <- match(continuum, selected)
  li <- match(line, selected)
  Xc <- cbind(1, wavelength[continuum] - centre)
  Xl <- cbind(1, wavelength[line] - centre)
  CiX <- tryCatch(solve(C[ci, ci, drop = FALSE], Xc),
                  error = function(e) NULL)
  if (is.null(CiX))
    stop("Continuum covariance is singular.", call. = FALSE)
  H <- tryCatch(solve(crossprod(Xc, CiX), t(CiX)),
                error = function(e) NULL)
  if (is.null(H))
    stop("Continuum sidebands cannot identify a linear continuum.", call. = FALSE)
  beta <- as.vector(H %*% flux[continuum])
  width <- pmax(0, pmin(edges[line + 1L], line_window[2L]) -
                   pmax(edges[line], line_window[1L]))
  excess <- sum(width * (flux[line] - as.vector(Xl %*% beta)))
  g <- numeric(length(selected))
  g[li] <- width
  g[ci] <- -as.vector(t(width) %*% Xl %*% H)
  h <- numeric(length(selected))
  h[ci] <- H[1L, ]
  M <- cbind(g, h)
  V <- crossprod(M, C %*% M)
  dimnames(V) <- dimnames(result$measurement_covariance)
  continuum_at_centre <- beta[1L]
  result$line_flux <- excess
  result$line_flux_se <- sqrt(max(V[1L, 1L], 0))
  result$continuum_density <- continuum_at_centre
  result$continuum_density_se <- sqrt(max(V[2L, 2L], 0))
  result$measurement_covariance <- V
  if (continuum_at_centre > 0) {
    result$equivalent_width <- excess / continuum_at_centre
    gradient <- g / continuum_at_centre -
      excess * h / continuum_at_centre^2
    result$equivalent_width_se <- sqrt(max(
      as.numeric(crossprod(gradient, C %*% gradient)), 0))
  }
  structure(result, class = "capivara_bandpass_measurement")
}
