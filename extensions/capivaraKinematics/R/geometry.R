#' Estimate basic disc geometry for Capivara kinematic modelling
#'
#' @param spaxels Spaxel-level data frame.
#' @param geometry List with optional x0, y0, pa_deg, inc_deg, and vsys values.
#' @param allow_placeholder_inclination If TRUE and inc_deg is missing, use
#'   `placeholder_inc_deg` and flag this in `geometry_status`.
#' @param placeholder_inc_deg Placeholder inclination in degrees.
#' @return A list with geometry values and status strings.
#' @export
estimate_disc_geometry <- function(spaxels,
                                   geometry = list(),
                                   allow_placeholder_inclination = FALSE,
                                   placeholder_inc_deg = 60) {
  get <- function(name) {
    value <- geometry[[name]]
    if (is.null(value) || length(value) == 0L || is.na(value)) NA_real_ else as.numeric(value)
  }
  get_character <- function(name, default) {
    value <- geometry[[name]]
    if (is.null(value) || length(value) == 0L || is.na(value)) default else as.character(value[[1]])
  }

  status <- character()
  x0 <- get("x0")
  y0 <- get("y0")
  if (!is.finite(x0) || !is.finite(y0)) {
    if ("flux" %in% names(spaxels) && any(is.finite(spaxels$flux))) {
      hit <- which.max(ifelse(is.finite(spaxels$flux), spaxels$flux, -Inf))
      x0 <- spaxels$x[hit]
      y0 <- spaxels$y[hit]
      status <- c(status, "centre_estimated_from_flux_peak")
    } else {
      x0 <- stats::median(spaxels$x, na.rm = TRUE)
      y0 <- stats::median(spaxels$y, na.rm = TRUE)
      status <- c(status, "centre_estimated_from_map_centre")
    }
  } else {
    status <- c(status, "centre_supplied")
  }

  vsys <- get("vsys")
  if (!is.finite(vsys)) {
    vsys <- stats::median(spaxels$velocity[spaxels$valid], na.rm = TRUE)
    status <- c(status, "vsys_estimated_from_valid_median")
  } else {
    status <- c(status, "vsys_supplied")
  }

  pa_deg <- get("pa_deg")
  if (!is.finite(pa_deg)) {
    valid <- spaxels$valid & is.finite(spaxels$velocity)
    if (sum(valid) >= 5L) {
      fit <- stats::lm(velocity ~ x + y, data = spaxels[valid, , drop = FALSE])
      b <- stats::coef(fit)
      pa_deg <- (atan2(b[["x"]], b[["y"]]) * 180 / pi) %% 180
      status <- c(status, "pa_estimated_from_velocity_gradient")
    } else {
      stop("PA is missing and there are too few valid spaxels to estimate it.", call. = FALSE)
    }
  } else {
    status <- c(status, "pa_supplied")
  }

  inc_deg <- get("inc_deg")
  if (!is.finite(inc_deg)) {
    if (!isTRUE(allow_placeholder_inclination)) {
      stop(
        "Inclination `inc_deg` is missing. Supply it in the config or set ",
        "`allow_placeholder_inclination: true` for exploratory runs.",
        call. = FALSE
      )
    }
    inc_deg <- placeholder_inc_deg
    status <- c(status, "inclination_placeholder")
  } else {
    status <- c(status, "inclination_supplied")
  }

  if (!is.finite(inc_deg) || inc_deg <= 1 || inc_deg >= 89.5) {
    stop("Inclination must be finite and safely between 1 and 89.5 degrees.", call. = FALSE)
  }

  list(
    x0 = x0,
    y0 = y0,
    vsys = vsys,
    pa_deg = pa_deg,
    inc_deg = inc_deg,
    pa_rad = pa_deg * pi / 180,
    inc_rad = inc_deg * pi / 180,
    coordinate_convention = get_character("coordinate_convention", "nirvana"),
    geometry_status = paste(status, collapse = ";")
  )
}

#' Deproject IFU coordinates using the Capivara kinematic convention
#'
#' The default convention follows NIRVANA's `projected_polar`: PA is measured
#' so that theta = 0 lies along the supplied PA, matching the barred-galaxy
#' paper. Set `geometry$coordinate_convention = "capivara_legacy"` to use the
#' older exploratory convention.
#'
#' @param x,y Pixel coordinates.
#' @param geometry Geometry list from `estimate_disc_geometry()`.
#' @return A data frame with X, Y, Yd, R, and theta.
#' @export
deproject_coordinates <- function(x, y, geometry) {
  dx <- x - geometry$x0
  dy <- y - geometry$y0
  pa <- geometry$pa_rad
  inc <- geometry$inc_rad
  convention <- geometry$coordinate_convention
  if (is.null(convention) || !length(convention)) {
    convention <- "nirvana"
  }
  convention <- tolower(as.character(convention[[1]]))
  if (identical(convention, "capivara_legacy")) {
    X <- -dx * sin(pa) + dy * cos(pa)
    Y <- -dx * cos(pa) - dy * sin(pa)
    Yd <- Y / cos(inc)
    theta <- atan2(Yd, X)
  } else if (identical(convention, "nirvana")) {
    X <- dx * sin(pa) + dy * cos(pa)
    Y <- dy * sin(pa) - dx * cos(pa)
    Yd <- Y / cos(inc)
    theta <- atan2(-Yd, X) %% (2 * pi)
  } else {
    stop("Unsupported coordinate_convention: ", convention, call. = FALSE)
  }
  data.frame(
    X = X,
    Y = Y,
    Yd = Yd,
    R = sqrt(X^2 + Yd^2),
    theta = theta
  )
}
