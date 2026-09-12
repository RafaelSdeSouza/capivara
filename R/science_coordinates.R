# Native image axes: x = column, y = row; axial PA from +y towards +x.
.capivara_axial_angle <- function(x) ((x + 90) %% 180) - 90

.capivara_image_to_disc_angle <- function(pa_image_deg, geometry) {
  if (length(pa_image_deg) != 1L || !is.finite(pa_image_deg)) {
    stop("Supply one finite pa_image_deg.", call. = FALSE)
  }
  a <- pa_image_deg * pi / 180
  p <- deproject_coordinates(geometry$x0 + sin(a), geometry$y0 + cos(a), geometry)
  .capivara_axial_angle(p$theta * 180 / pi)
}

# The callback must use the actual WCS and return RA/Dec in degrees.
.capivara_image_to_sky_angle <- function(pa_image_deg, x, y, pixel_to_sky) {
  if (!is.function(pixel_to_sky)) stop("A pixel-to-sky WCS function is required.")
  a <- pa_image_deg * pi / 180
  sky <- as.matrix(pixel_to_sky(c(x, x + sin(a)), c(y, y + cos(a))))
  if (!identical(dim(sky), c(2L, 2L)) || any(!is.finite(sky))) {
    stop("pixel_to_sky must return a finite 2x2 RA/Dec matrix in degrees.")
  }
  # Spherical initial bearing; works across RA=0 and with either WCS parity.
  dra <- (sky[2, 1] - sky[1, 1]) * pi / 180
  d0 <- sky[1, 2] * pi / 180
  d1 <- sky[2, 2] * pi / 180
  east <- sin(dra) * cos(d1)
  north <- cos(d0) * sin(d1) - sin(d0) * cos(d1) * cos(dra)
  if (abs(east) + abs(north) < 1e-15) stop("WCS maps the axis to zero length.")
  unname((atan2(east, north) * 180 / pi) %% 180)
}

.capivara_vetted_bar_geometry <- function(bar, geometry, dims, science = TRUE) {
  if (!is.list(bar)) stop("bar_geometry must be a named list.", call. = FALSE)
  if (science && !isTRUE(bar$vetted)) {
    stop("Scientific mode requires bar_geometry$vetted = TRUE.", call. = FALSE)
  }
  if (is.null(bar$source) || length(bar$source) != 1L ||
      is.na(bar$source) || !nzchar(bar$source)) {
    stop("bar_geometry requires a non-empty source.", call. = FALSE)
  }
  phi <- bar$phi_bar_disc_deg
  image_pa <- bar$pa_image_deg
  if (!is.null(image_pa)) {
    derived <- .capivara_image_to_disc_angle(image_pa, geometry)
    if (!is.null(phi) && (length(phi) != 1L || !is.finite(phi) ||
        abs(.capivara_axial_angle(phi - derived)) > 1e-6)) {
      stop("Supplied image and disc-plane bar angles disagree.", call. = FALSE)
    }
    phi <- derived
  }
  if (is.null(phi) || length(phi) != 1L || !is.finite(phi)) {
    stop("Supply pa_image_deg or phi_bar_disc_deg; sky PA alone needs a WCS conversion.")
  }
  if (!is.null(bar$axis_ratio) && (length(bar$axis_ratio) != 1L ||
      !is.finite(bar$axis_ratio) || bar$axis_ratio <= 0 || bar$axis_ratio > 1)) {
    stop("axis_ratio must be minor/major in (0, 1].")
  }
  if (!is.null(bar$bar_mask) && (!is.logical(bar$bar_mask) ||
      !identical(dim(bar$bar_mask), dims) || anyNA(bar$bar_mask))) {
    stop("bar_mask must be a non-missing logical matrix on the native image grid.")
  }
  bar$phi_bar_disc_deg <- .capivara_axial_angle(phi)
  bar$phi_bar_disc_rad <- bar$phi_bar_disc_deg * pi / 180
  if (is.null(bar$pa_sky_deg)) bar$pa_sky_deg <- NA_real_
  bar$bar_status <- if (isTRUE(bar$vetted)) "vetted_geometry_supplied" else "unvetted_preview_geometry"
  bar
}
