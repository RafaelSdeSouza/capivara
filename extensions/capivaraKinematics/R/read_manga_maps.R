#' Read MaNGA DAP MAPS velocity products
#'
#' @param maps_file MaNGA DAP MAPS FITS file. `.fits.gz` files are decompressed
#'   to a temporary file for FITSio.
#' @param velocity_component `stellar`, `halpha`, or `ha`.
#' @return A list with velocity, velocity_error, mask, flux, snr, header, and
#'   HDU metadata.
#' @export
read_manga_maps <- function(maps_file, velocity_component = "stellar") {
  .capivara_require("FITSio")

  if (!file.exists(maps_file)) {
    stop("MaNGA MAPS file does not exist: ", maps_file, call. = FALSE)
  }

  fits_path <- .fits_temp_unzip(maps_file)
  on.exit(if (!identical(fits_path, maps_file)) unlink(fits_path), add = TRUE)
  hdu_index <- .fits_hdu_index(fits_path)

  component <- tolower(velocity_component)
  line_channel <- NA_integer_
  if (component %in% c("stellar", "star", "stars")) {
    vel <- .fits_read_ext(fits_path, hdu_index, "STELLAR_VEL")
    ivar <- try(.fits_read_ext(fits_path, hdu_index, "STELLAR_VEL_IVAR"), silent = TRUE)
    mask <- try(.fits_read_ext(fits_path, hdu_index, "STELLAR_VEL_MASK"), silent = TRUE)
    velocity <- .matrix_from_fits_image(vel)
    velocity_ivar <- if (inherits(ivar, "try-error")) NULL else .matrix_from_fits_image(ivar)
    velocity_mask <- if (inherits(mask, "try-error")) NULL else .matrix_from_fits_image(mask)
    header <- vel$hdr
  } else if (component %in% c("halpha", "ha", "h-alpha", "gas_halpha")) {
    vel <- .fits_read_ext(fits_path, hdu_index, "EMLINE_GVEL")
    ch <- .fits_channel_index(vel$hdr, c("^Ha-?656", "Halpha", "H-alpha"))
    if (!is.finite(ch)) {
      stop("Could not identify the H-alpha channel in EMLINE_GVEL.", call. = FALSE)
    }
    line_channel <- ch
    velocity <- vel$imDat[, , ch]
    ivar <- try(.fits_read_ext(fits_path, hdu_index, "EMLINE_GVEL_IVAR"), silent = TRUE)
    mask <- try(.fits_read_ext(fits_path, hdu_index, "EMLINE_GVEL_MASK"), silent = TRUE)
    velocity_ivar <- if (inherits(ivar, "try-error")) NULL else ivar$imDat[, , ch]
    velocity_mask <- if (inherits(mask, "try-error")) NULL else mask$imDat[, , ch]
    header <- vel$hdr
  } else {
    stop("Unsupported velocity_component: ", velocity_component, call. = FALSE)
  }

  flux <- try(.fits_read_ext(fits_path, hdu_index, "SPX_MFLUX"), silent = TRUE)
  snr <- try(.fits_read_ext(fits_path, hdu_index, "SPX_SNR"), silent = TRUE)
  line_flux <- NULL
  line_sigma <- NULL
  if (is.finite(line_channel)) {
    eflux <- try(.fits_read_ext(fits_path, hdu_index, "EMLINE_GFLUX"), silent = TRUE)
    esigma <- try(.fits_read_ext(fits_path, hdu_index, "EMLINE_GSIGMA"), silent = TRUE)
    if (!inherits(eflux, "try-error")) {
      line_flux <- eflux$imDat[, , line_channel]
    }
    if (!inherits(esigma, "try-error")) {
      line_sigma <- esigma$imDat[, , line_channel]
    }
  }
  velocity_error <- NULL
  if (!is.null(velocity_ivar)) {
    velocity_error <- matrix(NA_real_, nrow(velocity_ivar), ncol(velocity_ivar))
    ok <- is.finite(velocity_ivar) & velocity_ivar > 0
    velocity_error[ok] <- 1 / sqrt(velocity_ivar[ok])
  }

  list(
    velocity = as.matrix(velocity),
    velocity_error = velocity_error,
    velocity_ivar = velocity_ivar,
    mask = velocity_mask,
    flux = if (inherits(flux, "try-error")) NULL else .matrix_from_fits_image(flux),
    snr = if (inherits(snr, "try-error")) NULL else .matrix_from_fits_image(snr),
    line_flux = line_flux,
    line_sigma = line_sigma,
    header = header,
    hdu_index = hdu_index,
    maps_file = maps_file,
    velocity_component = velocity_component
  )
}
