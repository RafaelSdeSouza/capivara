#' Read native MaNGA pre- and post-pixelization LSF cubes
#'
#' Reads the native DRP image sequence (FLUX, IVAR, MASK, post LSF, pre LSF,
#' WAVE), including that sequence in a MEGACUBE container. A later derived
#' extension named DISP is never used. Missing or conflicting native products
#' cause an error; no scalar resolution is substituted. Zero/nonfinite LSF
#' samples are marked missing, without interpolation.
#' @param path Native MaNGA cube or audited MEGACUBE FITS path.
#' @param flux Optional already-read FLUX object, in FITSio order (x,y,wavelength).
#' @return A list with `lsf_sigma_angstrom_pre`, `lsf_sigma_angstrom_post`,
#'   native `wavelength`, and provenance. Sigma arrays have exactly FLUX's order.
#' @export
read_manga_lsf <- function(path, flux = NULL) {
  if (is.null(flux)) flux <- .capivara_read_fits(path, hdu = 1)
  value <- function(h, k) unname(.fits_header_value(h$hdr, k, ""))
  if (value(flux, "INSTRUME") != "MaNGA" || value(flux, "EXTNAME") != "FLUX" ||
      !value(flux, "CTYPE3") %in% c("WAVE", "WAVE-LOG")) stop("Expected native MaNGA spectral FLUX.")
  drp <- value(flux, "VERSDRP3")
  if (!grepl("^v[23]_", drp)) stop("Unverified MaNGA DRP version: ", drp)
  post <- .capivara_read_fits(path, hdu = 4)
  pre <- .capivara_read_fits(path, hdu = 5)
  wave <- .capivara_read_fits(path, hdu = 6)
  names_native <- c(pre = value(pre, "EXTNAME"), post = value(post, "EXTNAME"))
  expected <- if (startsWith(drp, "v3_")) c(pre = "LSFPRE", post = "LSFPOST") else c(pre = "PREDISP", post = "DISP")
  if (!identical(names_native, expected)) stop("Native LSF extension sequence disagrees with DRP version.")
  if (value(wave, "EXTNAME") != "WAVE" || length(wave$imDat) != dim(flux$imDat)[3]) stop("LSF/FLUX WAVE length mismatch.")
  w <- as.numeric(wave$imDat)
  if (!identical(w, .wavelength_axis(flux, length(w)))) stop("LSF WAVE differs from FLUX wavelength axis.")
  .canonical_manga_lsf(flux, pre, post, w, drp, normalizePath(path), names_native)
}

.canonical_manga_lsf <- function(flux, pre, post, wavelength, drp_version, path, native_names) {
  value <- function(h, k) unname(.fits_header_value(h$hdr, k, ""))
  dims <- dim(flux$imDat)
  if (length(dims) != 3L || length(wavelength) != dims[3] ||
      any(!is.finite(wavelength)) || any(diff(wavelength) <= 0)) stop("Invalid LSF wavelength sampling.")
  arrays <- list(pre = pre, post = post)
  missing <- integer(2); names(missing) <- c("pre", "post")
  units <- character(2); names(units) <- names(missing)
  for (kind in names(arrays)) {
    h <- arrays[[kind]]
    if (!identical(dim(h$imDat), dims)) stop("Native LSF dimensions/orientation differ from FLUX.")
    units[kind] <- value(h, "BUNIT")
    if (nzchar(units[kind]) && !tolower(trimws(units[kind])) %in% c("angstrom", "angstroms", "ang", "a")) {
      stop("LSF units conflict with the native sigma_lambda Angstrom contract.")
    }
    # Some containers strip all LSF WCS cards. Any surviving cards must agree.
    for (key in c("CTYPE1", "CTYPE2", "CTYPE3", "CRVAL1", "CRVAL2", "CRVAL3",
                  "CRPIX1", "CRPIX2", "CRPIX3", "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CD3_3")) {
      a <- value(h, key); b <- value(flux, key)
      if (nzchar(a) && !identical(a, b)) stop("LSF WCS differs from FLUX: ", key)
    }
    bad <- !is.finite(h$imDat) | h$imDat <= 0
    missing[kind] <- sum(bad)
    if (all(bad)) stop("Native LSF contains no positive finite sigma_lambda.")
    h$imDat[bad] <- NA_real_
    arrays[[kind]] <- h$imDat
  }
  list(lsf_sigma_angstrom_pre = arrays$pre, lsf_sigma_angstrom_post = arrays$post,
    wavelength = wavelength,
    provenance = list(native_extension = native_names, native_hdu = c(pre = 5L, post = 4L),
      drp_version = drp_version, source_path = path, native_bunit = units,
      units = "Angstrom", quantity = "Gaussian sigma_lambda", wavelength_medium = "vacuum",
      wavelength_frame = "observed", orientation = "FITSio (x, y, wavelength); identical to FLUX",
      dimensions = dims, wavelength_hdu = "WAVE", sampling = "native WAVE, no resampling",
      pixelization = c(pre = "before detector pixel integration", post = "includes detector pixel integration"),
      invalid_samples = missing, scalar_fallback = FALSE,
      units_source = "SDSS native DRP data model; checked against BUNIT when present",
      reference = "https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/"))
}

#' Select a MaNGA LSF for an explicitly specified fitting method
#'
#' A Gaussian evaluated at pixel centres uses POST. A template model that
#' includes pixel integration uses PRE. Template convolution alone does not
#' establish that pixel integration has been included. This function never
#' chooses a default fitting method or a scalar FWHM.
#' @param lsf Result from `read_manga_lsf()`.
#' @param fitting_method `direct_profile` or `template_convolution`.
#' @param model_pixel_sampling `point_sampled` or `pixel_integrated`, required.
#' @return Selected sigma array, wavelength, and the explicit fitting contract.
#' @export
select_manga_lsf <- function(lsf, fitting_method, model_pixel_sampling) {
  fitting_method <- match.arg(fitting_method, c("direct_profile", "template_convolution"))
  model_pixel_sampling <- match.arg(model_pixel_sampling, c("point_sampled", "pixel_integrated"))
  kind <- if (model_pixel_sampling == "pixel_integrated") "pre" else "post"
  sigma <- lsf[[paste0("lsf_sigma_angstrom_", kind)]]
  if (is.null(sigma) || !any(is.finite(sigma) & sigma > 0)) stop("Selected wavelength-dependent LSF is unavailable.")
  list(sigma_angstrom = sigma, wavelength = lsf$wavelength,
    provenance = c(lsf$provenance, list(fitting_method = fitting_method,
      model_pixel_sampling = model_pixel_sampling, selected_lsf = kind,
      selected_native_extension = unname(lsf$provenance$native_extension[kind]),
      correction_applied = FALSE)))
}
