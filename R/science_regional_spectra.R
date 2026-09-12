#' Prepare native and resolution-controlled regional spectra
#'
#' Keeps the native finite-flux sum separate from a spectrum homogenized before
#' summing. This experimental stellar-only contract uses the Gaussian POST LSF
#' with a point-sampled model. It does not certify stellar populations. Input
#' variances assume independent native pixels and spaxels; induced spectral
#' covariance is retained. No spatial covariance correction is inferred.
#'
#' Summed regional spectra alone cannot supply the required spaxel responses.
#' Invalid fitting samples and their kernel halos are masked; the native sum
#' retains the original finite-flux convention. Requires Python numpy/scipy
#' through reticulate. No dependency is installed automatically.
#' @param segmentation A segmentation result containing `cluster_map`.
#' @param flux Native flux array in (x,y,wavelength) order.
#' @param variance Native variance array of the same dimensions, not IVAR.
#' @param lsf Result of `read_manga_lsf()`.
#' @param mask Logical fitting eligibility array, with the same dimensions.
#' @param galaxy_id Observation identifier.
#' @param segmentation_id Immutable partition identifier.
#' @param redshift Systemic redshift.
#' @param target_sigma Optional common target vector in observed vacuum Angstrom
#'   sigma, no narrower than any valid input. NULL selects each regional maximum.
#' @return A `capivara_regional_spectra` list with a serializable list of regions,
#'   native and fitting spectra, sparse spectral covariance and full provenance.
#' @export
prepare_segment_spectra <- function(segmentation, flux, variance, lsf, mask,
                                    galaxy_id, segmentation_id, redshift,
                                    target_sigma = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE) ||
      !requireNamespace("Matrix", quietly = TRUE)) {
    stop("Regional spectroscopy requires reticulate and Matrix.", call. = FALSE)
  }
  dims <- dim(flux)
  if (length(dims) != 3L || !identical(dim(variance), dims) ||
      !identical(dim(mask), dims) || !is.logical(mask) || anyNA(mask) ||
      !identical(dim(segmentation$cluster_map), dims[1:2])) {
    stop("Flux, variance, logical mask and partition dimensions must agree.")
  }
  selected <- select_manga_lsf(lsf, "template_convolution", "point_sampled")
  if (!identical(dim(selected$sigma_angstrom), dims) ||
      !identical(dim(lsf$lsf_sigma_angstrom_pre), dims)) stop("LSF dimensions differ from flux.")
  module <- reticulate::import_from_path("capivara_resolution",
    system.file("python", package = "capivara", mustWork = TRUE), convert = FALSE)
  ids <- sort(unique(as.vector(segmentation$cluster_map)))
  ids <- ids[is.finite(ids)]
  if (!length(ids) || any(ids <= 0) || any(ids != as.integer(ids))) stop("Positive integer segment IDs required.")
  flatten <- function(x) matrix(x, ncol = dims[3])
  f <- flatten(flux); v <- flatten(variance); m <- flatten(mask)
  pre <- flatten(lsf$lsf_sigma_angstrom_pre); post <- flatten(selected$sigma_angstrom)
  regions <- lapply(ids, function(id) {
    ix <- which(as.vector(segmentation$cluster_map) == id)
    result <- module$prepare_region(lsf$wavelength, f[ix,,drop=FALSE], v[ix,,drop=FALSE],
      pre[ix,,drop=FALSE], post[ix,,drop=FALSE], m[ix,,drop=FALSE], galaxy_id,
      as.integer(id), segmentation_id, redshift, lsf$provenance,
      strategy = if (is.null(target_sigma)) "local_max" else "common", target = target_sigma)
    reticulate::py_to_r(result)
  })
  names(regions) <- as.character(ids)
  structure(list(regions = regions, galaxy_id = galaxy_id,
    segmentation_id = segmentation_id, cluster_map = segmentation$cluster_map,
    package_version = as.character(utils::packageVersion("capivara")),
    scope = "Experimental regional response; stellar inference requires independent validation"),
    class = "capivara_regional_spectra")
}
