#' Segment spectra using a declared scientific representation
#'
#' Use the observed-entry spatial Ward hierarchy. Unlike the historical
#' interfaces, this entry point requires an explicit representation and
#' requested region count. Unsupported merges can leave more regions.
#' @param cube FITS-like cube or numeric three-dimensional array.
#' @param n_regions Positive integer requested hierarchy cut.
#' @param representation A non-null representation declaration, such as
#'   [full_flux_representation()] or [manga_sandra_spectral_shape()].
#' @param variance Optional aligned variance cube, not inverse variance.
#' @param sample_validity Optional aligned logical voxel eligibility array.
#' @param support Optional explicit spatial support object.
#' @param redshift Systemic redshift, required where rest-frame selection needs it.
#' @param wavelength_frame Native wavelength frame, if not declared in cube metadata.
#' @return The unchanged CAPIVARA segmentation object returned by [segment_large()].
#' @export
segment_spectra <- function(cube, n_regions, representation,
                            variance = NULL, sample_validity = NULL,
                            support = NULL, redshift = NA_real_,
                            wavelength_frame = NULL) {
  if (missing(n_regions) || !is.numeric(n_regions) || length(n_regions) != 1L ||
      !is.finite(n_regions) || n_regions < 1 || n_regions != floor(n_regions) ||
      n_regions > .Machine$integer.max)
    stop("Supply a positive integer n_regions.", call. = FALSE)
  if (missing(representation) || is.null(representation))
    stop("Choose an explicit representation; see full_flux_representation().",
         call. = FALSE)
  segment_large(input = cube, Ncomp = as.integer(n_regions),
                representation = representation, var_cube = variance,
                sample_validity = sample_validity, support = support,
                redshift = redshift, wavelength_frame = wavelength_frame)
}

#' Extract and visualize regional spectral products
#'
#' These convenience interfaces retain the existing algorithms, return
#' structures and defaults. Spectral summaries are sums over finite
#' contributors: inspect finite_counts before treating them as fixed-aperture
#' spectra. Variance propagation is diagonal, with an optional scalar
#' inflation, not a full covariance calculation.
#' @param segmentation A CAPIVARA segmentation result.
#' @param variance Optional aligned variance cube.
#' @param variance_inflation Positive variance multiplier, default one.
#' @param ... Arguments to the corresponding existing plotting function.
#' @return For extraction, the original spectral summary list; for plotting,
#'   the original ggplot object.
#' @export
extract_region_spectra <- function(segmentation, variance = NULL,
                                   variance_inflation = 1) {
  summarize_cluster_spectra(cluster_result = segmentation, var_cube = variance,
                            variance_inflation = variance_inflation)
}

#' @rdname extract_region_spectra
#' @export
plot_regions <- function(segmentation, ...) {
  plot_cluster(cluster_data = segmentation, ...)
}

#' @rdname extract_region_spectra
#' @export
plot_region_spectra <- function(segmentation, ...) {
  plot_cluster_spectra(cluster_result = segmentation, ...)
}

#' Inspect region geometry and the optical line registry
#'
#' @param segmentation A CAPIVARA segmentation result.
#' @param coords Coordinate convention passed to [regions_sf()].
#' @param wavelength_medium Laboratory wavelength convention.
#' @return An sf geometry object for region_geometry, or a line registry
#'   data frame for list_emission_lines.
#' @export
region_geometry <- function(segmentation,
                            coords = c("auto", "native", "sky", "physical")) {
  regions_sf(x = segmentation, coords = match.arg(coords))
}

#' @rdname region_geometry
#' @export
list_emission_lines <- function(wavelength_medium = c("vacuum", "air")) {
  emission_lines(wavelength_medium = match.arg(wavelength_medium))
}
