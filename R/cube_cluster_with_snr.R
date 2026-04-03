#' Legacy Target-SNR Wrapper for Capivara Segmentation
#'
#' This function is kept for backward compatibility and now forwards to
#' \code{\link{choose_ncomp_by_snr}} using the original flux-only SNR
#' approximation. For variance-aware SNR cuts, use
#' \code{\link{choose_ncomp_by_snr}} directly.
#'
#' @param input A 3D array representing the IFU data cube (dimensions: rows, columns, wavelengths).
#'   Alternatively, this can be a `readFITS` object containing the data cube and its metadata.
#' @param target_snr Numeric. The desired signal-to-noise ratio (SNR) for the final bins.
#' @param redshift Numeric. The redshift value used for any wavelength corrections. Default is 0.
#' @param scale_fn A function to scale each spaxel spectrum. Defaults to `median_scale`,
#'   which centers each spectrum by subtracting its median value.
#' @param sn_func A function to compute the signal-to-noise ratio for a group of spaxels.
#'   Defaults to a standard implementation that sums the signal and noise across spaxels.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{cluster_map}}{A matrix of cluster assignments (dimensions: rows x columns).}
#'   \item{\code{header}}{The FITS header associated with the input data cube.}
#'   \item{\code{axDat}}{The axis metadata from the input data cube.}
#' }
#'
#' @importFrom FITSio readFITS
#' @importFrom fastcluster hclust
#' @importFrom stats cutree
#'
#' @seealso \code{\link{choose_ncomp_by_snr}}, \code{\link{segment}},
#'   \code{\link{segment}}
#'
#' @examples
#' if (torch::torch_is_installed()) {
#'   # Read the FITS cube and cluster based on a target SNR
#'   fits_data <- FITSio::readFITS("example_cube.fits")
#'   result <- cube_cluster_with_snr(
#'     input = fits_data,
#'     target_snr = 50,
#'     redshift = 0,
#'     scale_fn = median_scale
#'   )
#'   cluster_map <- result$cluster_map
#'   dim(cluster_map) # Dimensions match the spatial size of the FITS cube
#' }
#'
#' @export

cube_cluster_with_snr <- function(input, target_snr, redshift = 0, scale_fn = median_scale) {
  choose_ncomp_by_snr(
    input = input,
    target_snr = target_snr,
    redshift = redshift,
    scale_fn = scale_fn,
    na_safe = FALSE
  )
}
