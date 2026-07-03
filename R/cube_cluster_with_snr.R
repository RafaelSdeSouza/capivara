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
#' @param redshift Numeric redshift placeholder kept for API compatibility.
#' @param scale_fn A function to scale each spaxel spectrum. Defaults to `median_scale`,
#'   which centers each spectrum by subtracting its median value.
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
#' @seealso \code{\link{choose_ncomp_by_snr}}, \code{\link{segment}}
#'
#' @examples
#' cube <- list(imDat = array(runif(4 * 4 * 8, min = 1), dim = c(4, 4, 8)))
#' result <- cube_cluster_with_snr(
#'   input = cube,
#'   target_snr = 1,
#'   redshift = 0,
#'   scale_fn = median_scale
#' )
#' cluster_map <- result$cluster_map
#' dim(cluster_map)
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
