#' Cluster a 2D Representation of a Data Cube
#'
#' This function processes a data cube (such as an IFU cube) by flattening it into
#' rows (spatial pixels) and columns (spectral variables), scaling each row, computing
#' pairwise distances using \code{\link{torch_dist}}, and performing hierarchical
#' clustering. The resulting clusters are rearranged back into a 2D grid consistent
#' with the original spatial dimensions. The function also retains the original data cube
#' for reference and post-processing.
#'
#' Missing spectral channels are handled automatically: non-finite values produced
#' during row-wise scaling are replaced with zero before distances are computed.
#' This keeps the standard exact workflow usable for masked cubes without a
#' separate public entry point.
#'
#' @param input A FITS object representing the input data cube. Typically, this is an IFU data cube.
#' @param Ncomp Integer, the number of clusters to form.
#' @param redshift Numeric, the redshift to apply for wavelength correction. Defaults to 0 (no correction).
#' @param scale_fn A function used to scale each row of the 2D representation of the data cube.
#'   Defaults to \code{\link[base]{scale}}. If you have a custom scaling function, pass it here.
#' @param target_snr Optional minimum accepted SNR per cluster. When supplied,
#'   Capivara chooses the largest number of clusters whose minimum cluster SNR
#'   remains above this threshold.
#' @param var_cube Optional variance cube matching the input cube. Used only
#'   when \code{target_snr} is supplied.
#' @param k_values Optional candidate cluster counts tested when
#'   \code{target_snr} is supplied.
#' @param wavelength_range Optional wavelength interval used to compute SNR when
#'   \code{target_snr} is supplied.
#' @param snr_stat Either integrated SNR or median per-wavelength SNR when
#'   \code{target_snr} is supplied.
#' @param variance_inflation Multiplicative factor applied to propagated
#'   variances when \code{target_snr} is supplied.
#'
#' @details
#' Steps performed by the function:
#' \enumerate{
#'   \item Reads the input FITS data cube.
#'   \item Converts the cube into a 2D matrix (spatial pixels x spectral variables).
#'   \item Scales the data row-wise using \code{scale_fn}.
#'   \item Computes pairwise distances between rows using \code{\link{torch_dist}}.
#'   \item Performs hierarchical clustering using Ward's D2 method via \code{\link[fastcluster]{hclust}}.
#'   \item Cuts the dendrogram into \code{Ncomp} clusters and reshapes the results into a 2D cluster map,
#'   or, when \code{target_snr} is supplied, chooses the largest cut whose minimum cluster SNR
#'   remains above the requested threshold.
#'   \item Calculates the signal-to-noise ratio (SNR) for each cluster.
#' }
#'
#' @return A list containing:
#' \item{cluster_map}{A \code{n_rows x n_cols} matrix of cluster assignments (integers), where each element corresponds to a spatial pixel in the original cube layout.}
#' \item{header}{The header metadata from the input FITS file.}
#' \item{axDat}{Axis information (e.g., spatial and spectral axes) from the input FITS file.}
#' \item{cluster_snr}{A numeric vector containing the signal-to-noise ratio (SNR) for each cluster.}
#' \item{original_cube}{The original FITS data cube as input to the function, for reference and post-processing.}
#'
#' This process is often used in IFU data analysis, where clustering is applied to grouped
#' spectral profiles of spatial pixels to identify regions with similar characteristics.
#'
#' @seealso \code{\link{segment_big_cube}}, \code{\link{segment_starlet}},
#'   \code{\link{torch_dist}}, \code{\link[fastcluster]{hclust}},
#'   \code{\link[stats]{cutree}}
#'
#' @examples
#' if (torch::torch_is_installed()) {
#'   # Read a FITS cube and cluster it
#'   input_cube <- FITSio::readFITS("manga_7443_12703_LOGCUBE.fits")
#'   clustering_result <- segment(input = input_cube, Ncomp = 5)
#'
#'   # Access the cluster map and original cube
#'   cluster_map <- clustering_result$cluster_map
#'   original_cube <- clustering_result$original_cube
#' }
#'
#' @export
segment <- function(input,
                    Ncomp = 5,
                    redshift = 0,
                    scale_fn = median_scale,
                    target_snr = NULL,
                    var_cube = NULL,
                    k_values = NULL,
                    wavelength_range = NULL,
                    snr_stat = c("integrated", "median_per_wavelength"),
                    variance_inflation = 1) {
  if (!is.null(target_snr)) {
    if (!missing(Ncomp)) {
      stop("Specify either `Ncomp` or `target_snr`, not both.")
    }

    return(choose_ncomp_by_snr(
      input = input,
      target_snr = target_snr,
      var_cube = var_cube,
      k_values = k_values,
      wavelength_range = wavelength_range,
      redshift = redshift,
      scale_fn = scale_fn,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation,
      na_safe = TRUE
    ))
  }

  .segment_core(
    input = input,
    Ncomp = Ncomp,
    redshift = redshift,
    scale_fn = scale_fn,
    na_to_zero = TRUE
  )
}
