#' Compute Median Spectra and Process Clustered Data from a Data Cube
#'
#' This function processes the result of a clustering procedure to:
#' \enumerate{
#'   \item Compute the median spectrum for each cluster,
#'   \item Generate a new 3D array (or "cube") with the median spectrum assigned to each pixel, and
#'   \item Create a table linking each pixel's original spectrum to its assigned cluster.
#' }
#'
#' @param cluster_result A list containing:
#'   \itemize{
#'     \item \code{original_cube}: A FITS object (or similar) representing the input data cube.
#'     \item \code{cluster_map}: A matrix or array indicating each pixel's cluster assignment.
#'     \item \code{axDat}, \code{header}, \code{cluster_snr}: Additional metadata (if available) from the clustering process.
#'   }
#'
#' @details
#' **Steps performed by the function**:
#' \enumerate{
#'   \item Reads the input data cube (from \code{cluster_result$original_cube}) and the cluster map.
#'   \item Reshapes the data cube into a 2D matrix, aligning each pixel with its spectrum.
#'   \item Computes the median spectrum for each cluster by grouping pixels based on their cluster assignment.
#'   \item Creates a new 3D array where each pixel's spectrum is replaced by the median spectrum of its cluster.
#'   \item Generates a table linking each pixel's original spectrum to its assigned cluster.
#' }
#'
#' This legacy helper is retained for compatibility. For new workflows, prefer
#' \code{\link{reconstruct_cluster_cube}} or
#' \code{\link{reconstruct_flux_preserving_cube}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{MeanSpec}}{A data frame with the median spectrum for each cluster (rows = clusters, columns = wavelength).}
#'   \item{\code{reshaped_cube}}{A 3D array containing the new cube with median spectra assigned to each pixel.}
#'   \item{\code{post_process_table}}{A data frame linking each pixel's original spectrum to its assigned cluster.}
#'   \item{\code{cluster_summary}}{A richer summary object returned by \code{\link{summarize_cluster_spectra}}.}
#' }
#'
#' @examples
#' # Example clustering result (pseudo-code)
#' cluster_result <- list(
#'   original_cube = FITSio::readFITS("manga-11749-12701-LOGCUBE.fits"),
#'   cluster_map = some_cluster_map,
#'   axDat = some_axDat
#' )
#'
#' # Run SpecMean
#' result <- SpecMean(cluster_result)
#'
#' # Access outputs
#' median_specs       <- result$MeanSpec
#' median_spectrum_cube <- result$reshaped_cube
#' pixel_cluster_table  <- result$post_process_table
#'
#' @import FITSio
#' @importFrom dplyr group_by summarise mutate ungroup
#' @importFrom reshape2 melt dcast
#' @importFrom utils write.csv
#' @export
SpecMean <- function(cluster_result) {
  cubedat <- .as_cubedat(cluster_result$original_cube)
  cluster_map <- cluster_result$cluster_map
  summary <- summarize_cluster_spectra(cluster_result)
  wavelength <- summary$wavelength
  recon <- reconstruct_cluster_cube(
    cluster_result = cluster_result,
    template = "median",
    preserve_mask = FALSE,
    fill_mode = "na",
    return_residual = FALSE
  )

  IFU2D <- cube_to_matrix(cubedat)
  class2D <- as.vector(cluster_map)

  if (length(class2D) != nrow(IFU2D)) {
    stop("Cluster map dimensions do not match the input cube dimensions.")
  }

  MeanSpec <- data.frame(
    class = summary$cluster_ids,
    summary$median_spectra,
    check.names = FALSE
  )

  reshaped_cube <- recon$model_cube$imDat

  post_process_table <- as.data.frame(IFU2D)
  post_process_table$class <- class2D
  colnames(post_process_table) <- c(as.character(wavelength), "class")

  return(list(
    MeanSpec           = MeanSpec,
    reshaped_cube      = reshaped_cube,
    post_process_table = post_process_table,
    cluster_summary    = summary
  ))
}
