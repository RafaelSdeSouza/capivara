#' Compute Median Spectra and Process Clustered Data from a Data Cube
#'
#' This function processes a data cube and the corresponding clustering output to compute
#' the median spectrum for each cluster, generate a new cube with median spectra assigned
#' to each pixel, and create a table linking each pixel's original spectrum to its cluster.
#'
#' @param input_cube A FITS file or object representing the input data cube.
#'   Typically, this is the original IFU data cube used in clustering.
#' @param cluster_output A list returned by the \code{\link{cube_cluster}} function,
#'   containing the cluster map (\code{cluster_map}), metadata (\code{header}, \code{axDat}),
#'   and SNR estimates (\code{cluster_snr}).
#' @param output_dir A character string specifying the directory where the output files
#'   will be saved. Defaults to \code{"output_results"}.
#'
#' @details
#' Steps performed by the function:
#' \enumerate{
#'   \item Reads the input data cube and clustering output.
#'   \item Reshapes the data cube and cluster map into 2D formats.
#'   \item Computes the median spectrum for each cluster by grouping pixels based on their cluster assignment.
#'   \item Creates a new cube where each pixel's spectrum is replaced by the median spectrum of its cluster.
#'   \item Generates a table linking each pixel's original spectrum to its assigned cluster.
#'   \item Saves the median spectrum cube and the table to the specified directory.
#' }
#'
#' This function is designed to facilitate post-processing and analysis of clustered
#' spectral data, commonly used in IFU data analysis to study regions with similar
#' spectral properties.
#'
#' @return A list containing:
#' \item{MeanSpec}{A data frame with the median spectrum for each cluster, where rows correspond to clusters and columns correspond to wavelengths.}
#' \item{median_cube_file}{The file path to the FITS file containing the new cube with median spectra assigned to each pixel.}
#' \item{post_process_table}{The file path to the CSV file containing the table of original spectra and their associated clusters.}
#'
#' @seealso \code{\link{cube_cluster}}, \code{\link{FITSio}}, \code{\link{reshape2}}, \code{\link{dplyr}}
#'
#' @examples
#' # Example usage with a FITS file and cluster output
#' cube_cluster_result <- cube_cluster(
#'   input = FITSio::readFITS("manga-11749-12701-LOGCUBE.fits"),
#'   Ncomp = 5
#' )
#'
#' Spec_mean_result <- Spec_mean(
#'   input_cube = FITSio::readFITS("manga-11749-12701-LOGCUBE.fits"),
#'   cluster_output = cube_cluster_result
#' )
#'
#' # Access outputs
#' MeanSpec <- Spec_mean_result$MeanSpec
#' median_cube_file <- Spec_mean_result$median_cube_file
#' post_process_table <- Spec_mean_result$post_process_table
#'
#' @import FITSio
#' @importFrom dplyr group_by summarise mutate ungroup
#' @importFrom reshape2 melt dcast
#' @importFrom utils write.csv
#'
#' @export
# Compute Median Spectra and Process Clustered Data
Spec_mean <- function(cluster_result) {
  # Step 1: Extract components from the cluster result
  cubedat <- cluster_result$original_cube
  cluster_map <- cluster_result$cluster_map
  wavelength <- FITSio::axVec(3, cubedat$axDat)  # Wavelength axis

  # Step 2: Convert the original cube to a 2D matrix using cube_to_matrix
  IFU2D <- cube_to_matrix(cubedat)

  # Step 3: Flatten the cluster map
  class2D <- as.vector(cluster_map)

  # Ensure dimensions of class2D match IFU2D
  if (length(class2D) != nrow(IFU2D)) {
    stop("Cluster map dimensions do not match the input cube dimensions.")
  }

  # Step 4: Prepare data frame with spectra and their corresponding clusters
  IFU2Db <- as.data.frame(IFU2D)
  IFU2Db$class <- class2D
  colnames(IFU2Db) <- c(wavelength, "class")

  # Step 5: Compute the median spectrum for each cluster
  MeanSpec <- IFU2Db %>%
    reshape2::melt(id = "class") %>%
    dplyr::group_by(class, variable) %>%
    dplyr::summarise(MedianSpec = median(value, na.rm = TRUE), .groups = "drop") %>%
    reshape2::dcast(class ~ variable, value.var = "MedianSpec")

  # Step 6: Create a new cube with median spectra
  median_cube <- matrix(NA, nrow = nrow(IFU2D), ncol = ncol(IFU2D))
  unique_groups <- unique(class2D)

  for (group in unique_groups) {
    if (!is.na(group)) {
      # Extract the median spectrum for the current group
      median_spectrum <- MeanSpec %>%
        dplyr::filter(class == group) %>%
        dplyr::select(-class) %>%
        unlist() %>%
        as.numeric()

      # Assign the median spectrum to the corresponding rows in the cube
      group_indices <- which(class2D == group)
      median_cube[group_indices, ] <- matrix(median_spectrum, nrow = length(group_indices), ncol = ncol(IFU2D), byrow = TRUE)
    }
  }
  reshaped_cube <- array(median_cube, dim = dim(cubedat$imDat))

  # Step 7: Generate a table linking each pixel's spectrum to its cluster group
  post_process_table <- as.data.frame(IFU2D)
  post_process_table$class <- class2D
  colnames(post_process_table) <- c(wavelength, "class")

  # Step 8: Return structured outputs
  return(list(
    MeanSpec = MeanSpec,
    reshaped_cube = reshaped_cube,
    post_process_table = post_process_table
  ))
}
