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
#' @return A list containing:
#' \describe{
#'   \item{\code{MeanSpec}}{A data frame with the median spectrum for each cluster (rows = clusters, columns = wavelength).}
#'   \item{\code{reshaped_cube}}{A 3D array containing the new cube with median spectra assigned to each pixel.}
#'   \item{\code{post_process_table}}{A data frame linking each pixel's original spectrum to its assigned cluster.}
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
  # Step 1: Extract components from the cluster result
  cubedat     <- cluster_result$original_cube
  cluster_map <- cluster_result$cluster_map
  wavelength  <- FITSio::axVec(3, cubedat$axDat)  # Wavelength axis

  # Step 2: Convert the original cube to a 2D matrix
  IFU2D <- cube_to_matrix(cubedat)

  # Step 3: Flatten the cluster map
  class2D <- as.vector(cluster_map)

  # Ensure dimensions match
  if (length(class2D) != nrow(IFU2D)) {
    stop("Cluster map dimensions do not match the input cube dimensions.")
  }

  # Step 4: Prepare data frame with spectra + cluster assignments
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
      # Extract the median spectrum for this cluster
      median_spectrum <- MeanSpec %>%
        dplyr::filter(class == group) %>%
        dplyr::select(-class) %>%
        unlist() %>%
        as.numeric()

      # Assign to matching pixels
      group_indices <- which(class2D == group)
      median_cube[group_indices, ] <- matrix(
        median_spectrum,
        nrow = length(group_indices),
        ncol = ncol(IFU2D),
        byrow = TRUE
      )
    }
  }

  reshaped_cube <- array(median_cube, dim = dim(cubedat$imDat))

  # Step 7: Table linking each pixel's original spectrum to its cluster
  post_process_table <- as.data.frame(IFU2D)
  post_process_table$class <- class2D
  colnames(post_process_table) <- c(wavelength, "class")

  # Step 8: Return outputs
  return(list(
    MeanSpec           = MeanSpec,
    reshaped_cube      = reshaped_cube,
    post_process_table = post_process_table
  ))
}

