#' Compute Median Spectra and Process Clustered Data from a Data Cube
#'
#' This function processes an input data cube and its corresponding clustering output to compute
#' the median spectrum for each cluster, generate a new cube where each pixel's spectrum is replaced
#' by the median spectrum of its cluster, and create a table linking each pixel's original spectrum to its cluster.
#'
#' @param input_cube A FITS file or object representing the input data cube.
#' @param cluster_output A list returned by the \code{\link{cube_cluster}} function, containing the
#'   cluster map (\code{cluster_map}), metadata (\code{header}, \code{axDat}), and SNR estimates (\code{cluster_snr}).
#' @param output_dir A character string specifying the directory where the output files will be saved.
#'
#' @return A list containing:
#' \item{median_spectra}{A data frame with the median spectrum for each cluster, with rows corresponding to clusters and columns to wavelengths.}
#' \item{median_cube}{The FITS cube where each pixel's spectrum is replaced by its cluster's median spectrum.}
#' \item{cluster_assignment_table}{A CSV file linking each pixel's original spectrum to its assigned cluster.}
#'
#' @export
compute_cluster_median_spectra <- function(input_cube, cluster_output, output_dir = "output_results") {
  # Extract components from the cluster output
  cube_data <- cluster_output$original_cube
  cluster_map <- cluster_output$cluster_map
  wavelength <- FITSio::axVec(3, cube_data$axDat)  # Wavelength axis

  # Convert the original cube to a 2D matrix
  spectral_matrix <- cube_to_matrix(cube_data)

  # Flatten the cluster map to get cluster assignments
  cluster_assignments <- as.vector(cluster_map)

  if (length(cluster_assignments) != nrow(spectral_matrix)) {
    stop("Cluster map dimensions do not match the input cube dimensions.")
  }

  # Prepare data frame with spectra and their corresponding clusters
  spectra_df <- as.data.frame(spectral_matrix)
  spectra_df$cluster <- cluster_assignments
  colnames(spectra_df) <- c(wavelength, "cluster")

  # Compute the median spectrum for each cluster
  median_spectra <- spectra_df %>%
    reshape2::melt(id = "cluster") %>%
    dplyr::group_by(cluster, variable) %>%
    dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = "drop") %>%
    reshape2::dcast(cluster ~ variable, value.var = "median_value")

  # Create a new cube with median spectra for each cluster
  median_cube <- matrix(NA, nrow = nrow(spectral_matrix), ncol = ncol(spectral_matrix))
  unique_groups <- unique(cluster_assignments)

  for (group in unique_groups) {
    if (!is.na(group)) {
      median_spectrum <- median_spectra %>%
        dplyr::filter(cluster == group) %>%
        dplyr::select(-cluster) %>%
        unlist() %>%
        as.numeric()

      group_indices <- which(cluster_assignments == group)
      median_cube[group_indices, ] <- matrix(median_spectrum, nrow = length(group_indices), ncol = ncol(spectral_matrix), byrow = TRUE)
    }
  }
  reshaped_cube <- array(median_cube, dim = dim(cube_data$imDat))

  # Generate a table linking each pixel's spectrum to its cluster
  cluster_assignment_table <- as.data.frame(spectral_matrix)
  cluster_assignment_table$cluster <- cluster_assignments
  colnames(cluster_assignment_table) <- c(wavelength, "cluster")

  # Save outputs as needed using output_dir, e.g., using write.csv or FITSio::writeFITS

  return(list(
    median_spectra = median_spectra,
    median_cube = reshaped_cube,
    cluster_assignment_table = cluster_assignment_table
  ))
}

