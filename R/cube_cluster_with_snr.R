#' Cluster a 2D Representation of an IFU Cube Based on Target SNR
#'
#' This function takes an IFU data cube, computes the signal-to-noise ratio (SNR)
#' for each spaxel, and performs clustering using hierarchical methods to group
#' spaxels into regions with similar spectral properties, ensuring that the bins
#' achieve the desired target SNR.
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
#' @details
#' Steps performed by the function:
#' \enumerate{
#'   \item Converts the 3D data cube into a 2D matrix where rows correspond to spatial pixels
#'         and columns to spectral channels.
#'   \item Calculates the signal (sum of flux across wavelengths) and noise (square root of the signal).
#'   \item Initializes bins starting from the spaxel with the highest SNR and iteratively
#'         adds neighboring spaxels until the target SNR is met for each bin.
#'   \item Performs hierarchical clustering with Ward's D2 method on the scaled spectral matrix.
#'   \item Reshapes the resulting cluster assignments into a 2D map consistent with the original
#'         spatial dimensions of the data cube.
#' }
#'
#' This approach is suitable for analyzing IFU data to identify regions with
#' similar spectral characteristics while meeting SNR constraints for robust analysis.
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
#' @import torch
#'
#' @seealso \code{\link[FITSio]{readFITS}}, \code{\link[fastcluster]{hclust}},
#'   \code{\link[stats]{cutree}}, \code{\link{median_scale}}
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
  # Step 1: Read the FITS cube
  cubedat <- input

  # Step 2: Extract the cube dimensions
  n_row <- dim(cubedat$imDat)[1]
  n_col <- dim(cubedat$imDat)[2]
  n_wave <- dim(cubedat$imDat)[3]

  # Step 3: Convert the cube to a 2D matrix using cube_to_matrix
  IFU2D <- cube_to_matrix(cubedat)

  # Step 4: Calculate signal and noise
  signal <- rowSums(IFU2D, na.rm = TRUE)  # Total signal per spatial pixel
  noise <- sqrt(signal)  # Assuming Poisson noise

  # Handle missing or invalid data
  signal[is.na(signal) | signal <= 0] <- 0
  noise[is.na(noise) | noise == 0] <- Inf  # Avoid division errors

  # Step 5: Filter valid pixels (based on signal > 0)
  valid_indices <- which(signal > 0)
  if (length(valid_indices) == 0) stop("No valid pixels after filtering.")

  IFU2D_valid <- IFU2D[valid_indices, ]
  signal_valid <- signal[valid_indices]
  noise_valid <- noise[valid_indices]

  # Step 6: Scale the valid data row-wise
  scaled_data <- t(apply(IFU2D_valid, 1, scale_fn))

  # Step 7: Compute pairwise distances using torch_dist
  distance_matrix <- torch_dist(scaled_data)

  # Step 8: Perform hierarchical clustering with Ward.D2 method
  hc <- fastcluster::hclust(distance_matrix, method = "ward.D2")

  # Step 9: Define the number of groups to meet the SNR target, searching backwards
  sn_func <- function(indices) {
    sum(signal_valid[indices]) / sqrt(sum(noise_valid[indices]^2))
  }

  # Start with the maximum number of clusters and work backwards
  Nmax <- length(valid_indices)
  Ncomp <- NA  # To store the final number of components
  cluster_snr <- NA  # To store SNR of the clusters

  for (k in seq(Nmax, 1)) {  # Reverse search
    clusters <- cutree(hc, k = k)

    # Efficiently calculate SNR for each cluster
    cluster_snr <- vapply(unique(clusters), function(cl) {
      indices <- which(clusters == cl)
      if (length(indices) > 0) {
        sn_func(indices)
      } else {
        NA  # Return NA for empty clusters
      }
    }, numeric(1))

    # Check if all clusters meet the target SNR
    if (all(cluster_snr >= target_snr, na.rm = TRUE)) {
      Ncomp <- k
      break
    }
  }

  if (is.na(Ncomp)) {
    stop("No clustering configuration satisfies the target SNR.")
  }

  # Step 10: Map bin IDs back to the 2D grid
  clusters <- cutree(hc, k = Ncomp)
  cluster_map <- matrix(NA, nrow = n_row, ncol = n_col)
  cluster_map[valid_indices] <- clusters

  # Return the cluster map and metadata
  return(list(
    cluster_map = cluster_map,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    Ncomp = Ncomp,
    cluster_snr = cluster_snr
  ))
}
