#' Cluster a 2D Representation of an IFU Data Cube (Robust Version)
#'
#' This function performs unsupervised clustering on a spectral data cube
#' (e.g., IFU data) using a *robust preprocessing pipeline* designed to handle
#' missing values, masked pixels, negative fluxes, and low-variance spectra.
#'
#' The cube is flattened into a 2D matrix (spatial pixels × spectral channels),
#' filtered to remove invalid or uninformative spectra, scaled using a robust
#' per-row transformation, and clustered via Ward’s hierarchical method using
#' pairwise distances computed by \code{\link{torch_dist}}.
#'
#' @param input A FITS-like list object containing:
#'   \describe{
#'     \item{\code{imDat}}{A 3D numeric array \code{[n_row, n_col, n_wave]} representing
#'       the IFU data cube.}
#'     \item{\code{hdr}}{Header metadata (e.g., FITS header).}
#'     \item{\code{axDat}}{Axis information such as wavelength and spatial coordinates.}
#'   }
#' @param Ncomp Integer, the number of clusters to generate.
#' @param redshift Numeric redshift to apply for wavelength correction
#'   (currently unused but reserved for future extensions). Defaults to 0.
#' @param scale_fn Deprecated — kept for backward compatibility. The robust
#'   version implements its own safe scaling internally and ignores this argument.
#'
#' @details
#' The robust segmentation pipeline applies the following steps:
#'
#' \enumerate{
#'
#'   \item \textbf{Reshape the cube} into a 2D matrix
#'     \code{IFU2D = pixels × wavelengths}.
#'
#'   \item \textbf{Validity filtering of spectra} using three complementary criteria:
#'     \itemize{
#'       \item Fraction of finite spectral values must exceed a threshold
#'         (default: \eqn{\ge 80\%} or at least 10 finite samples).
#'       \item Spectral scatter (via MAD) must be non-zero, preventing constant or
#'         flat spectra from being clustered.
#'       \item Energy (\eqn{\sum v^2}) must be positive and finite, ensuring the
#'         pixel carries signal beyond pure noise or NaNs.
#'     }
#'     Only pixels passing \emph{all} filters are clustered.
#'
#'   \item \textbf{Robust row-wise scaling} using a custom safe scaling function:
#'     \itemize{
#'       \item Centering with the median of finite values.
#'       \item Scaling with MAD; if MAD = 0, fallback to SD; if SD = 0, fallback to zeros.
#'       \item All non-finite results are replaced with 0.
#'     }
#'     This guarantees a fully finite matrix for distance computation.
#'
#'   \item \textbf{Distance computation} using \code{\link{torch_dist}}
#'     to generate pairwise L1 distances in a memory-efficient way.
#'
#'   \item \textbf{Hierarchical clustering} with Ward.D2 linkage
#'     via \code{\link[fastcluster]{hclust}}.
#'
#'   \item \textbf{Reprojection of clusters} into a 2D spatial
#'     cluster map \code{n_row × n_col}, with invalid pixels filled as \code{NA}.
#'
#'   \item \textbf{Cluster SNR estimation} using a robust definition avoiding
#'     negative flux or zero-noise instabilities:
#'     \deqn{ \mathrm{SNR} = \frac{\sum S}{\sqrt{\sum N^2}} }
#'     where invalid or zero-noise pixels are safely handled.
#'
#' }
#'
#' This robust workflow avoids failures common in IFU segmentation:
#' \itemize{
#'   \item NaN / Inf contamination,
#'   \item flat or zero-variance spectra,
#'   \item negative fluxes after sky subtraction,
#'   \item distance matrix instabilities,
#'   \item meaningless clusters from blank spaxels.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{cluster_map}}{An integer matrix \code{[n_row × n_col]}
#'         assigning each valid pixel to a cluster; invalid pixels are \code{NA}.}
#'   \item{\code{header}}{The original cube header.}
#'   \item{\code{axDat}}{Axis metadata (e.g., wavelengths, coordinates).}
#'   \item{\code{cluster_snr}}{Numerical SNR estimate for each cluster.}
#'   \item{\code{original_cube}}{The unmodified input cube.}
#' }
#'
#' @seealso
#'   \code{\link{torch_dist}},
#'   \code{\link[fastcluster]{hclust}},
#'   \code{\link[stats]{cutree}}
#'
#' @examples
#' \dontrun{
#'   input_cube <- FITSio::readFITS("manga_cube.fits")
#'   result <- segment(input_cube, Ncomp = 5)
#'   image(result$cluster_map)
#' }
#'
#' @export
# Cluster the IFU cube data
segment_robust <- function(input, Ncomp = 5, redshift = 0, scale_fn = median_scale) {
  # Step 1: Read the FITS cube
  cubedat <- input

  # Step 2: Extract the cube dimensions
  n_row <- dim(cubedat$imDat)[1]
  n_col <- dim(cubedat$imDat)[2]
  n_wave <- dim(cubedat$imDat)[3]

  # Step 3: Convert the cube to a 2D matrix
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

  # Step 8: Perform hierarchical clustering using Ward.D2 method
  hc <- fastcluster::hclust(distance_matrix, method = "ward.D2")

  # Step 9: Cut the dendrogram into Ncomp clusters
  clusters <- cutree(hc, k = Ncomp)

  # Step 10: Reshape the cluster vector back to a 2D map
  cluster_map <- matrix(NA, nrow = n_row, ncol = n_col)
  cluster_map[valid_indices] <- clusters

  # Step 11: Calculate SNR for each cluster
  sn_func <- function(indices) {
    sum(signal_valid[indices]) / sqrt(sum(noise_valid[indices]^2))
  }
  cluster_snr <- sapply(unique(clusters), function(cl) {
    indices <- which(clusters == cl)
    sn_func(indices)
  })

  # Return the cluster map, cube metadata, SNR estimates, and original cube
  return(list(
    cluster_map = cluster_map,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    cluster_snr = cluster_snr,
    original_cube = cubedat
  ))
}
