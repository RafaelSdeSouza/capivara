#' Cluster a 2D Representation of a Data Cube
#'
#' This function takes a 2D representation of a data cube (such as an IFU cube that has
#' been flattened into rows = spatial pixels and columns = spectral variables), scales
#' each row, computes pairwise distances using \code{\link{torch_dist}}, and then
#' performs hierarchical clustering. The resulting clusters are rearranged back
#' into a 2D grid consistent with the original spatial dimensions.
#'
#' @param cube_data A numeric matrix representing the 2D flattened data cube.
#'   Each row typically corresponds to a spatial pixel, and each column represents
#'   a spectral variable (e.g., wavelength).
#' @param n_rows Integer, the number of rows in the original cube’s spatial dimension.
#' @param n_cols Integer, the number of columns in the original cube’s spatial dimension.
#' @param n_clusters Integer, the number of clusters to form.
#' @param scale_fn A function used to scale each row of \code{cube_data}. Defaults to
#'   \code{\link[base]{scale}}. If you have a custom scaling function, pass it here.
#'
#' @details
#' Steps performed by the function:
#' \enumerate{
#'   \item Applies the scaling function \code{scale_fn} to each row of \code{cube_data}.
#'   \item Computes Manhattan distances between rows using \code{\link{torch_dist}}.
#'   \item Performs hierarchical clustering with Ward's D2 method via \code{\link[fastcluster]{hclust}}.
#'   \item Cuts the resulting dendrogram into \code{n_clusters} clusters.
#'   \item Reshapes the cluster assignments into a \code{n_rows x n_cols} matrix, restoring
#'   the original spatial layout.
#' }
#'
#' This process is often used in IFU data analysis, where clustering is applied to
#' grouped spectral profiles of spatial pixels to identify regions with similar
#' characteristics.
#'
#' @return A \code{n_rows x n_cols} matrix of cluster assignments (integers), where each
#'   element corresponds to a spatial pixel in the original cube layout.
#'
#' @import torch
#' @importFrom fastcluster hclust
#' @importFrom stats cutree
#'
#' @seealso \code{\link{torch_dist}}, \code{\link[fastcluster]{hclust}}, \code{\link[stats]{cutree}}
#'
#' @examples
#' if (torch::torch_is_installed()) {
#'   # Suppose we have a 10x10 spatial grid and a spectrum of length 50
#'   cube_data_example <- matrix(rnorm(10*10*50), nrow = 100, ncol = 50)
#'   cluster_map <- cube_cluster(
#'     cube_data = cube_data_example,
#'     n_rows = 10,
#'     n_cols = 10,
#'     n_clusters = 5,
#'     scale_fn = scale
#'   )
#'   dim(cluster_map) # Should be c(10, 10)
#' }
#'
#' @export
# Cluster the IFU cube data
cube_cluster <- function(input, Ncomp = 5, redshift = 0, scale_fn = median_scale) {
  # Step 1: Read the FITS cube
#  cubedat <- FITSio::readFITS(input, hdu = 1)
  cubedat <- input

  # Step 2: Extract the cube dimensions
  n_row <- dim(cubedat$imDat)[1]
  n_col <- dim(cubedat$imDat)[2]
  n_wave <- dim(cubedat$imDat)[3]

  # Step 3: Convert the cube to a 2D matrix using cube_to_matrix
  IFU2D <- cube_to_matrix(cubedat)

  # Step 4: Scale the data row-wise
  scaled_data <- t(apply(IFU2D, 1, scale_fn))

  # Step 5: Compute pairwise distances using torch_dist
  distance_matrix <- torch_dist(scaled_data)

  # Step 6: Perform hierarchical clustering using Ward.D2 method
  hc <- fastcluster::hclust(distance_matrix, method = "ward.D2")

  # Step 7: Cut the dendrogram into Ncomp clusters
  clusters <- cutree(hc, k = Ncomp)

  # Step 8: Reshape the cluster vector back to a 2D map
  cluster_map <- matrix(clusters, nrow = n_row, ncol = n_col)

  # Return the cluster map and cube metadata
  return(list(cluster_map = cluster_map, header = cubedat$hdr, axDat = cubedat$axDat))
}

