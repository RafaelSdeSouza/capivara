#' Compute Pairwise Distances Using Torch
#'
#' This function computes pairwise distances between observations in a numeric matrix using
#' the `torch` package's \code{\link[torch]{nnf_pdist}} function. By default, it
#' computes L1 (Manhattan) distances. The result is returned as a \code{dist}
#' object, making it directly compatible with various clustering and multivariate
#' analysis functions in R.
#'
#' @param x A numeric matrix. Rows correspond to observations and columns correspond
#'   to variables.
#'
#' @details
#' This function:
#' \enumerate{
#'   \item Ensures \code{x} is a numeric matrix.
#'   \item Converts \code{x} into a torch tensor.
#'   \item Uses \code{\link[torch]{nnf_pdist}} with \code{p = 1}, resulting in Manhattan distances.
#'   \item Fills a symmetric distance matrix using these pairwise distances, and returns
#'   the result via \code{\link[stats]{as.dist}}.
#' }
#'
#' If you have a GPU and want to leverage it, you can modify the code to place
#' \code{x_ten} on a CUDA device. Ensure \code{torch::cuda_is_available()} returns
#' \code{TRUE} before doing so.
#'
#' @return A \code{dist} object representing the pairwise Manhattan distances among rows of \code{x}.
#'
#' @import torch
#' @seealso \code{\link[stats]{dist}}, \code{\link[torch]{nnf_pdist}}
#'
#' @examples
#' if (torch::torch_is_installed()) {
#'   x <- matrix(rnorm(10 * 3), nrow = 10, ncol = 3)
#'   d <- torch_dist(x)
#'   print(d)
#'   # Use the resulting dist object with clustering functions, e.g.:
#'   # hclust(d)
#' }
#'
#' @export
torch_dist <- function(x) {
  # Ensure x is a numeric matrix
  x <- as.matrix(x)

  N <- nrow(x)

  # Initialize an NxN matrix for the distances
  mat <- matrix(0, nrow = N, ncol = N)

  # If GPU is desired and available, you could do:
  # device <- if (torch::cuda_is_available()) torch::torch_device("cuda:0") else torch::torch_device("cpu")
  # x_ten <- torch::torch_tensor(x, dtype = torch::torch_float(), device = device)

  # For simplicity, we stick to CPU here
  x_ten <- torch::torch_tensor(x, dtype = torch::torch_float())

  # Compute pairwise distances (Manhattan distance)
  pd <- torch::nnf_pdist(x_ten, p = 1)

  # Convert the torch tensor to an R numeric vector
  pd <- as.matrix(pd)

  # Fill the lower triangle of mat with distances
  # pd is given in order: (1,2), (1,3), ..., (1,N), (2,3), ..., (N-1,N)
  # This matches the order as.dist() expects from the lower triangle of a distance matrix.
  mat[lower.tri(mat, diag = FALSE)] <- pd
  # Return as a dist object
  as.dist(mat)
}
