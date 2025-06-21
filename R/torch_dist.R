#' Compute Pairwise Distances Using Torch
#'
#' This function computes pairwise distances between rows of a numeric matrix using
#' the `torch` backend. It leverages the highly efficient `nnf_pdist()` function,
#' optionally using GPU if available.
#'
#' @param x A numeric matrix. Rows are observations, columns are features.
#' @param p The order of the norm (default is 1 for Manhattan; use 2 for Euclidean).
#' @param device A torch device, e.g., `"cpu"` or `"cuda"`. If NULL, selects `"cuda"`
#'   if available, otherwise `"cpu"`.
#'
#' @return An object of class `dist`, compatible with R clustering and multivariate analysis functions.
#'
#' @examples
#' if (torch::torch_is_installed()) {
#'   x <- matrix(rnorm(100 * 5), nrow = 100)
#'   d <- torch_dist(x, p = 2)
#'   h <- hclust(d)
#'   plot(h)
#' }
#'
#' @export
torch_dist <- function(x) {
  # Ensure x is a numeric matrix
  x <- as.matrix(x)

  N <- nrow(x)

  if (N < 2) {
    stop("Input matrix 'x' must have at least two rows.")
  }

  # Expected length of condensed distance vector
  expected_length <- N * (N - 1) / 2

  # Convert to torch tensor
  x_ten <- torch::torch_tensor(x, dtype = torch::torch_float())

  # Compute pairwise distances (Manhattan)
  pd <- torch::nnf_pdist(x_ten, p = 1)
  pd <- as.numeric(pd)

  # Check if the output length matches expectation
  if (length(pd) != expected_length) {
    stop(sprintf(
      "Mismatch in distance vector length: expected %d, got %d.",
      expected_length, length(pd)
    ))
  }

  # Initialize distance matrix
  mat <- matrix(0, nrow = N, ncol = N)

  # Fill lower triangle
  mat[lower.tri(mat, diag = FALSE)] <- pd

  # Mirror to upper triangle
  mat <- mat + t(mat)

  # Return as dist object
  as.dist(mat)
}
