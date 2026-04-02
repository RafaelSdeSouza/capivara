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
torch_dist <- function(x, p = 1, device = NULL) {
  # Ensure x is a numeric matrix
  x <- as.matrix(x)

  N <- nrow(x)

  if (N < 2) {
    stop("Input matrix 'x' must have at least two rows.")
  }

  # Expected length of condensed distance vector
  expected_length <- N * (N - 1) / 2

  if (!requireNamespace("torch", quietly = TRUE) || !isTRUE(torch::torch_is_installed())) {
    method <- if (p == 1) "manhattan" else if (p == 2) "euclidean" else NULL
    if (is.null(method)) {
      stop("Fallback distance only supports p = 1 or p = 2 when torch is unavailable.")
    }
    return(stats::dist(x, method = method))
  }

  # Convert to torch tensor
  if (is.null(device)) {
    x_ten <- torch::torch_tensor(x, dtype = torch::torch_float())
  } else {
    x_ten <- torch::torch_tensor(x, dtype = torch::torch_float(), device = device)
  }

  # Compute pairwise distances
  pd <- torch::nnf_pdist(x_ten, p = p)
  pd <- as.numeric(pd)

  # Check if the output length matches expectation
  if (length(pd) != expected_length) {
    stop(sprintf(
      "Mismatch in distance vector length: expected %d, got %d.",
      expected_length, length(pd)
    ))
  }

  structure(
    pd,
    Size = N,
    Diag = FALSE,
    Upper = FALSE,
    method = if (p == 1) "manhattan" else if (p == 2) "euclidean" else paste0("minkowski_p_", p),
    class = "dist"
  )
}
