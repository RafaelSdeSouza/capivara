#' Estimate segmentation memory requirements
#'
#' Estimate the RAM pressure from the exact Ward backend before allocating the
#' all-pairs distance vector. This is a lower bound: clustering and R object
#' overhead require additional memory.
#'
#' @param input A FITS-like object with \code{imDat}, a raw 3-D cube array, or
#'   an integer giving the number of valid pixels.
#' @param valid_mode Valid-pixel rule used when \code{input} is a cube:
#'   \code{"signal"} matches \code{\link{segment}}, while \code{"finite"}
#'   requires every spectral channel to be finite.
#' @param knn_k Number of nearest neighbours used for the sparse graph estimate.
#' @param overhead_factor Multiplicative factor used for a conservative exact
#'   Ward RAM estimate beyond the condensed distance vector alone.
#'
#' @return A one-row data frame with valid-pixel count and approximate RAM
#'   estimates in GB.
#' @export
estimate_segment_memory <- function(input,
                                    valid_mode = c("signal", "finite"),
                                    knn_k = 20,
                                    overhead_factor = 2) {
  valid_mode <- match.arg(valid_mode)

  if (length(input) == 1L && is.numeric(input)) {
    n_valid <- as.integer(input)
  } else {
    cubedat <- .as_cubedat(input)
    cube <- cubedat$imDat

    if (!is.array(cube) || length(dim(cube)) != 3L) {
      stop("`input$imDat` must be a 3D array with dimensions (n_row, n_col, n_wave).")
    }

    IFU2D <- cube_to_matrix(cubedat)
    if (valid_mode == "signal") {
      signal <- rowSums(IFU2D, na.rm = TRUE)
      signal[!is.finite(signal)] <- 0
      n_valid <- sum(signal > 0)
    } else {
      n_valid <- sum(rowSums(is.finite(IFU2D)) == ncol(IFU2D))
    }
  }

  if (!is.finite(n_valid) || n_valid < 0L) {
    stop("Could not determine a valid non-negative pixel count.")
  }

  knn_k <- max(0L, min(as.integer(knn_k), max(0L, n_valid - 1L)))
  exact_distance_gb <- n_valid * (n_valid - 1) / 2 * 8 / 1024^3
  sparse_graph_gb <- n_valid * knn_k * (4 + 8) / 1024^3

  data.frame(
    valid_pixels = n_valid,
    exact_distance_gb = exact_distance_gb,
    exact_conservative_gb = exact_distance_gb * overhead_factor,
    sparse_knn_graph_gb = sparse_graph_gb,
    knn_k = knn_k,
    row.names = NULL
  )
}
