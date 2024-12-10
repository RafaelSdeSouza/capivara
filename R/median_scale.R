#' Median Scale
#'
#' This function scales a numeric vector by subtracting its median value.
#' It is useful for normalizing data while preserving relative variation.
#'
#' @param x A numeric vector to be scaled.
#'
#' @return A numeric vector of the same length as `x`, where each element is the difference between the original value and the median of `x`.
#'
#' @examples
#' median_scale(c(1, 2, 3, 4, 5))  # Returns c(-2, -1, 0, 1, 2)
#' median_scale(c(10, 15, 20))     # Returns c(-5, 0, 5)
#' median_scale(c(5, NA, 15, 25))  # Handles NA and returns c(-10, NA, 0, 10)
#'
#' @export
median_scale <- function(x) {
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector.")
  }

  x - median(x, na.rm = TRUE)
}
