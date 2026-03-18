#' Backward-compatible Alias for Missing-data-safe Segmentation
#'
#' \code{segment_robust()} now forwards to \code{\link{segment_masked}}.
#' This keeps the exported API stable while making the function actually safe
#' for masked cubes and spectra with missing channels.
#'
#' @inheritParams segment_masked
#' @return A list with the same fields returned by \code{\link{segment}}.
#' @seealso \code{\link{segment_masked}}, \code{\link{segment_starlet}}
#' @export
segment_robust <- function(input, Ncomp = 5, redshift = 0, scale_fn = median_scale) {
  segment_masked(
    input = input,
    Ncomp = Ncomp,
    redshift = redshift,
    scale_fn = scale_fn
  )
}
