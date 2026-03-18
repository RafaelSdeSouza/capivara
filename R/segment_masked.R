#' Cluster a Data Cube While Preserving `segment()` Semantics for Missing Data
#'
#' This function is the missing-data-safe counterpart of \code{\link{segment}}.
#' It uses the same signal filtering, row-wise scaling, Manhattan distances, and
#' Ward.D2 clustering as \code{segment()}, but replaces non-finite scaled values
#' with zero before computing distances. This makes it suitable for masked cubes
#' or spectra containing isolated \code{NA} / \code{Inf} channels.
#'
#' @param input A FITS-like object with an \code{imDat} cube, or a raw
#'   3-D numeric array.
#' @param Ncomp Integer, the number of clusters to form.
#' @param redshift Numeric redshift placeholder kept for API compatibility.
#' @param scale_fn Row-wise scaling function. Defaults to \code{\link{median_scale}}.
#'
#' @return A list with the same fields returned by \code{\link{segment}}.
#'
#' @seealso \code{\link{segment}}, \code{\link{segment_starlet}}
#'
#' @examples
#' cube <- array(runif(6 * 6 * 20), dim = c(6, 6, 20))
#' cube[1, 1, 1:3] <- NA
#' res <- segment_masked(list(imDat = cube), Ncomp = 3)
#'
#' @export
segment_masked <- function(input, Ncomp = 5, redshift = 0, scale_fn = median_scale) {
  .segment_core(
    input = input,
    Ncomp = Ncomp,
    redshift = redshift,
    scale_fn = scale_fn,
    na_to_zero = TRUE
  )
}
