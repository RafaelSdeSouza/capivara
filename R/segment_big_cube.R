#' Scalable Segmentation for Very Large IFU Cubes
#'
#' This is the user-facing large-cube entry point for Capivara. It forwards to
#' \code{\link{segment_blockward}}, which approximates Ward-style clustering by
#' aggregating spectra in spatial blocks before clustering and then assigning
#' labels back at full spatial resolution.
#'
#' Use this function when the exact \code{\link{segment}} /
#' \code{\link{segment_masked}} workflow would require too much RAM because it
#' must build an all-pairs distance object over all valid spaxels.
#'
#' @inheritParams segment_blockward
#' @return A list with the same fields returned by \code{\link{segment_blockward}}.
#' @seealso \code{\link{segment}}, \code{\link{segment_masked}},
#'   \code{\link{segment_blockward}}, \code{\link{segment_starlet}}
#' @export
segment_big_cube <- function(input,
                             Ncomp = 5,
                             redshift = 0,
                             scale_fn = NULL,
                             block_size = NULL,
                             ram_gb = 32,
                             frac_for_dist = 0.25,
                             m_cap = NULL,
                             agg_fn = c("median", "mean"),
                             valid_mode = c("signal", "finite_frac"),
                             min_finite_frac = 0.80,
                             min_block_finite_frac = 0.80,
                             use_pca = FALSE,
                             pca_k = 30,
                             seed = 42,
                             scale_first = TRUE,
                             full_res_assign = TRUE,
                             dist_method = c("capivara_l1", "euclidean"),
                             pixel_assign = c("l1_nearest_center", "l2_nearest_center"),
                             polish_iters = 0L,
                             refine = TRUE,
                             refine_radius = 1L,
                             refine_max_pixels = 50000L,
                             verbose = TRUE) {
  segment_blockward(
    input = input,
    Ncomp = Ncomp,
    redshift = redshift,
    scale_fn = scale_fn,
    block_size = block_size,
    ram_gb = ram_gb,
    frac_for_dist = frac_for_dist,
    m_cap = m_cap,
    agg_fn = agg_fn,
    valid_mode = valid_mode,
    min_finite_frac = min_finite_frac,
    min_block_finite_frac = min_block_finite_frac,
    use_pca = use_pca,
    pca_k = pca_k,
    seed = seed,
    scale_first = scale_first,
    full_res_assign = full_res_assign,
    dist_method = dist_method,
    pixel_assign = pixel_assign,
    polish_iters = polish_iters,
    refine = refine,
    refine_radius = refine_radius,
    refine_max_pixels = refine_max_pixels,
    verbose = verbose
  )
}
