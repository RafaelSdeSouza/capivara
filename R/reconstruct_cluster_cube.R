#' Reconstruct a Model Cube from Cluster Representative Spectra
#'
#' This function builds a model cube from a segmentation result by assigning one
#' representative spectrum to every pixel in each cluster. It is useful for
#' visualization, denoising, and downstream spectral analysis where a
#' cluster-level spectral template is desired.
#'
#' @param cluster_result A list returned by a segmentation function.
#' @param template Representative spectrum to assign to each cluster.
#'   \code{"median"} is robust, \code{"mean"} is the arithmetic cluster mean,
#'   and \code{"weighted_mean"} uses inverse-variance weighting.
#' @param var_cube Optional variance cube. Required when
#'   \code{template = "weighted_mean"}.
#' @param variance_inflation Multiplicative factor applied to propagated
#'   variances when \code{var_cube} is supplied.
#' @param preserve_mask Logical; if \code{TRUE}, channels that were non-finite in
#'   the original cube remain masked in the reconstructed cube.
#' @param fill_mode Either \code{"na"} or \code{"zero"} for unassigned pixels and,
#'   when \code{preserve_mask = TRUE}, for masked spectral channels.
#' @param return_residual Logical; if \code{TRUE}, also return the residual cube.
#'
#' @return A list with the reconstructed cube, optional residual cube, the
#'   template spectra, a cluster summary object, and a global flux comparison
#'   between the original and reconstructed cubes.
#' @export
reconstruct_cluster_cube <- function(cluster_result,
                                     template = c("median", "mean", "weighted_mean"),
                                     var_cube = NULL,
                                     variance_inflation = 1,
                                     preserve_mask = TRUE,
                                     fill_mode = c("na", "zero"),
                                     return_residual = TRUE) {
  template <- match.arg(template)
  fill_mode <- match.arg(fill_mode)

  if (identical(template, "weighted_mean") && is.null(var_cube)) {
    stop("`var_cube` must be supplied when `template = \"weighted_mean\"`.")
  }

  cubedat <- .as_cubedat(cluster_result$original_cube)
  flux_mat <- cube_to_matrix(cubedat)
  cluster_vec <- as.vector(cluster_result$cluster_map)

  if (length(cluster_vec) != nrow(flux_mat)) {
    stop("Cluster map dimensions do not match the input cube dimensions.")
  }

  summary <- summarize_cluster_spectra(
    cluster_result = cluster_result,
    var_cube = var_cube,
    variance_inflation = variance_inflation
  )

  template_spectra <- switch(
    template,
    median = summary$median_spectra,
    mean = summary$mean_spectra,
    weighted_mean = summary$weighted_mean_spectra
  )

  model_mat <- matrix(NA_real_, nrow = nrow(flux_mat), ncol = ncol(flux_mat))

  for (i in seq_along(summary$cluster_ids)) {
    idx <- which(cluster_vec == summary$cluster_ids[i])
    if (!length(idx)) {
      next
    }

    model_mat[idx, ] <- matrix(
      template_spectra[i, ],
      nrow = length(idx),
      ncol = ncol(flux_mat),
      byrow = TRUE
    )
  }

  if (isTRUE(preserve_mask)) {
    orig_finite <- is.finite(flux_mat)
    if (fill_mode == "na") {
      model_mat[!orig_finite & !is.na(cluster_vec)] <- NA_real_
    } else {
      model_mat[!orig_finite & !is.na(cluster_vec)] <- 0
    }
  }

  unassigned <- is.na(cluster_vec)
  if (fill_mode == "zero") {
    model_mat[unassigned, ] <- 0
  }

  original_sum_spectrum <- colSums(flux_mat[!unassigned, , drop = FALSE], na.rm = TRUE)
  model_sum_spectrum <- colSums(model_mat[!unassigned, , drop = FALSE], na.rm = TRUE)
  flux_difference <- model_sum_spectrum - original_sum_spectrum
  flux_check <- list(
    original_sum_spectrum = original_sum_spectrum,
    model_sum_spectrum = model_sum_spectrum,
    difference = flux_difference,
    max_abs_difference = max(abs(flux_difference), na.rm = TRUE)
  )

  model_cube <- cubedat
  model_cube$imDat <- array(model_mat, dim = dim(cubedat$imDat))

  template_df <- data.frame(
    cluster = summary$cluster_ids,
    template_spectra,
    check.names = FALSE
  )

  out <- list(
    model_cube = model_cube,
    template = template,
    template_spectra = template_df,
    cluster_summary = summary,
    flux_check = flux_check
  )

  if (isTRUE(return_residual)) {
    residual_cube <- cubedat
    residual_cube$imDat <- array(flux_mat - model_mat, dim = dim(cubedat$imDat))
    out$residual_cube <- residual_cube
  }

  out
}

#' Reconstruct a Flux-preserving Model Cube from Cluster Spectra
#'
#' This is a convenience wrapper around \code{\link{reconstruct_cluster_cube}}
#' that uses the arithmetic cluster mean while preserving the original spectral
#' mask. In this mode, the reconstructed cube preserves the summed flux spectrum
#' of the segmented cube on the observed support, which makes it a sensible
#' default for later spectral fitting.
#'
#' @param cluster_result A list returned by a segmentation function.
#' @param fill_mode Either \code{"na"} or \code{"zero"} for unassigned pixels
#'   and masked channels.
#' @param return_residual Logical; if \code{TRUE}, also return the residual cube.
#'
#' @return The same structure returned by \code{\link{reconstruct_cluster_cube}}.
#' @export
reconstruct_flux_preserving_cube <- function(cluster_result,
                                             fill_mode = c("na", "zero"),
                                             return_residual = TRUE) {
  reconstruct_cluster_cube(
    cluster_result = cluster_result,
    template = "mean",
    preserve_mask = TRUE,
    fill_mode = fill_mode,
    return_residual = return_residual
  )
}
