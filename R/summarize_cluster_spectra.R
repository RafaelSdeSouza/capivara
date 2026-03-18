#' Summarize Cluster Spectra with Optional Variance Propagation
#'
#' This is the recommended post-segmentation summary layer for Capivara. It keeps
#' the segmentation step separate from downstream spectral products and can return
#' median spectra, arithmetic mean spectra, summed spectra, and
#' inverse-variance-weighted mean spectra.
#'
#' @param cluster_result A list returned by a segmentation function.
#' @param var_cube Optional variance cube matching the dimensions of the original
#'   flux cube. It can be a FITS-like list with \code{imDat} or a raw 3-D array.
#' @param variance_inflation Multiplicative factor applied to propagated
#'   variances. Use this to account for covariance if needed.
#'
#' @return A list with wavelength coordinates, cluster ids, cluster sizes,
#'   per-wavelength finite counts, median spectra, mean spectra, summed spectra,
#'   and, when \code{var_cube} is supplied, propagated sum variances and
#'   inverse-variance-weighted means.
#' @export
summarize_cluster_spectra <- function(cluster_result,
                                      var_cube = NULL,
                                      variance_inflation = 1) {
  cubedat <- .as_cubedat(cluster_result$original_cube)
  flux_mat <- cube_to_matrix(cubedat)
  cluster_vec <- as.vector(cluster_result$cluster_map)

  if (length(cluster_vec) != nrow(flux_mat)) {
    stop("Cluster map dimensions do not match the input cube dimensions.")
  }

  cluster_ids <- sort(unique(cluster_vec[!is.na(cluster_vec)]))
  wavelengths <- if (!is.null(cubedat$axDat)) {
    FITSio::axVec(3, cubedat$axDat)
  } else {
    seq_len(ncol(flux_mat))
  }
  wave_names <- as.character(wavelengths)

  median_spectra <- matrix(
    NA_real_,
    nrow = length(cluster_ids),
    ncol = ncol(flux_mat),
    dimnames = list(as.character(cluster_ids), wave_names)
  )
  mean_spectra <- median_spectra
  sum_spectra <- median_spectra
  finite_counts <- matrix(
    0L,
    nrow = length(cluster_ids),
    ncol = ncol(flux_mat),
    dimnames = list(as.character(cluster_ids), wave_names)
  )
  n_spaxels <- integer(length(cluster_ids))

  has_var <- !is.null(var_cube)
  if (has_var) {
    var_input <- .as_cubedat(var_cube)
    if (!identical(dim(var_input$imDat), dim(cubedat$imDat))) {
      stop("`var_cube` must have the same dimensions as `cluster_result$original_cube$imDat`.")
    }
    var_mat <- cube_to_matrix(var_input)

    sum_variance <- median_spectra
    weighted_mean_spectra <- median_spectra
    weighted_mean_variance <- median_spectra
  }

  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    idx <- which(cluster_vec == cluster_id)
    X <- flux_mat[idx, , drop = FALSE]

    n_spaxels[i] <- length(idx)
    median_spectra[i, ] <- apply(X, 2, stats::median, na.rm = TRUE)
    sum_spectra[i, ] <- colSums(X, na.rm = TRUE)
    finite_counts[i, ] <- colSums(is.finite(X))
    mean_spectra[i, ] <- sum_spectra[i, ] / finite_counts[i, ]
    mean_spectra[i, finite_counts[i, ] == 0] <- NA_real_

    if (has_var) {
      V <- var_mat[idx, , drop = FALSE]
      V[!is.finite(V) | V <= 0] <- NA_real_

      sum_variance[i, ] <- variance_inflation * colSums(V, na.rm = TRUE)

      w <- 1 / V
      w[!is.finite(w)] <- 0
      weight_sum <- colSums(w)
      weighted_mean <- colSums(w * X, na.rm = TRUE) / weight_sum
      weighted_var <- variance_inflation / weight_sum

      weighted_mean[!is.finite(weighted_mean)] <- NA_real_
      weighted_var[!is.finite(weighted_var)] <- NA_real_

      weighted_mean_spectra[i, ] <- weighted_mean
      weighted_mean_variance[i, ] <- weighted_var
    }
  }

  out <- list(
    wavelength = wavelengths,
    cluster_ids = cluster_ids,
    n_spaxels = stats::setNames(n_spaxels, cluster_ids),
    finite_counts = finite_counts,
    median_spectra = median_spectra,
    mean_spectra = mean_spectra,
    sum_spectra = sum_spectra
  )

  if (has_var) {
    out <- c(out, list(
      sum_variance = sum_variance,
      weighted_mean_spectra = weighted_mean_spectra,
      weighted_mean_variance = weighted_mean_variance
    ))
  }

  out
}
