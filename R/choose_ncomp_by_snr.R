.compute_cluster_snr_from_variance <- function(flux_mat,
                                               var_mat,
                                               clusters,
                                               wave_idx,
                                               snr_stat = c("integrated", "median_per_wavelength"),
                                               variance_inflation = 1) {
  snr_stat <- match.arg(snr_stat)

  vapply(sort(unique(clusters)), function(cluster_id) {
    idx <- which(clusters == cluster_id)

    bin_flux <- colSums(flux_mat[idx, wave_idx, drop = FALSE], na.rm = TRUE)
    bin_var <- variance_inflation * colSums(var_mat[idx, wave_idx, drop = FALSE], na.rm = TRUE)
    bin_snr <- bin_flux / sqrt(bin_var)
    bin_snr[!is.finite(bin_snr)] <- NA_real_

    if (snr_stat == "integrated") {
      total_flux <- sum(bin_flux, na.rm = TRUE)
      total_var <- sum(bin_var, na.rm = TRUE)
      if (!is.finite(total_var) || total_var <= 0) {
        return(NA_real_)
      }
      total_flux / sqrt(total_var)
    } else {
      stats::median(bin_snr, na.rm = TRUE)
    }
  }, numeric(1))
}

#' Choose the Number of Clusters from a Target SNR Cut
#'
#' This function computes the Capivara dendrogram once, evaluates a grid of cuts,
#' and chooses the largest number of clusters whose minimum bin SNR remains above
#' the requested threshold.
#'
#' @param input A FITS-like object with an \code{imDat} cube, or a raw 3-D array.
#' @param target_snr Minimum accepted SNR per cluster.
#' @param var_cube Optional variance cube matching the flux cube dimensions. When
#'   omitted, a simple Poisson-like approximation \code{var = pmax(flux, 0)} is used.
#' @param k_values Optional integer vector of candidate cluster counts. By default,
#'   all valid counts are tested from the maximum down to 1.
#' @param wavelength_range Optional numeric vector of length 2 selecting the
#'   wavelength interval used to compute SNR.
#' @param redshift Numeric redshift placeholder kept for API compatibility.
#' @param scale_fn Row-wise scaling function used during segmentation.
#' @param snr_stat Either integrated SNR across the chosen window or the median
#'   per-wavelength SNR inside that window.
#' @param variance_inflation Multiplicative factor applied to propagated variances.
#' @param na_safe Logical; when \code{TRUE}, use the missing-data-safe scaling path.
#'
#' @return A segmentation result with the chosen \code{Ncomp}, the per-cluster
#'   SNRs at the selected cut, and a table showing the SNR screen across tested
#'   values of \code{k}.
#' @export
choose_ncomp_by_snr <- function(input,
                                target_snr,
                                var_cube = NULL,
                                k_values = NULL,
                                wavelength_range = NULL,
                                redshift = 0,
                                scale_fn = median_scale,
                                snr_stat = c("integrated", "median_per_wavelength"),
                                variance_inflation = 1,
                                na_safe = TRUE) {
  snr_stat <- match.arg(snr_stat)
  cubedat <- .as_cubedat(input)

  details <- .segment_core(
    input = cubedat,
    Ncomp = NULL,
    redshift = redshift,
    scale_fn = scale_fn,
    na_to_zero = na_safe,
    return_details = TRUE
  )

  flux_mat <- cube_to_matrix(cubedat)[details$valid_indices, , drop = FALSE]
  if (is.null(var_cube)) {
    var_mat <- pmax(flux_mat, 0)
  } else {
    var_input <- .as_cubedat(var_cube)
    if (!identical(dim(var_input$imDat), dim(cubedat$imDat))) {
      stop("`var_cube` must have the same dimensions as the input cube.")
    }
    var_mat <- cube_to_matrix(var_input)[details$valid_indices, , drop = FALSE]
  }
  var_mat[!is.finite(var_mat) | var_mat < 0] <- NA_real_

  wavelengths <- if (!is.null(cubedat$axDat)) {
    FITSio::axVec(3, cubedat$axDat)
  } else {
    seq_len(ncol(flux_mat))
  }

  if (is.null(wavelength_range)) {
    wave_idx <- seq_len(ncol(flux_mat))
  } else {
    if (length(wavelength_range) != 2) {
      stop("`wavelength_range` must have length 2.")
    }
    wave_idx <- which(wavelengths >= min(wavelength_range) & wavelengths <= max(wavelength_range))
    if (!length(wave_idx)) {
      stop("No wavelengths fall inside `wavelength_range`.")
    }
  }

  n_valid <- length(details$valid_indices)
  if (is.null(k_values)) {
    k_values <- seq(from = n_valid, to = 1L, by = -1L)
  } else {
    k_values <- sort(unique(as.integer(k_values)), decreasing = TRUE)
    k_values <- k_values[k_values >= 1L & k_values <= n_valid]
  }
  if (!length(k_values)) {
    stop("No valid `k_values` to test.")
  }

  snr_grid <- lapply(k_values, function(k) {
    clusters <- stats::cutree(details$hclust, k = k)
    cluster_snr <- .compute_cluster_snr_from_variance(
      flux_mat = flux_mat,
      var_mat = var_mat,
      clusters = clusters,
      wave_idx = wave_idx,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation
    )

    data.frame(
      Ncomp = k,
      min_cluster_snr = min(cluster_snr, na.rm = TRUE),
      all_clusters_pass = all(cluster_snr >= target_snr, na.rm = TRUE)
    )
  })
  snr_grid <- do.call(rbind, snr_grid)

  ok <- which(snr_grid$all_clusters_pass)
  if (!length(ok)) {
    stop("No clustering configuration satisfies the target SNR.")
  }

  best_k <- snr_grid$Ncomp[ok[1]]
  best_clusters <- stats::cutree(details$hclust, k = best_k)
  best_cluster_snr <- .compute_cluster_snr_from_variance(
    flux_mat = flux_mat,
    var_mat = var_mat,
    clusters = best_clusters,
    wave_idx = wave_idx,
    snr_stat = snr_stat,
    variance_inflation = variance_inflation
  )

  cluster_map <- matrix(NA_integer_, nrow = dim(cubedat$imDat)[1], ncol = dim(cubedat$imDat)[2])
  cluster_map[details$valid_indices] <- best_clusters

  list(
    cluster_map = cluster_map,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    Ncomp = best_k,
    cluster_snr = best_cluster_snr,
    snr_grid = snr_grid,
    hclust = details$hclust,
    original_cube = cubedat
  )
}
