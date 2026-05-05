.hgc_robust_col_scale <- function(X) {
  X <- as.matrix(X)

  center <- apply(X, 2, stats::median, na.rm = TRUE)
  scale <- apply(X, 2, stats::mad, center = center, constant = 1.4826, na.rm = TRUE)

  fallback <- stats::median(scale[is.finite(scale) & scale > 0], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 0) {
    fallback <- 1
  }

  scale[!is.finite(scale) | scale <= 0] <- fallback

  X <- sweep(X, 2, center, "-")
  X <- sweep(X, 2, scale, "/")
  X[!is.finite(X)] <- 0

  X
}

.hgc_snn_cluster_matrix <- function(features,
                                    Ncomp,
                                    knn_k = 20,
                                    auto_k = TRUE,
                                    max_k = NULL,
                                    verbose = FALSE) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for the HGC-SNN backend.")
  }

  features <- as.matrix(features)
  features[!is.finite(features)] <- 0

  n <- nrow(features)

  if (n < 1L) {
    stop("`features` has no rows.")
  }

  if (!is.finite(Ncomp) || Ncomp < 1L || Ncomp > n) {
    stop("`Ncomp` must be between 1 and the number of feature rows.")
  }

  if (Ncomp == n) {
    return(list(
      labels = seq_len(n),
      knn_k = 0L,
      requested_Ncomp = Ncomp,
      actual_Ncomp = n,
      disconnected = FALSE
    ))
  }

  knn_k <- as.integer(knn_k)
  if (!is.finite(knn_k) || knn_k < 2L) {
    stop("`knn_k` must be an integer >= 2.")
  }

  knn_k <- min(knn_k, n - 1L)

  if (is.null(max_k)) {
    max_k <- min(n - 1L, max(knn_k, 100L))
  } else {
    max_k <- min(as.integer(max_k), n - 1L)
  }

  if (max_k < knn_k) {
    max_k <- knn_k
  }

  labels <- NULL
  actual_k <- NA_integer_
  used_k <- knn_k

  repeat {
    if (verbose) {
      message(sprintf("Building HGC-SNN kNN graph with k = %d...", used_k))
    }

    nn <- RANN::nn2(features, k = used_k + 1L)$nn.idx
    nn <- nn[, -1L, drop = FALSE]

    labels <- capivara_hgc_snn_cut_cpp(nn, as.integer(Ncomp))
    actual_k <- length(unique(labels))

    if (verbose) {
      message(sprintf(
        "HGC-SNN returned %d clusters for requested Ncomp = %d.",
        actual_k,
        Ncomp
      ))
    }

    if (actual_k == Ncomp || !isTRUE(auto_k) || used_k >= max_k) {
      break
    }

    used_k <- min(max_k, max(used_k + 1L, ceiling(1.5 * used_k)))
  }

  disconnected <- actual_k > Ncomp

  if (disconnected) {
    warning(
      "HGC-SNN graph remained disconnected: requested ",
      Ncomp,
      " clusters but obtained ",
      actual_k,
      ". Try increasing `knn_k`, increasing `max_k`, or adding a small `spatial_weight`."
    )
  }

  list(
    labels = as.integer(labels),
    knn_k = used_k,
    requested_Ncomp = as.integer(Ncomp),
    actual_Ncomp = as.integer(actual_k),
    disconnected = disconnected
  )
}

.hgc_snn_wavelength_index <- function(cubedat, flux_mat, wavelength_range = NULL) {
  wavelengths <- if (!is.null(cubedat$axDat)) {
    FITSio::axVec(3, cubedat$axDat)
  } else {
    seq_len(ncol(flux_mat))
  }

  if (is.null(wavelength_range)) {
    return(seq_len(ncol(flux_mat)))
  }

  if (length(wavelength_range) != 2) {
    stop("`wavelength_range` must have length 2.")
  }

  wave_idx <- which(wavelengths >= min(wavelength_range) & wavelengths <= max(wavelength_range))
  if (!length(wave_idx)) {
    stop("No wavelengths fall inside `wavelength_range`.")
  }

  wave_idx
}

#' Fast Capivara segmentation using SNN graph clustering
#'
#' This backend replaces the all-pairs distance matrix plus Ward clustering
#' with an approximate shared-nearest-neighbor graph cut. It is intended for
#' large cubes where \code{\link{segment}} becomes memory-limited.
#'
#' @param input A FITS-like object with `imDat`, or a raw 3D array.
#' @param Ncomp Integer, the number of clusters to form. Defaults to `15`.
#' @param redshift Kept for API compatibility with `segment()`.
#' @param scale_fn Row-wise spectral scaling function. Defaults to `median_scale`.
#' @param target_snr Optional minimum accepted SNR per cluster. When supplied,
#'   Capivara chooses the largest number of clusters whose minimum cluster SNR
#'   remains above this threshold.
#' @param var_cube Optional variance cube matching the input cube. Used only
#'   when \code{target_snr} is supplied.
#' @param k_values Optional candidate cluster counts tested when
#'   \code{target_snr} is supplied. For this scalable backend, the default grid
#'   is capped at 50 clusters.
#' @param wavelength_range Optional wavelength interval used to compute SNR when
#'   \code{target_snr} is supplied.
#' @param snr_stat Either integrated SNR or median per-wavelength SNR when
#'   \code{target_snr} is supplied.
#' @param variance_inflation Multiplicative factor applied to propagated
#'   variances when \code{target_snr} is supplied.
#' @param use_starlet_mask Logical; if \code{TRUE}, build a Sagui-style
#'   white-light starlet mask before clustering.
#' @param collapse_fn Function used to collapse the cube to white light when
#'   \code{use_starlet_mask = TRUE}.
#' @param starlet_J Number of starlet scales when \code{use_starlet_mask = TRUE}.
#' @param starlet_scales Integer vector of scales kept in the reconstruction
#'   when \code{use_starlet_mask = TRUE}.
#' @param include_coarse Logical; include the coarse starlet plane when
#'   \code{use_starlet_mask = TRUE}.
#' @param denoise_k Optional denoising threshold in MAD units when
#'   \code{use_starlet_mask = TRUE}.
#' @param starlet_mode Thresholding mode for the starlet reconstruction when
#'   \code{use_starlet_mask = TRUE}.
#' @param positive_only Logical; keep only positive reconstructed values when
#'   \code{use_starlet_mask = TRUE}.
#' @param mask_mode Either \code{"na"} or \code{"zero"} for masked spaxels
#'   when \code{use_starlet_mask = TRUE}.
#' @param knn_k Number of nearest neighbours for the SNN graph.
#' @param auto_k If TRUE, increase `knn_k` when the graph is too disconnected.
#' @param max_k Maximum k allowed when `auto_k = TRUE`.
#' @param feature_scale Optional column-wise feature scaling after row scaling.
#' @param spatial_weight Optional weight for appending normalized x/y coordinates.
#' @param mask Optional logical spatial mask with dimensions n_row x n_col.
#' @param valid_mode Valid-pixel rule: `"signal"` or `"finite"`.
#' @param return_details Return features, labels, and diagnostics.
#' @param verbose Print progress messages.
#'
#' @return A Capivara-style segmentation object with the same core fields as
#'   \code{\link{segment}}: \code{cluster_map}, \code{header}, \code{axDat},
#'   \code{cluster_snr}, and \code{original_cube}. Scalable-backend diagnostics
#'   are stored under \code{backend_info}.
#' @export
segment_snn <- function(input,
                        Ncomp = 15,
                        redshift = 0,
                        scale_fn = median_scale,
                        target_snr = NULL,
                        var_cube = NULL,
                        k_values = NULL,
                        wavelength_range = NULL,
                        snr_stat = c("integrated", "median_per_wavelength"),
                        variance_inflation = 1,
                        use_starlet_mask = FALSE,
                        collapse_fn = collapse_white_light,
                        starlet_J = 5,
                        starlet_scales = 2:5,
                        include_coarse = FALSE,
                        denoise_k = 0,
                        starlet_mode = c("soft", "hard"),
                        positive_only = TRUE,
                        mask_mode = c("na", "zero"),
                        knn_k = 20,
                        auto_k = TRUE,
                        max_k = NULL,
                        feature_scale = c("none", "robust_col"),
                        spatial_weight = 0,
                        mask = NULL,
                        valid_mode = c("signal", "finite"),
                        return_details = FALSE,
                        verbose = TRUE) {
  feature_scale <- match.arg(feature_scale)
  valid_mode <- match.arg(valid_mode)
  snr_stat <- match.arg(snr_stat)
  starlet_mode <- match.arg(starlet_mode)
  mask_mode <- match.arg(mask_mode)

  if (!is.null(target_snr) && !missing(Ncomp)) {
    stop("Specify either `Ncomp` or `target_snr`, not both.")
  }

  starlet_prep <- .apply_starlet_support(
    input = input,
    use_starlet_mask = use_starlet_mask,
    collapse_fn = collapse_fn,
    starlet_J = starlet_J,
    starlet_scales = starlet_scales,
    include_coarse = include_coarse,
    denoise_k = denoise_k,
    mode = starlet_mode,
    positive_only = positive_only,
    mask_mode = mask_mode
  )

  cubedat <- .as_cubedat(starlet_prep$input)
  cube <- cubedat$imDat

  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("`input$imDat` must be a 3D array with dimensions (n_row, n_col, n_wave).")
  }

  if (!is.function(scale_fn)) {
    stop("`scale_fn` must be a function.")
  }

  if (!is.finite(spatial_weight) || spatial_weight < 0) {
    stop("`spatial_weight` must be >= 0.")
  }

  n_row <- dim(cube)[1]
  n_col <- dim(cube)[2]
  n_wave <- dim(cube)[3]

  IFU2D <- cube_to_matrix(cubedat)
  sn <- .compute_signal_noise(IFU2D)

  if (valid_mode == "signal") {
    valid <- sn$signal > 0
  } else {
    valid <- rowSums(is.finite(IFU2D)) == n_wave
  }

  if (!is.null(mask)) {
    if (!identical(dim(mask), c(n_row, n_col))) {
      stop("`mask` must have dimensions n_row x n_col.")
    }

    valid <- valid & as.vector(mask)
  }

  valid_indices <- which(valid)

  if (!length(valid_indices)) {
    stop("No valid pixels after filtering.")
  }

  if (!is.null(Ncomp) && length(valid_indices) < Ncomp) {
    stop("`Ncomp` is larger than the number of valid pixels.")
  }

  if (verbose) {
    message(sprintf("Valid pixels for HGC-SNN: %d", length(valid_indices)))
  }

  IFU2D_valid <- IFU2D[valid_indices, , drop = FALSE]
  signal_valid <- sn$signal[valid_indices]
  noise_valid <- sn$noise[valid_indices]

  features <- .scale_rows(
    IFU2D_valid,
    scale_fn = scale_fn,
    na_to_zero = TRUE
  )

  if (feature_scale == "robust_col") {
    features <- .hgc_robust_col_scale(features)
  }

  if (spatial_weight > 0) {
    ij <- arrayInd(valid_indices, .dim = c(n_row, n_col))

    x_sd <- stats::sd(ij[, 2])
    y_sd <- stats::sd(ij[, 1])

    if (!is.finite(x_sd) || x_sd <= 0) x_sd <- 1
    if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1

    xy <- cbind(
      x = (ij[, 2] - mean(ij[, 2])) / x_sd,
      y = (ij[, 1] - mean(ij[, 1])) / y_sd
    )

    features <- cbind(features, spatial_weight * xy)
  }

  snr_grid <- NULL

  if (is.null(target_snr)) {
    fit <- .hgc_snn_cluster_matrix(
      features = features,
      Ncomp = Ncomp,
      knn_k = knn_k,
      auto_k = auto_k,
      max_k = max_k,
      verbose = verbose
    )
    clusters <- fit$labels
    cluster_snr <- .compute_cluster_snr(
      clusters = clusters,
      signal_valid = signal_valid,
      noise_valid = noise_valid
    )
    chosen_k <- fit$actual_Ncomp
  } else {
    flux_mat_valid <- IFU2D_valid
    if (is.null(var_cube)) {
      var_mat_valid <- pmax(flux_mat_valid, 0)
    } else {
      var_input <- .as_cubedat(var_cube)
      if (!identical(dim(var_input$imDat), dim(cubedat$imDat))) {
        stop("`var_cube` must have the same dimensions as the input cube.")
      }
      var_mat_valid <- cube_to_matrix(var_input)[valid_indices, , drop = FALSE]
    }
    var_mat_valid[!is.finite(var_mat_valid) | var_mat_valid < 0] <- NA_real_

    wave_idx <- .hgc_snn_wavelength_index(
      cubedat = cubedat,
      flux_mat = flux_mat_valid,
      wavelength_range = wavelength_range
    )

    n_valid <- length(valid_indices)
    if (is.null(k_values)) {
      k_values <- seq(from = min(n_valid, 50L), to = 1L, by = -1L)
    } else {
      k_values <- sort(unique(as.integer(k_values)), decreasing = TRUE)
      k_values <- k_values[k_values >= 1L & k_values <= n_valid]
    }
    if (!length(k_values)) {
      stop("No valid `k_values` to test.")
    }

    fits <- vector("list", length(k_values))
    snr_rows <- vector("list", length(k_values))
    for (i in seq_along(k_values)) {
      k <- k_values[i]
      fit_k <- .hgc_snn_cluster_matrix(
        features = features,
        Ncomp = k,
        knn_k = knn_k,
        auto_k = auto_k,
        max_k = max_k,
        verbose = FALSE
      )
      fits[[i]] <- fit_k
      cluster_snr_k <- .compute_cluster_snr_from_variance(
        flux_mat = flux_mat_valid,
        var_mat = var_mat_valid,
        clusters = fit_k$labels,
        wave_idx = wave_idx,
        snr_stat = snr_stat,
        variance_inflation = variance_inflation
      )
      screen <- .evaluate_cluster_snr_screen(cluster_snr_k, target_snr)

      snr_rows[[i]] <- data.frame(
        Ncomp = k,
        actual_Ncomp = fit_k$actual_Ncomp,
        min_cluster_snr = screen$min_cluster_snr,
        all_clusters_pass = screen$all_clusters_pass
      )
    }
    snr_grid <- do.call(rbind, snr_rows)

    ok <- which(snr_grid$all_clusters_pass)
    if (!length(ok)) {
      stop("No clustering configuration satisfies the target SNR.")
    }

    fit <- fits[[ok[1]]]
    clusters <- fit$labels
    chosen_k <- fit$actual_Ncomp
    cluster_snr <- .compute_cluster_snr_from_variance(
      flux_mat = flux_mat_valid,
      var_mat = var_mat_valid,
      clusters = clusters,
      wave_idx = wave_idx,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation
    )
  }

  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  cluster_map[valid_indices] <- clusters

  out <- list(
    cluster_map = cluster_map,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    cluster_snr = cluster_snr,
    Ncomp = chosen_k,
    snr_grid = snr_grid,
    requested_Ncomp = fit$requested_Ncomp,
    original_cube = cubedat,
    backend = "snn",
    backend_info = list(
      algorithm = "hgc_snn",
      knn_k = fit$knn_k,
      requested_Ncomp = fit$requested_Ncomp,
      actual_Ncomp = fit$actual_Ncomp,
      disconnected = fit$disconnected,
      feature_scale = feature_scale,
      spatial_weight = spatial_weight,
      valid_mode = valid_mode,
      valid_pixels = length(valid_indices)
    )
  )

  if (!is.null(starlet_prep$starlet_info)) {
    out$starlet_info <- starlet_prep$starlet_info
  }

  if (return_details) {
    out <- c(out, list(
      valid_indices = valid_indices,
      signal = sn$signal,
      noise = sn$noise,
      features = features,
      labels = clusters
    ))
  }

  out
}

#' @rdname segment_snn
#' @export
segment_hgc_snn <- segment_snn
