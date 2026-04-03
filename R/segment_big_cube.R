.segment_bigcube_safe_scale <- function(v) {
  ok <- is.finite(v)
  vv <- v[ok]
  if (length(vv) < 2L) {
    return(rep(0, length(v)))
  }

  med <- stats::median(vv)
  sc <- stats::mad(vv, center = med, constant = 1)
  if (!is.finite(sc) || sc == 0) {
    sc <- stats::sd(vv)
  }
  if (!is.finite(sc) || sc == 0) {
    return(rep(0, length(v)))
  }

  out <- (v - med) / sc
  out[!is.finite(out)] <- 0
  out
}

.segment_bigcube_choose_block_size <- function(n_row,
                                               n_col,
                                               ram_gb = 32,
                                               frac_for_dist = 0.25,
                                               m_cap = NULL) {
  dist_bytes <- ram_gb * 1e9 * frac_for_dist
  m_budget <- floor(sqrt(dist_bytes / 4))
  if (is.null(m_cap)) {
    m_cap <- if (ram_gb <= 64) 15000L else 50000L
  }
  m_target <- min(m_budget, m_cap)

  for (b in 2:128) {
    m <- ceiling(n_row / b) * ceiling(n_col / b)
    if (m <= m_target) {
      return(list(block_size = b, n_blocks = m, m_target = m_target))
    }
  }

  b <- 128L
  m <- ceiling(n_row / b) * ceiling(n_col / b)
  list(block_size = b, n_blocks = m, m_target = m_target)
}

.segment_bigcube_boundary_pixels <- function(lbl_map, r = 1L) {
  nr <- nrow(lbl_map)
  nc <- ncol(lbl_map)
  out <- matrix(FALSE, nr, nc)
  shifts <- expand.grid(dr = -r:r, dc = -r:r)
  shifts <- shifts[!(shifts$dr == 0 & shifts$dc == 0), , drop = FALSE]

  for (k in seq_len(nrow(shifts))) {
    dr <- shifts$dr[k]
    dc <- shifts$dc[k]
    r1 <- max(1, 1 + dr)
    r2 <- min(nr, nr + dr)
    c1 <- max(1, 1 + dc)
    c2 <- min(nc, nc + dc)

    a <- lbl_map[r1:r2, c1:c2, drop = FALSE]
    b <- lbl_map[(r1 - dr):(r2 - dr), (c1 - dc):(c2 - dc), drop = FALSE]
    diff <- (a != b) & !is.na(a) & !is.na(b)
    out[r1:r2, c1:c2] <- out[r1:r2, c1:c2] | diff
  }

  which(out, arr.ind = TRUE)
}

.segment_bigcube_medoid_index <- function(M, metric = c("l1", "l2")) {
  metric <- match.arg(metric)

  if (nrow(M) <= 1L) {
    return(1L)
  }

  scores <- numeric(nrow(M))
  for (i in seq_len(nrow(M))) {
    ref <- matrix(M[i, ], nrow = nrow(M), ncol = ncol(M), byrow = TRUE)
    if (metric == "l2") {
      scores[i] <- sum(rowSums((M - ref)^2))
    } else {
      scores[i] <- sum(rowSums(abs(M - ref)))
    }
  }

  which.min(scores)
}

.segment_bigcube_assign_nearest <- function(X, centers, metric = c("l1", "l2")) {
  metric <- match.arg(metric)
  n <- nrow(X)
  K <- nrow(centers)
  out <- integer(n)
  chunk <- 2000L

  for (i in seq(1, n, by = chunk)) {
    j <- min(i + chunk - 1L, n)
    XX <- X[i:j, , drop = FALSE]

    if (metric == "l2") {
      c2 <- rowSums(centers^2)
      x2 <- rowSums(XX^2)
      dots <- XX %*% t(centers)
      D2 <- matrix(rep(x2, times = K), ncol = K) +
        matrix(rep(c2, each = nrow(XX)), ncol = K) -
        2 * dots
      out[i:j] <- max.col(-D2)
    } else {
      D <- matrix(0, nrow = nrow(XX), ncol = K)
      for (k in seq_len(K)) {
        D[, k] <- rowSums(abs(XX - matrix(centers[k, ], nrow = nrow(XX), ncol = ncol(XX), byrow = TRUE)))
      }
      out[i:j] <- max.col(-D)
    }
  }

  out
}

.segment_bigcube_compute_cluster_medoids <- function(X,
                                                     labels,
                                                     K,
                                                     metric = c("l1", "l2"),
                                                     fallback_centers = NULL) {
  metric <- match.arg(metric)
  centers <- matrix(NA_real_, nrow = K, ncol = ncol(X))

  for (k in seq_len(K)) {
    idx <- which(labels == k)
    if (!length(idx)) {
      if (!is.null(fallback_centers)) {
        centers[k, ] <- fallback_centers[k, ]
      }
      next
    }
    mk <- .segment_bigcube_medoid_index(X[idx, , drop = FALSE], metric = metric)
    centers[k, ] <- X[idx[mk], ]
  }

  centers
}

.segment_bigcube_assign_map <- function(feat_pixels,
                                        centers,
                                        valid_pix,
                                        n_row,
                                        n_col,
                                        center_metric = c("l1", "l2"),
                                        polish_iters = 0L,
                                        refine = TRUE,
                                        refine_radius = 1L,
                                        refine_max_pixels = 50000L,
                                        seed = 42,
                                        verbose = FALSE) {
  center_metric <- match.arg(center_metric)

  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  pix_lbl <- .segment_bigcube_assign_nearest(feat_pixels, centers, metric = center_metric)

  if (polish_iters > 0L) {
    if (verbose) {
      message(sprintf("Polishing %d iterations in full pixel space...", polish_iters))
    }
    for (it in seq_len(polish_iters)) {
      centers <- .segment_bigcube_compute_cluster_medoids(
        feat_pixels,
        labels = pix_lbl,
        K = nrow(centers),
        metric = center_metric,
        fallback_centers = centers
      )
      pix_lbl <- .segment_bigcube_assign_nearest(feat_pixels, centers, metric = center_metric)
    }
  }

  cluster_map[valid_pix] <- pix_lbl

  if (refine) {
    if (verbose) {
      message("Boundary refinement...")
    }
    ref_centers <- .segment_bigcube_compute_cluster_medoids(
      feat_pixels,
      labels = pix_lbl,
      K = nrow(centers),
      metric = center_metric,
      fallback_centers = centers
    )
    bp <- .segment_bigcube_boundary_pixels(cluster_map, r = refine_radius)
    if (nrow(bp) > refine_max_pixels) {
      set.seed(seed)
      bp <- bp[sample.int(nrow(bp), refine_max_pixels), , drop = FALSE]
    }
    if (nrow(bp) > 0L) {
      lin <- bp[, 1] + (bp[, 2] - 1L) * n_row
      pos <- match(lin, valid_pix)
      ok <- !is.na(pos)
      if (any(ok)) {
        pos <- pos[ok]
        lin <- lin[ok]
        new_lbl <- .segment_bigcube_assign_nearest(
          feat_pixels[pos, , drop = FALSE],
          ref_centers,
          metric = center_metric
        )
        cluster_map[lin] <- new_lbl
        pix_lbl[pos] <- new_lbl
        centers <- ref_centers
      }
    }
  }

  list(
    cluster_map = cluster_map,
    pixel_labels = pix_lbl,
    centers = centers
  )
}

.segment_bigcube_vector_cap <- function(ram_gb = 32) {
  if (ram_gb <= 16) {
    return(4000L)
  }
  if (ram_gb <= 32) {
    return(8000L)
  }
  12000L
}

.segment_bigcube_vector_cluster <- function(feat_pixels,
                                            valid_pix,
                                            n_row,
                                            n_col,
                                            Ncomp = NULL,
                                            target_snr = NULL,
                                            flux_mat_valid,
                                            var_mat_valid = NULL,
                                            wavelengths,
                                            wavelength_range = NULL,
                                            k_values = NULL,
                                            snr_stat = c("integrated", "median_per_wavelength"),
                                            variance_inflation = 1,
                                            verbose = FALSE) {
  snr_stat <- match.arg(snr_stat)

  if (verbose) {
    message(sprintf(
      "Using fastcluster::hclust.vector on %d valid pixels...",
      nrow(feat_pixels)
    ))
  }

  hc <- fastcluster::hclust.vector(
    feat_pixels,
    method = "ward",
    metric = "euclidean"
  )

  if (is.null(target_snr)) {
    chosen_k <- Ncomp
    pix_lbl <- stats::cutree(hc, k = chosen_k)
    snr_grid <- NULL
  } else {
    if (is.null(var_mat_valid)) {
      stop("`var_mat_valid` must be supplied when `target_snr` is used.")
    }

    if (is.null(wavelength_range)) {
      wave_idx <- seq_len(ncol(flux_mat_valid))
    } else {
      if (length(wavelength_range) != 2) {
        stop("`wavelength_range` must have length 2.")
      }
      wave_idx <- which(
        wavelengths >= min(wavelength_range) &
          wavelengths <= max(wavelength_range)
      )
      if (!length(wave_idx)) {
        stop("No wavelengths fall inside `wavelength_range`.")
      }
    }

    max_k <- nrow(feat_pixels)
    if (is.null(k_values)) {
      k_values <- seq(from = min(max_k, 50L), to = 1L, by = -1L)
    } else {
      k_values <- sort(unique(as.integer(k_values)), decreasing = TRUE)
      k_values <- k_values[k_values >= 1L & k_values <= max_k]
    }
    if (!length(k_values)) {
      stop("No valid `k_values` to test.")
    }

    snr_grid <- lapply(k_values, function(k) {
      pix_lbl_k <- stats::cutree(hc, k = k)
      cluster_snr_k <- .compute_cluster_snr_from_variance(
        flux_mat = flux_mat_valid,
        var_mat = var_mat_valid,
        clusters = pix_lbl_k,
        wave_idx = wave_idx,
        snr_stat = snr_stat,
        variance_inflation = variance_inflation
      )
      screen <- .evaluate_cluster_snr_screen(cluster_snr_k, target_snr)

      data.frame(
        Ncomp = k,
        min_cluster_snr = screen$min_cluster_snr,
        all_clusters_pass = screen$all_clusters_pass
      )
    })
    snr_grid <- do.call(rbind, snr_grid)

    ok <- which(snr_grid$all_clusters_pass)
    if (!length(ok)) {
      stop("No clustering configuration satisfies the target SNR.")
    }

    chosen_k <- snr_grid$Ncomp[ok[1]]
    pix_lbl <- stats::cutree(hc, k = chosen_k)
  }

  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  cluster_map[valid_pix] <- pix_lbl

  list(
    cluster_map = cluster_map,
    pixel_labels = pix_lbl,
    hclust = hc,
    Ncomp = chosen_k,
    snr_grid = snr_grid
  )
}

#' Scalable Segmentation for Very Large IFU Cubes via Block Medoids
#'
#' This is the user-facing large-cube entry point for Capivara. It keeps the
#' exact \code{\link{segment}} workflow philosophy, but replaces the full
#' all-pairs pixel clustering with a scalable approximation based on spatial
#' blocks and medoid spectra. Each block is represented by a real pixel spectrum
#' (its medoid) rather than by a coordinate-wise aggregated spectrum.
#'
#' This medoid backend behaves more faithfully on compact structures such as
#' galaxy nuclei, bars, and HII-like knots than block-averaging approaches.
#'
#' Use this function when the exact \code{\link{segment}} workflow would require
#' too much RAM because it must build an all-pairs distance object over all
#' valid spaxels.
#'
#' @param input A FITS-like object with an \code{imDat} 3-D numeric cube, or a
#'   raw 3-D numeric array.
#' @param Ncomp Integer, the number of clusters to form. Defaults to `15`.
#' @param redshift Numeric redshift placeholder kept for API compatibility.
#' @param scale_fn Row-wise scaling function. If \code{NULL}, an internal
#'   robust median/MAD scaler is used.
#' @param block_size Spatial block size in pixels. If \code{NULL}, it is chosen
#'   automatically from \code{ram_gb}. When left as \code{NULL},
#'   \code{segment_big_cube()} first tries a memory-saving
#'   \code{fastcluster::hclust.vector()} Ward backend on the full set of valid
#'   pixels for medium-large cubes, then falls back to block medoids only when
#'   the cube is too large for that path.
#' @param ram_gb Approximate RAM budget used when choosing \code{block_size}.
#' @param frac_for_dist Fraction of \code{ram_gb} allocated to the block-level
#'   condensed distance object.
#' @param m_cap Optional cap on the number of block representatives kept for
#'   clustering.
#' @param target_snr Optional minimum accepted SNR per cluster. When supplied,
#'   Capivara chooses the largest number of clusters whose minimum cluster SNR
#'   remains above this threshold.
#' @param var_cube Optional variance cube matching the input cube. Used only
#'   when \code{target_snr} is supplied.
#' @param k_values Optional candidate cluster counts tested when
#'   \code{target_snr} is supplied. If omitted, the scalable backend screens up
#'   to the largest 50 cuts.
#' @param wavelength_range Optional wavelength interval used to compute SNR when
#'   \code{target_snr} is supplied.
#' @param valid_mode Either \code{"signal"} or \code{"finite_frac"}.
#' @param min_finite_frac Minimum finite fraction required when
#'   \code{valid_mode = "finite_frac"}.
#' @param use_pca Logical; if \code{TRUE}, run PCA on block medoid spectra
#'   before clustering.
#' @param pca_k Number of PCA components if \code{use_pca = TRUE}.
#' @param seed Integer random seed used when boundary refinement must sample a
#'   subset of candidate pixels.
#' @param dist_method Either \code{"euclidean"} or \code{"capivara_l1"} for
#'   block-level clustering. The default \code{"euclidean"} is the most
#'   coherent choice with Ward's D2 linkage in the scalable backend.
#' @param pixel_assign Either \code{"l2_nearest_center"} or
#'   \code{"l1_nearest_center"} for full-resolution reassignment.
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
#' @param polish_iters Number of medoid-polishing iterations in full pixel
#'   space after the initial assignment.
#' @param refine Logical; if \code{TRUE}, only boundary pixels are reassigned
#'   after medoid polishing.
#' @param refine_radius Neighborhood radius in pixels used to detect boundary
#'   pixels.
#' @param refine_max_pixels Maximum number of boundary pixels to refine.
#' @param verbose Logical.
#' @return A list containing the full-resolution \code{cluster_map}, the
#'   chosen \code{Ncomp}, a per-cluster SNR summary, axis/header data, the
#'   original cube, and \code{backend} indicating whether the scalable result
#'   used the exact-ish vector Ward backend or the block-medoid approximation.
#'   Scalable-engine diagnostics are stored under \code{backend_info}. When
#'   \code{use_starlet_mask = TRUE}, starlet products are stored under
#'   \code{starlet_info}.
#' @seealso \code{\link{segment}}, \code{\link{build_starlet_mask}}
#' @export
segment_big_cube <- function(input,
                             Ncomp = 15,
                             redshift = 0,
                             scale_fn = NULL,
                             block_size = NULL,
                             ram_gb = 32,
                             frac_for_dist = 0.25,
                             m_cap = NULL,
                             target_snr = NULL,
                             var_cube = NULL,
                             k_values = NULL,
                             wavelength_range = NULL,
                             valid_mode = c("signal", "finite_frac"),
                             min_finite_frac = 0.80,
                             use_pca = FALSE,
                             pca_k = 30,
                             seed = 42,
                             dist_method = c("euclidean", "capivara_l1"),
                             pixel_assign = c("l2_nearest_center", "l1_nearest_center"),
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
                             polish_iters = 0L,
                             refine = TRUE,
                             refine_radius = 1L,
                             refine_max_pixels = 50000L,
                             verbose = TRUE) {
  valid_mode <- match.arg(valid_mode)
  dist_method <- match.arg(dist_method)
  pixel_assign <- match.arg(pixel_assign)
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
  stopifnot(!is.null(cube), length(dim(cube)) == 3L)

  n_row <- dim(cube)[1]
  n_col <- dim(cube)[2]
  n_wave <- dim(cube)[3]

  if (is.null(scale_fn)) {
    scale_fn <- .segment_bigcube_safe_scale
  }
  block_size_user_supplied <- !is.null(block_size)

  IFU2D <- cube_to_matrix(cubedat)

  if (valid_mode == "signal") {
    signal <- rowSums(IFU2D, na.rm = TRUE)
    signal[!is.finite(signal)] <- 0
    finite_ok <- signal > 0
  } else {
    n_finite <- rowSums(is.finite(IFU2D))
    finite_ok <- n_finite >= max(10L, floor(min_finite_frac * n_wave))
  }

  valid_pix <- which(finite_ok)
  if (!length(valid_pix)) {
    stop("No valid pixels after filtering.")
  }
  if (is.null(target_snr) && !is.null(Ncomp) && length(valid_pix) < Ncomp) {
    stop("`Ncomp` is larger than the number of valid pixels.")
  }

  X_valid <- IFU2D[valid_pix, , drop = FALSE]
  X_pix <- .scale_rows(X_valid, scale_fn = scale_fn, na_to_zero = TRUE)
  flux_mat_valid <- IFU2D[valid_pix, , drop = FALSE]

  center_metric <- if (pixel_assign == "l2_nearest_center") "l2" else "l1"
  wavelengths <- if (!is.null(cubedat$axDat)) {
    FITSio::axVec(3, cubedat$axDat)
  } else {
    seq_len(ncol(flux_mat_valid))
  }

  var_mat_valid <- NULL
  if (!is.null(target_snr)) {
    if (is.null(var_cube)) {
      var_mat_valid <- pmax(flux_mat_valid, 0)
    } else {
      var_input <- .as_cubedat(var_cube)
      if (!identical(dim(var_input$imDat), dim(cubedat$imDat))) {
        stop("`var_cube` must have the same dimensions as the input cube.")
      }
      var_mat_valid <- cube_to_matrix(var_input)[valid_pix, , drop = FALSE]
    }
    var_mat_valid[!is.finite(var_mat_valid) | var_mat_valid < 0] <- NA_real_
  }

  vector_cap <- .segment_bigcube_vector_cap(ram_gb = ram_gb)
  use_vector_backend <- !block_size_user_supplied && length(valid_pix) <= vector_cap

  block_map <- NULL
  n_blocks <- NA_integer_
  valid_blocks <- integer()
  bl <- integer()
  block_medoid_pos <- integer()
  feat_blocks <- NULL
  pca <- NULL
  centers <- NULL

  if (use_vector_backend) {
    feat_pixels <- X_pix

    if (use_pca) {
      if (!requireNamespace("irlba", quietly = TRUE)) {
        stop("Package 'irlba' required for PCA. Install it or set use_pca = FALSE.")
      }
      pca_k_use <- min(pca_k, ncol(feat_pixels), max(1L, nrow(feat_pixels) - 1L))
      if (pca_k_use < 1L) {
        stop("Not enough valid pixels to compute PCA for `segment_big_cube()`.")
      }
      if (verbose) {
        message(sprintf("PCA with k=%d on full valid pixels...", pca_k_use))
      }
      pca <- irlba::prcomp_irlba(feat_pixels, n = pca_k_use, center = FALSE, scale. = FALSE)
      feat_pixels <- pca$x
    }

    vector_fit <- .segment_bigcube_vector_cluster(
      feat_pixels = feat_pixels,
      valid_pix = valid_pix,
      n_row = n_row,
      n_col = n_col,
      Ncomp = Ncomp,
      target_snr = target_snr,
      flux_mat_valid = flux_mat_valid,
      var_mat_valid = var_mat_valid,
      wavelengths = wavelengths,
      wavelength_range = wavelength_range,
      k_values = k_values,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation,
      verbose = verbose
    )

    cluster_map <- vector_fit$cluster_map
    pix_lbl <- vector_fit$pixel_labels
    hc <- vector_fit$hclust
    chosen_k <- vector_fit$Ncomp
    snr_grid <- vector_fit$snr_grid
    backend_name <- "vector_ward"
  } else {
    if (is.null(block_size)) {
      sel <- .segment_bigcube_choose_block_size(
        n_row = n_row,
        n_col = n_col,
        ram_gb = ram_gb,
        frac_for_dist = frac_for_dist,
        m_cap = m_cap
      )
      block_size <- sel$block_size
      if (verbose) {
        message(sprintf(
          "Chosen block_size=%d (target blocks<=%d, expected n_blocks≈%d)",
          block_size, sel$m_target, sel$n_blocks
        ))
      }
    } else if (verbose) {
      message(sprintf("Using user block_size=%d", block_size))
    }

    br <- ceiling(seq_len(n_row) / block_size)
    bc <- ceiling(seq_len(n_col) / block_size)
    nb_r <- max(br)
    nb_c <- max(bc)
    n_blocks <- nb_r * nb_c
    block_map <- outer(br, bc, FUN = function(r, c) (r - 1L) * nb_c + c)

    block_vec <- as.vector(block_map)
    block_vec_valid <- block_vec[valid_pix]
    by_block <- split(seq_along(valid_pix), block_vec_valid)

    block_rep <- matrix(NA_real_, nrow = n_blocks, ncol = ncol(X_pix))
    block_medoid_pos <- rep(NA_integer_, n_blocks)

    if (verbose) {
      message("Selecting medoid representatives per block...")
    }
    for (b in names(by_block)) {
      idx <- by_block[[b]]
      M <- X_pix[idx, , drop = FALSE]
      midx <- .segment_bigcube_medoid_index(M, metric = center_metric)
      b_int <- as.integer(b)
      block_rep[b_int, ] <- M[midx, ]
      block_medoid_pos[b_int] <- idx[midx]
    }

    valid_blocks <- which(!is.na(block_medoid_pos))
    if (!length(valid_blocks)) {
      stop("No valid blocks after medoid selection.")
    }
    if (is.null(target_snr) && !is.null(Ncomp) && Ncomp > length(valid_blocks)) {
      stop("`Ncomp` is larger than the number of valid block medoids.")
    }

    feat_blocks <- block_rep[valid_blocks, , drop = FALSE]
    feat_pixels <- X_pix

    if (use_pca) {
      if (!requireNamespace("irlba", quietly = TRUE)) {
        stop("Package 'irlba' required for PCA. Install it or set use_pca = FALSE.")
      }
      pca_k_use <- min(pca_k, ncol(feat_blocks), max(1L, nrow(feat_blocks) - 1L))
      if (pca_k_use < 1L) {
        stop("Not enough valid blocks to compute PCA for `segment_big_cube()`.")
      }
      if (verbose) {
        message(sprintf("PCA with k=%d on medoid block representatives...", pca_k_use))
      }
      pca <- irlba::prcomp_irlba(feat_blocks, n = pca_k_use, center = FALSE, scale. = FALSE)
      feat_blocks <- pca$x
      feat_pixels <- X_pix %*% pca$rotation
    }

    if (verbose) {
      message(sprintf("Computing distances for %d block medoids...", nrow(feat_blocks)))
    }
    d <- if (dist_method == "capivara_l1") {
      torch_dist(feat_blocks, p = 1)
    } else {
      torch_dist(feat_blocks, p = 2)
    }

    hc <- fastcluster::hclust(d, method = "ward.D2")

    if (is.null(target_snr)) {
      chosen_k <- Ncomp
      bl <- stats::cutree(hc, k = chosen_k)
      centers <- .segment_bigcube_compute_cluster_medoids(
        feat_blocks,
        labels = bl,
        K = chosen_k,
        metric = center_metric
      )
      assigned <- .segment_bigcube_assign_map(
        feat_pixels = feat_pixels,
        centers = centers,
        valid_pix = valid_pix,
        n_row = n_row,
        n_col = n_col,
        center_metric = center_metric,
        polish_iters = polish_iters,
        refine = refine,
        refine_radius = refine_radius,
        refine_max_pixels = refine_max_pixels,
        seed = seed,
        verbose = verbose
      )
      cluster_map <- assigned$cluster_map
      pix_lbl <- assigned$pixel_labels
      centers <- assigned$centers
      snr_grid <- NULL
    } else {
      if (is.null(wavelength_range)) {
        wave_idx <- seq_len(ncol(flux_mat_valid))
      } else {
        if (length(wavelength_range) != 2) {
          stop("`wavelength_range` must have length 2.")
        }
        wave_idx <- which(wavelengths >= min(wavelength_range) & wavelengths <= max(wavelength_range))
        if (!length(wave_idx)) {
          stop("No wavelengths fall inside `wavelength_range`.")
        }
      }

      max_k <- nrow(feat_blocks)
      if (is.null(k_values)) {
        k_values <- seq(from = min(max_k, 50L), to = 1L, by = -1L)
      } else {
        k_values <- sort(unique(as.integer(k_values)), decreasing = TRUE)
        k_values <- k_values[k_values >= 1L & k_values <= max_k]
      }
      if (!length(k_values)) {
        stop("No valid `k_values` to test.")
      }

      snr_grid <- lapply(k_values, function(k) {
        labels_k <- stats::cutree(hc, k = k)
        centers_k <- .segment_bigcube_compute_cluster_medoids(
          feat_blocks,
          labels = labels_k,
          K = k,
          metric = center_metric
        )
        assigned_k <- .segment_bigcube_assign_map(
          feat_pixels = feat_pixels,
          centers = centers_k,
          valid_pix = valid_pix,
          n_row = n_row,
          n_col = n_col,
          center_metric = center_metric,
          polish_iters = polish_iters,
          refine = refine,
          refine_radius = refine_radius,
          refine_max_pixels = refine_max_pixels,
          seed = seed,
          verbose = FALSE
        )
        cluster_snr_k <- .compute_cluster_snr_from_variance(
          flux_mat = flux_mat_valid,
          var_mat = var_mat_valid,
          clusters = assigned_k$pixel_labels,
          wave_idx = wave_idx,
          snr_stat = snr_stat,
          variance_inflation = variance_inflation
        )
        screen <- .evaluate_cluster_snr_screen(cluster_snr_k, target_snr)

        data.frame(
          Ncomp = k,
          min_cluster_snr = screen$min_cluster_snr,
          all_clusters_pass = screen$all_clusters_pass
        )
      })
      snr_grid <- do.call(rbind, snr_grid)

      ok <- which(snr_grid$all_clusters_pass)
      if (!length(ok)) {
        stop("No clustering configuration satisfies the target SNR.")
      }

      chosen_k <- snr_grid$Ncomp[ok[1]]
      bl <- stats::cutree(hc, k = chosen_k)
      centers <- .segment_bigcube_compute_cluster_medoids(
        feat_blocks,
        labels = bl,
        K = chosen_k,
        metric = center_metric
      )
      assigned <- .segment_bigcube_assign_map(
        feat_pixels = feat_pixels,
        centers = centers,
        valid_pix = valid_pix,
        n_row = n_row,
        n_col = n_col,
        center_metric = center_metric,
        polish_iters = polish_iters,
        refine = refine,
        refine_radius = refine_radius,
        refine_max_pixels = refine_max_pixels,
        seed = seed,
        verbose = verbose
      )
      cluster_map <- assigned$cluster_map
      pix_lbl <- assigned$pixel_labels
      centers <- assigned$centers
    }

    backend_name <- "medoid"
  }

  if (is.null(target_snr)) {
    sn <- .compute_signal_noise(IFU2D)
    cluster_snr <- .compute_cluster_snr(
      clusters = pix_lbl,
      signal_valid = sn$signal[valid_pix],
      noise_valid = sn$noise[valid_pix]
    )
  } else {
    cluster_snr <- .compute_cluster_snr_from_variance(
      flux_mat = flux_mat_valid,
      var_mat = var_mat_valid,
      clusters = pix_lbl,
      wave_idx = wave_idx,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation
    )
  }

  out <- list(
    cluster_map = cluster_map,
    cluster_snr = cluster_snr,
    Ncomp = chosen_k,
    snr_grid = snr_grid,
    hclust = hc,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    original_cube = cubedat,
    backend = backend_name
  )

  if (identical(backend_name, "vector_ward")) {
    out$backend_info <- list(
      strategy = "fastcluster::hclust.vector",
      vector_cap = vector_cap,
      pca = pca,
      valid_pix = valid_pix
    )
  } else {
    out$backend_info <- list(
      strategy = "block_medoids",
      block_map = block_map,
      block_size = block_size,
      n_blocks = n_blocks,
      valid_blocks = valid_blocks,
      block_labels_valid = bl,
      block_features = feat_blocks,
      block_medoid_pos = block_medoid_pos,
      pca = pca,
      centers = centers,
      valid_pix = valid_pix
    )
  }

  if (!is.null(starlet_prep$starlet_info)) {
    out$starlet_info <- starlet_prep$starlet_info
  }

  out
}
