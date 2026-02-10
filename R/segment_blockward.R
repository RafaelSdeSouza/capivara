#' Segment an IFU cube via spatial block aggregation + Ward.D2
#'
#' This is a Ward-like approximation designed for very large IFU cutouts where
#' full pairwise distances over all spaxels are impossible (O(n^2) memory/time).
#' The cube is first aggregated in the *spatial* domain into macro-spaxels
#' (blocks) of size block_size x block_size. A robust representative spectrum is
#' computed per block (median across pixels, ignoring NA/Inf). Ward.D2 clustering
#' is then performed on block spectra (optionally after PCA), and the resulting
#' block labels are expanded back to the original spatial resolution.
#'
#' Optional boundary refinement reassigns only pixels near cluster borders using
#' a Ward-consistent ΔSSE rule, yielding sharper edges while keeping runtime low.
#'
#' @param input A FITS-like object with field \code{imDat} as a 3-D numeric array
#'   [n_row, n_col, n_wave]. Compatible with \code{cube_to_matrix()}.
#' @param Ncomp Integer, number of clusters.
#' @param redshift Numeric, redshift for wavelength correction (placeholder here; not applied).
#' @param scale_fn Function for robust row-wise scaling. If NULL, uses internal \code{safe_scale()}.
#' @param block_size Integer, block size in pixels. If NULL, chosen from RAM via \code{choose_block_size()}.
#' @param ram_gb Numeric, available RAM used for selecting block_size when \code{block_size=NULL}.
#' @param frac_for_dist Fraction of RAM to budget for the condensed distance vector.
#' @param m_cap Practical cap on number of blocks to avoid time blowups; defaults depend on ram_gb.
#' @param agg_fn Aggregation for block spectra: "median" (default) or "mean".
#' @param min_finite_frac Minimum fraction of finite values required per pixel spectrum to be valid.
#'   Used for boundary refinement (and for optional pixel validity masks).
#' @param min_block_finite_frac Minimum fraction of finite values required per *block spectrum*.
#' @param use_pca Logical. If TRUE, run PCA on scaled block spectra before Ward.
#' @param pca_k Integer, number of PCA components if use_pca=TRUE.
#' @param seed Integer random seed (only used if needed; currently not used).
#' @param refine Logical. If TRUE, refine only boundary pixels at full resolution.
#' @param refine_radius Integer, neighborhood radius (in pixels) for detecting boundary pixels.
#' @param refine_max_pixels Integer cap on number of pixels to refine (safety).
#' @param block_assign_space "pca" or "scaled". Space used for refinement distances.
#' @param verbose Logical.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{cluster_map}: [n_row x n_col] integer map (NA for invalid/unassigned).
#'   \item \code{block_map}: [n_row x n_col] integer block id per pixel.
#'   \item \code{block_labels}: integer vector of length n_blocks.
#'   \item \code{block_size}: chosen block size.
#'   \item \code{n_blocks}: number of clustered blocks.
#'   \item \code{pca}: PCA object (if use_pca).
#'   \item \code{centers}: cluster means in the chosen feature space (for refinement).
#'   \item \code{cluster_sizes}: cluster sizes in blocks.
#' }
#'
#' @export
segment_blockward <- function(input,
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
                              # capivara-closer:
                              scale_first = TRUE,
                              # do full-res assignment after block Ward
                              full_res_assign = TRUE,
                              # distance settings (capivara uses L1 via torch_dist)
                              dist_method = c("capivara_l1", "euclidean"),
                              # pixel assignment metric (match L1 if capivara_l1)
                              pixel_assign = c("l1_nearest_center", "l2_nearest_center"),
                              # optional polishing in assignment space
                              polish_iters = 0L,
                              # optional boundary refinement
                              refine = TRUE,
                              refine_radius = 1L,
                              refine_max_pixels = 50000L,
                              verbose = TRUE) {
  
  agg_fn <- match.arg(agg_fn)
  valid_mode <- match.arg(valid_mode)
  dist_method <- match.arg(dist_method)
  pixel_assign <- match.arg(pixel_assign)
  
  cubedat <- input
  stopifnot(!is.null(cubedat$imDat), length(dim(cubedat$imDat)) == 3L)
  n_row  <- dim(cubedat$imDat)[1]
  n_col  <- dim(cubedat$imDat)[2]
  n_wave <- dim(cubedat$imDat)[3]
  
  # ---- helpers -------------------------------------------------------------
  
  safe_scale <- function(v) {
    ok <- is.finite(v)
    vv <- v[ok]
    if (length(vv) < 2) return(rep(0, length(v)))
    med <- stats::median(vv)
    sc  <- stats::mad(vv, center = med, constant = 1)
    if (!is.finite(sc) || sc == 0) sc <- stats::sd(vv)
    if (!is.finite(sc) || sc == 0) return(rep(0, length(v)))
    out <- (v - med) / sc
    out[!is.finite(out)] <- 0
    out
  }
  if (is.null(scale_fn)) scale_fn <- safe_scale
  
  choose_block_size <- function(n_row, n_col, ram_gb = 32, frac_for_dist = 0.25, m_cap = NULL) {
    dist_bytes <- ram_gb * 1e9 * frac_for_dist
    m_budget <- floor(sqrt(dist_bytes / 4))  # dist ~ 4*m^2 bytes
    if (is.null(m_cap)) m_cap <- if (ram_gb <= 64) 15000L else 50000L
    m_target <- min(m_budget, m_cap)
    
    for (b in 2:128) {
      m <- ceiling(n_row / b) * ceiling(n_col / b)
      if (m <= m_target) return(list(block_size = b, n_blocks = m, m_target = m_target))
    }
    b <- 128L
    m <- ceiling(n_row / b) * ceiling(n_col / b)
    list(block_size = b, n_blocks = m, m_target = m_target)
  }
  
  boundary_pixels <- function(lbl_map, r = 1L) {
    nr <- nrow(lbl_map); nc <- ncol(lbl_map)
    out <- matrix(FALSE, nr, nc)
    shifts <- expand.grid(dr = -r:r, dc = -r:r)
    shifts <- shifts[!(shifts$dr == 0 & shifts$dc == 0), , drop = FALSE]
    
    for (k in seq_len(nrow(shifts))) {
      dr <- shifts$dr[k]; dc <- shifts$dc[k]
      r1 <- max(1, 1 + dr); r2 <- min(nr, nr + dr)
      c1 <- max(1, 1 + dc); c2 <- min(nc, nc + dc)
      
      a  <- lbl_map[r1:r2, c1:c2, drop = FALSE]
      b  <- lbl_map[(r1 - dr):(r2 - dr), (c1 - dc):(c2 - dc), drop = FALSE]
      
      diff <- (a != b) & !is.na(a) & !is.na(b)
      out[r1:r2, c1:c2] <- out[r1:r2, c1:c2] | diff
    }
    which(out, arr.ind = TRUE)
  }
  
  # capivara torch_dist: L1 via torch::nnf_pdist(p=1)
  torch_dist_l1 <- function(x) {
    x <- as.matrix(x)
    N <- nrow(x)
    if (N < 2) stop("Input matrix must have at least two rows.")
    x_ten <- torch::torch_tensor(x, dtype = torch::torch_float())
    pd <- torch::nnf_pdist(x_ten, p = 1)
    pd <- as.numeric(pd)
    mat <- matrix(0, nrow = N, ncol = N)
    mat[lower.tri(mat, diag = FALSE)] <- pd
    mat <- mat + t(mat)
    as.dist(mat)
  }
  
  # nearest-center assignment under L1 or L2
  assign_nearest_center <- function(X, centers, metric = c("l1", "l2")) {
    metric <- match.arg(metric)
    n <- nrow(X); K <- nrow(centers)
    out <- integer(n)
    
    block <- 2000L  # keep RAM down
    for (i in seq(1, n, by = block)) {
      j <- min(i + block - 1L, n)
      XX <- X[i:j, , drop = FALSE]
      
      if (metric == "l2") {
        # squared L2: ||x||^2 + ||c||^2 - 2 x·c
        c2 <- rowSums(centers^2)
        x2 <- rowSums(XX^2)
        dots <- XX %*% t(centers)  # (block x K)
        D2 <- matrix(rep(x2, times = K), ncol = K) +
          matrix(rep(c2, each = nrow(XX)), ncol = K) -
          2 * dots
        out[i:j] <- max.col(-D2)
      } else {
        # L1: compute sum(abs(x - c)) for each center (block x K)
        # Do it center-by-center to avoid a 3D array.
        D <- matrix(0, nrow = nrow(XX), ncol = K)
        for (k in seq_len(K)) {
          D[, k] <- rowSums(abs(XX - matrix(centers[k, ], nrow = nrow(XX), ncol = ncol(XX), byrow = TRUE)))
        }
        out[i:j] <- max.col(-D)
      }
    }
    out
  }
  
  recompute_centers <- function(X, labels, K) {
    centers <- matrix(NA_real_, nrow = K, ncol = ncol(X))
    for (k in seq_len(K)) {
      idx <- which(labels == k)
      if (length(idx)) centers[k, ] <- colMeans(X[idx, , drop = FALSE])
    }
    centers
  }
  
  # ---- choose block size ---------------------------------------------------
  if (is.null(block_size)) {
    sel <- choose_block_size(n_row, n_col, ram_gb = ram_gb, frac_for_dist = frac_for_dist, m_cap = m_cap)
    block_size <- sel$block_size
    if (verbose) message(sprintf("Chosen block_size=%d (target blocks<=%d, expected n_blocks≈%d)",
                                 block_size, sel$m_target, sel$n_blocks))
  } else {
    if (verbose) message(sprintf("Using user block_size=%d", block_size))
  }
  
  # ---- block map -----------------------------------------------------------
  br <- ceiling(seq_len(n_row) / block_size)
  bc <- ceiling(seq_len(n_col) / block_size)
  nb_r <- max(br); nb_c <- max(bc)
  n_blocks <- nb_r * nb_c
  
  block_map <- outer(br, bc, FUN = function(r, c) (r - 1L) * nb_c + c)
  stopifnot(max(block_map) == n_blocks)
  
  if (verbose) message(sprintf("Block grid: %d x %d = %d blocks", nb_r, nb_c, n_blocks))
  
  # ---- flatten cube to pixels x wavelength --------------------------------
  IFU2D <- cube_to_matrix(cubedat)
  n_pix <- nrow(IFU2D)
  stopifnot(n_pix == n_row * n_col)
  
  # ---- validity mask -------------------------------------------------------
  if (valid_mode == "signal") {
    signal <- rowSums(IFU2D, na.rm = TRUE)
    signal[!is.finite(signal)] <- 0
    finite_ok <- signal > 0
  } else {
    n_finite <- rowSums(is.finite(IFU2D))
    finite_ok <- n_finite >= max(10L, floor(min_finite_frac * n_wave))
  }
  
  valid_pix <- which(finite_ok)
  if (!length(valid_pix)) stop("No valid pixels after filtering.")
  
  block_vec <- as.vector(block_map)
  block_vec_valid <- block_vec[valid_pix]
  
  # ---- capivara-like scaling (pixel-wise) ---------------------------------
  if (verbose) message("Scaling pixel spectra row-wise (capivara-like)...")
  Xp_scaled <- t(apply(IFU2D[valid_pix, , drop = FALSE], 1, scale_fn))  # (n_valid_pix x n_wave)
  
  # ---- aggregate scaled spectra per block ---------------------------------
  if (verbose) message("Aggregating scaled spectra per block...")
  block_rep <- matrix(NA_real_, nrow = n_blocks, ncol = n_wave)
  
  by_block <- split(seq_along(valid_pix), block_vec_valid)
  
  if (agg_fn == "median") {
    agg_col <- function(M) apply(M, 2, stats::median, na.rm = TRUE)
  } else {
    agg_col <- function(M) colMeans(M, na.rm = TRUE)
  }
  
  for (b in names(by_block)) {
    ids <- by_block[[b]]
    block_rep[as.integer(b), ] <- agg_col(Xp_scaled[ids, , drop = FALSE])
  }
  
  # block validity
  bn_finite <- rowSums(is.finite(block_rep))
  block_ok <- bn_finite >= max(10L, floor(min_block_finite_frac * n_wave))
  valid_blocks <- which(block_ok)
  if (!length(valid_blocks)) stop("No valid blocks after aggregation/filters.")
  
  if (verbose) message(sprintf("Valid blocks: %d / %d (%.1f%%)",
                               length(valid_blocks), n_blocks, 100 * length(valid_blocks) / n_blocks))
  
  Xb <- block_rep[valid_blocks, , drop = FALSE]  # already scaled-space
  
  # ---- optional PCA on blocks ---------------------------------------------
  feat_blocks <- Xb
  pca <- NULL
  if (use_pca) {
    if (!requireNamespace("irlba", quietly = TRUE)) {
      stop("Package 'irlba' required for PCA. Install it or set use_pca=FALSE.")
    }
    if (verbose) message(sprintf("PCA with k=%d (irlba) on block reps...", pca_k))
    pca <- irlba::prcomp_irlba(Xb, n = pca_k, center = FALSE, scale. = FALSE)
    feat_blocks <- pca$x
  }
  
  # ---- Ward clustering on blocks using capivara distance -------------------
  m <- nrow(feat_blocks)
  if (verbose) message(sprintf("Computing distances for %d blocks...", m))
  
  if (dist_method == "capivara_l1") {
    if (!requireNamespace("torch", quietly = TRUE)) stop("Need 'torch' installed for capivara_l1 distances.")
    d <- torch_dist_l1(feat_blocks)   # L1 (Manhattan)
  } else {
    d <- stats::dist(feat_blocks)     # Euclidean
  }
  
  hc <- fastcluster::hclust(d, method = "ward.D2")
  bl <- stats::cutree(hc, k = Ncomp)
  
  # ---- centers for pixel assignment ---------------------------------------
  if (use_pca) {
    # compute centers in PCA space from block features
    centers <- recompute_centers(feat_blocks, bl, Ncomp)
    # project pixels to PCA space if needed
    Zp <- Xp_scaled %*% pca$rotation
    metric <- if (pixel_assign == "l1_nearest_center") "l1" else "l2"
    pix_lbl <- assign_nearest_center(Zp, centers, metric = metric)
    # optional polishing
    if (polish_iters > 0L) {
      if (verbose) message(sprintf("Polishing %d iters in PCA space...", polish_iters))
      for (it in seq_len(polish_iters)) {
        centers <- recompute_centers(Zp, pix_lbl, Ncomp)
        pix_lbl <- assign_nearest_center(Zp, centers, metric = metric)
      }
    }
  } else {
    centers <- recompute_centers(Xb, bl, Ncomp)  # centers in scaled spectra space
    metric <- if (pixel_assign == "l1_nearest_center") "l1" else "l2"
    pix_lbl <- assign_nearest_center(Xp_scaled, centers, metric = metric)
    if (polish_iters > 0L) {
      if (verbose) message(sprintf("Polishing %d iters in scaled space...", polish_iters))
      for (it in seq_len(polish_iters)) {
        centers <- recompute_centers(Xp_scaled, pix_lbl, Ncomp)
        pix_lbl <- assign_nearest_center(Xp_scaled, centers, metric = metric)
      }
    }
  }
  
  # ---- output map ----------------------------------------------------------
  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  if (full_res_assign) {
    cluster_map[valid_pix] <- pix_lbl
  } else {
    # diagnostic block paint
    block_labels <- rep(NA_integer_, n_blocks)
    block_labels[valid_blocks] <- bl
    cluster_map[] <- block_labels[block_map]
  }
  
  # ---- optional boundary refinement ----------------------------------------
  if (refine && full_res_assign) {
    if (verbose) message("Boundary refinement...")
    bp <- boundary_pixels(cluster_map, r = refine_radius)
    if (nrow(bp) > refine_max_pixels) {
      set.seed(seed)
      bp <- bp[sample.int(nrow(bp), refine_max_pixels), , drop = FALSE]
    }
    if (nrow(bp) > 0) {
      lin <- bp[, 1] + (bp[, 2] - 1L) * n_row
      pos <- match(lin, valid_pix)
      ok <- !is.na(pos)
      if (any(ok)) {
        pos <- pos[ok]
        lin <- lin[ok]
        metric <- if (pixel_assign == "l1_nearest_center") "l1" else "l2"
        if (use_pca) {
          Zq <- (Xp_scaled[pos, , drop = FALSE]) %*% pca$rotation
          new_lbl <- assign_nearest_center(Zq, centers, metric = metric)
        } else {
          Xq <- Xp_scaled[pos, , drop = FALSE]
          new_lbl <- assign_nearest_center(Xq, centers, metric = metric)
        }
        cluster_map[lin] <- new_lbl
      }
    }
  }
  
  list(
    cluster_map = cluster_map,
    block_map = block_map,
    block_size = block_size,
    n_blocks = n_blocks,
    valid_blocks = valid_blocks,
    block_labels_valid = bl,
    block_features = feat_blocks,
    pca = pca,
    centers = centers,
    valid_pix = valid_pix,
    header = cubedat$hdr,
    axDat = cubedat$axDat,
    original_cube = cubedat
  )
}

