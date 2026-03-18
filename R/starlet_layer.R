#' B3-spline low-pass kernel (1D)
#'
#' @return Numeric vector of length 5.
#' @keywords internal
#' @noRd
b3_kernel <- function() {
  c(1, 4, 6, 4, 1) / 16
}

#' Insert zeros between kernel taps
#'
#' @param k Base kernel.
#' @param step Integer upsampling factor.
#'
#' @return Upsampled kernel.
#' @keywords internal
#' @noRd
upsample_kernel <- function(k, step) {
  if (step <= 1) {
    return(k)
  }

  out <- rep(0, length(k) + (length(k) - 1) * (step - 1))
  out[seq(1, length(out), by = step)] <- k
  out
}

.reflect_pad_vec <- function(v, pad) {
  n <- length(v)
  if (pad <= 0L || n == 0L) {
    return(v)
  }

  pad <- max(0L, min(pad, n - 1L))
  if (pad == 0L) {
    return(v)
  }

  left <- rev(v[seq_len(pad)])
  right <- rev(v[seq.int(n - pad + 1L, n)])
  c(left, v, right)
}

.conv_reflect_rows <- function(mat, k) {
  p <- floor(length(k) / 2)
  out <- matrix(0, nrow(mat), ncol(mat))

  for (r in seq_len(nrow(mat))) {
    pv <- .reflect_pad_vec(mat[r, ], p)
    for (x in seq_len(ncol(mat))) {
      out[r, x] <- sum(k * pv[x:(x + length(k) - 1L)])
    }
  }

  out
}

.conv_reflect_cols <- function(mat, k) {
  p <- floor(length(k) / 2)
  out <- matrix(0, nrow(mat), ncol(mat))

  for (c in seq_len(ncol(mat))) {
    pv <- .reflect_pad_vec(mat[, c], p)
    for (y in seq_len(nrow(mat))) {
      out[y, c] <- sum(k * pv[y:(y + length(k) - 1L)])
    }
  }

  out
}

.smooth_sep_raw <- function(img, k) {
  tmp <- .conv_reflect_rows(img, k)
  .conv_reflect_cols(tmp, k)
}

.smooth_sep_masked <- function(img, k, mask) {
  num <- .smooth_sep_raw(img * mask, k)
  den <- .smooth_sep_raw(mask, k)
  out <- num / pmax(den, 1e-12)
  out[den < 1e-12] <- NA_real_
  out
}

.soft_thresh <- function(x, t) {
  sign(x) * pmax(abs(x) - t, 0)
}

.estimate_sigma <- function(w) {
  stats::mad(w, center = 0, constant = 1, na.rm = TRUE)
}

#' Robust per-wavelength background and scatter
#'
#' @param cube 3-D numeric array with dimensions \code{[nx, ny, nlambda]}.
#'
#' @return A list with numeric vectors \code{bkg} and \code{mad}.
#' @keywords internal
#' @noRd
estimate_bkg_mad_per_lambda <- function(cube) {
  stopifnot(length(dim(cube)) == 3L)

  M <- matrix(cube, nrow = dim(cube)[1] * dim(cube)[2], ncol = dim(cube)[3])
  bkg <- matrixStats::colMedians(M, na.rm = TRUE)
  mad <- matrixStats::colMads(M, na.rm = TRUE)

  repl <- stats::median(mad[is.finite(mad) & mad > 0], na.rm = TRUE)
  if (!is.finite(repl) || repl <= 0) {
    repl <- 1
  }
  mad[!is.finite(mad) | mad <= 0] <- repl

  list(bkg = bkg, mad = mad)
}

#' Collapse a spectral cube into a Sagui-style white-light image
#'
#' This reproduces the white-light collapse used by \code{sagui} before the
#' starlet masking step.
#'
#' @param cube 3-D numeric array with dimensions \code{[nx, ny, nlambda]}.
#' @param kclip Numeric clipping threshold in MAD units.
#' @param use_weights Logical; if \code{TRUE}, use inverse-variance weights.
#'
#' @return A 2-D numeric matrix.
#' @keywords internal
#' @noRd
collapse_white_light <- function(cube, kclip = 2, use_weights = TRUE) {
  stopifnot(length(dim(cube)) == 3L)

  est <- estimate_bkg_mad_per_lambda(cube)
  bkg <- est$bkg
  mad <- est$mad

  for (k in seq_len(dim(cube)[3])) {
    cube[, , k] <- cube[, , k] - bkg[k]
  }

  for (k in seq_len(dim(cube)[3])) {
    thr <- -kclip * mad[k]
    sli <- cube[, , k]
    sli[!is.finite(sli)] <- 0
    sli[sli < thr] <- thr
    cube[, , k] <- sli
  }

  if (isTRUE(use_weights)) {
    w <- 1 / (mad^2)
    w <- w / stats::median(w[is.finite(w) & w > 0], na.rm = TRUE)
  } else {
    w <- rep(1, length(mad))
  }

  wl <- as.numeric(matrix(cube, nrow = dim(cube)[1] * dim(cube)[2], ncol = dim(cube)[3]) %*% w)
  wl <- matrix(wl, nrow = dim(cube)[1], ncol = dim(cube)[2])
  wl[wl < 0] <- 0
  wl
}

#' Starlet (a trous) decomposition
#'
#' @param img 2-D numeric matrix.
#' @param J Number of starlet scales.
#'
#' @return An object of class \code{capivara_starlet}.
#' @keywords internal
#' @noRd
starlet_mask <- function(img, J = 5) {
  stopifnot(is.matrix(img))

  mask_valid <- is.finite(img)
  img0 <- img
  img0[!mask_valid] <- 0

  k0 <- b3_kernel()
  c_j <- img0
  wlist <- vector("list", J)

  for (j in seq_len(J)) {
    step <- 2^(j - 1)
    k_j <- upsample_kernel(k0, step)
    smooth <- .smooth_sep_masked(c_j, k_j, mask_valid)
    wlist[[j]] <- c_j - smooth
    c_j <- smooth
  }

  for (j in seq_len(J)) {
    wlist[[j]][!mask_valid] <- NA_real_
  }
  c_j[!mask_valid] <- NA_real_

  structure(
    list(w = wlist, cJ = c_j, mask = mask_valid),
    class = c("capivara_starlet", "list")
  )
}

#' Reconstruct selected starlet scales
#'
#' @param dec A \code{capivara_starlet} object.
#' @param keep_scales Integer vector of scales to include.
#' @param include_coarse Logical; include the coarse plane.
#' @param denoise_k Optional MAD-based threshold multiplier.
#' @param mode Thresholding mode.
#' @param na_policy How to handle non-finite coefficients while summing.
#'
#' @return A reconstructed 2-D image.
#' @keywords internal
#' @noRd
starlet_reconstruct <- function(dec,
                                keep_scales = c(2, 3),
                                include_coarse = FALSE,
                                denoise_k = NULL,
                                mode = c("soft", "hard"),
                                na_policy = c("preserve", "zero")) {
  mode <- match.arg(mode)
  na_policy <- match.arg(na_policy)

  wlist <- dec$w
  J <- length(wlist)

  stopifnot(
    is.numeric(keep_scales),
    all(keep_scales >= 1),
    all(keep_scales <= J)
  )

  if (!is.null(denoise_k)) {
    kvec <- if (length(denoise_k) == 1) rep(denoise_k, J) else denoise_k
    stopifnot(length(kvec) == J)

    for (j in seq_len(J)) {
      w <- wlist[[j]]
      idx <- is.finite(w)
      if (any(idx)) {
        sig <- .estimate_sigma(w[idx])
        thr <- kvec[j] * sig
        if (is.finite(thr) && thr > 0) {
          if (mode == "soft") {
            w[idx] <- .soft_thresh(w[idx], thr)
          } else {
            w[idx] <- w[idx] * (abs(w[idx]) >= thr)
          }
        }
      }
      wlist[[j]] <- w
    }
  }

  selected <- wlist[keep_scales]

  if (na_policy == "zero") {
    present <- Reduce(`+`, lapply(selected, function(w) as.numeric(is.finite(w))))
    selected0 <- lapply(selected, function(w) {
      w2 <- w
      w2[!is.finite(w2)] <- 0
      w2
    })
    rec <- Reduce(`+`, selected0)

    if (isTRUE(include_coarse)) {
      cJ <- dec$cJ
      present <- present + as.numeric(is.finite(cJ))
      cJ0 <- cJ
      cJ0[!is.finite(cJ0)] <- 0
      rec <- rec + cJ0
    }

    rec[present == 0] <- NA_real_
  } else {
    rec <- Reduce(`+`, selected)
    if (isTRUE(include_coarse)) {
      rec <- rec + dec$cJ
    }
  }

  if (!is.null(dec$mask)) {
    rec[!dec$mask] <- NA_real_
  }

  rec
}

#' Mask a 3-D cube by a spatial map
#'
#' @param cube 3-D numeric array.
#' @param labels 2-D matrix or logical mask.
#' @param mode Either \code{"zero"} or \code{"na"}.
#'
#' @return Masked cube with the same dimensions as the input.
#' @keywords internal
#' @noRd
mask_cube <- function(cube, labels, mode = c("zero", "na")) {
  mode <- match.arg(mode)
  stopifnot(is.array(cube), length(dim(cube)) == 3L, is.matrix(labels))
  stopifnot(identical(dim(labels), dim(cube)[1:2]))

  keep <- is.finite(labels) & (labels > 0)
  keep3 <- array(keep, dim = dim(cube))

  out <- cube
  if (mode == "zero") {
    out[!keep3] <- 0
  } else {
    storage.mode(out) <- "double"
    out[!keep3] <- NA_real_
  }

  out
}

#' Build a Sagui-style starlet mask from a spectral cube
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array.
#' @param collapse_fn Function used to collapse the cube to a 2-D image.
#' @param starlet_J Number of starlet scales.
#' @param starlet_scales Integer vector of scales kept in the reconstruction.
#' @param include_coarse Logical; include the coarse starlet plane.
#' @param denoise_k Optional denoising threshold in MAD units.
#' @param mode Thresholding mode passed to \code{\link{starlet_reconstruct}}.
#' @param positive_only Logical; keep only positive reconstructed values.
#'
#' @return A list with \code{collapsed}, \code{decomposition},
#'   \code{reconstruction}, and \code{mask}.
#' @export
build_starlet_mask <- function(input,
                               collapse_fn = collapse_white_light,
                               starlet_J = 5,
                               starlet_scales = 2:5,
                               include_coarse = FALSE,
                               denoise_k = 0,
                               mode = c("soft", "hard"),
                               positive_only = TRUE) {
  mode <- match.arg(mode)

  cube <- .as_cubedat(input)$imDat
  stopifnot(is.array(cube), length(dim(cube)) == 3L)

  collapsed <- collapse_fn(cube)
  decomposition <- starlet_mask(collapsed, J = starlet_J)
  reconstruction <- starlet_reconstruct(
    decomposition,
    keep_scales = starlet_scales,
    include_coarse = include_coarse,
    denoise_k = denoise_k,
    mode = mode
  )

  mask <- is.finite(reconstruction)
  if (isTRUE(positive_only)) {
    mask <- mask & reconstruction > 0
  }

  list(
    collapsed = collapsed,
    decomposition = decomposition,
    reconstruction = reconstruction,
    mask = mask
  )
}
