args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/export_white_m83_tophat_candidate_fits.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

fits_in <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "..", "capivara_experimental", "white_m83.fits")
}

results_rds <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_denoise_results.rds")
}

out_dir <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

resid_fits <- file.path(out_dir, "white_m83_tophat_sq7_residual.fits")
mask_fits <- file.path(out_dir, "white_m83_tophat_consensus_mask.fits")
masked_resid_fits <- file.path(out_dir, "white_m83_tophat_consensus_residual_only.fits")
masked_white_fits <- file.path(out_dir, "white_m83_tophat_consensus_white_only.fits")

suppressPackageStartupMessages({
  library(FITSio)
  library(RANN)
})

estimate_center_sigma <- function(mat, mask = is.finite(mat)) {
  vals <- as.numeric(mat[mask & is.finite(mat)])
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(list(center = 0, sigma = 1))
  }

  center <- stats::median(vals, na.rm = TRUE)
  neg <- vals[vals <= center]
  sigma <- stats::mad(neg, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- stats::mad(vals, center = center, constant = 1.4826, na.rm = TRUE)
  }
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- stats::sd(vals, na.rm = TRUE)
  }
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- 1
  }

  list(center = center, sigma = sigma)
}

shift_matrix <- function(mat, dx = 0L, dy = 0L, fill = -Inf) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  out <- matrix(fill, nrow = nr, ncol = nc)

  src_r <- seq_len(nr)
  src_c <- seq_len(nc)
  dst_r <- src_r + dx
  dst_c <- src_c + dy

  keep <- dst_r >= 1L & dst_r <= nr
  src_r <- src_r[keep]
  dst_r <- dst_r[keep]

  keep <- dst_c >= 1L & dst_c <= nc
  src_c <- src_c[keep]
  dst_c <- dst_c[keep]

  out[dst_r, dst_c] <- mat[src_r, src_c]
  out
}

local_maxima_mask <- function(mat, mask = is.finite(mat), threshold = -Inf) {
  x <- mat
  x[!mask | !is.finite(x)] <- -Inf
  neigh_max <- matrix(-Inf, nrow(x), ncol(x))

  for (dx in -1:1) {
    for (dy in -1:1) {
      neigh_max <- pmax(neigh_max, shift_matrix(x, dx = dx, dy = dy), na.rm = TRUE)
    }
  }

  mask & is.finite(mat) & (mat == neigh_max) & (mat > threshold)
}

extract_peaks <- function(mat, method, threshold = 5) {
  mask <- is.finite(mat)
  stats <- estimate_center_sigma(mat, mask = mask)
  snr <- (mat - stats$center) / stats$sigma
  peak_mask <- local_maxima_mask(snr, mask = mask, threshold = threshold)
  idx <- which(peak_mask, arr.ind = TRUE)

  peaks <- if (nrow(idx)) {
    data.frame(
      method = method,
      row = idx[, 1],
      col = idx[, 2],
      value = mat[peak_mask],
      snr = snr[peak_mask],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      method = character(),
      row = integer(),
      col = integer(),
      value = double(),
      snr = double(),
      stringsAsFactors = FALSE
    )
  }

  list(peaks = peaks, snr = snr)
}

match_radius <- function(query_df, ref_df, radius = 3) {
  if (!nrow(query_df) || !nrow(ref_df)) {
    return(list(match = rep(FALSE, nrow(query_df)), dist = rep(Inf, nrow(query_df))))
  }

  qmat <- as.matrix(query_df[, c("row", "col")])
  rmat <- as.matrix(ref_df[, c("row", "col")])
  nn <- RANN::nn2(data = rmat, query = qmat, k = 1)
  list(match = nn$nn.dists[, 1] <= radius, dist = nn$nn.dists[, 1])
}

patch_mask_from_peaks <- function(dim_xy, peaks, radius = 2L) {
  out <- matrix(0, nrow = dim_xy[1], ncol = dim_xy[2])
  if (!nrow(peaks)) {
    return(out)
  }

  for (i in seq_len(nrow(peaks))) {
    rr <- seq.int(max(1L, peaks$row[i] - radius), min(dim_xy[1], peaks$row[i] + radius))
    cc <- seq.int(max(1L, peaks$col[i] - radius), min(dim_xy[2], peaks$col[i] + radius))
    out[rr, cc] <- 1
  }
  out
}

orig <- FITSio::readFITS(fits_in)
x <- readRDS(results_rds)

t7 <- x$results$tophat_sq7$image
t11 <- x$results$tophat_sq11$image

t7_info <- extract_peaks(t7, "tophat_sq7", threshold = 5)
t11_info <- extract_peaks(t11, "tophat_sq11", threshold = 5)

match11 <- match_radius(t7_info$peaks, t11_info$peaks, radius = 3)
cand <- t7_info$peaks[match11$match, , drop = FALSE]
cand <- cand[order(-cand$snr), , drop = FALSE]

candidate_mask <- patch_mask_from_peaks(dim(t7), cand, radius = 2L)
valid_mask <- is.finite(orig$imDat)
candidate_mask[!valid_mask] <- 0

masked_resid <- t7
masked_resid[candidate_mask == 0] <- NA_real_

masked_white <- orig$imDat
masked_white[candidate_mask == 0] <- NA_real_

FITSio::writeFITSim(t7, file = resid_fits, type = "double", axDat = orig$axDat)
FITSio::writeFITSim(candidate_mask, file = mask_fits, type = "double", axDat = orig$axDat)
FITSio::writeFITSim(masked_resid, file = masked_resid_fits, type = "double", axDat = orig$axDat)
FITSio::writeFITSim(masked_white, file = masked_white_fits, type = "double", axDat = orig$axDat)

cat("Saved residual FITS to:", resid_fits, "\n")
cat("Saved mask FITS to:", mask_fits, "\n")
cat("Saved masked residual FITS to:", masked_resid_fits, "\n")
cat("Saved masked white FITS to:", masked_white_fits, "\n")
cat(sprintf("Consensus candidates kept: %d\n", nrow(cand)))
