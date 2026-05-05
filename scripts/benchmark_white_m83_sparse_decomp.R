args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/benchmark_white_m83_sparse_decomp.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

fits_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "..", "capivara_experimental", "white_m83.fits")
}

out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_sparse_decomp")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_detect_path <- file.path(out_dir, "white_m83_sparse_detection.png")
png_comp_path <- file.path(out_dir, "white_m83_sparse_components.png")
csv_summary_path <- file.path(out_dir, "white_m83_sparse_summary.csv")
csv_peaks_path <- file.path(out_dir, "white_m83_sparse_top_peaks.csv")
rds_path <- file.path(out_dir, "white_m83_sparse_results.rds")

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(imager)
  library(irlba)
  library(reshape2)
  library(viridis)
})

quant_safely <- function(x, p, default = NA_real_) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(default)
  }
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 8, na.rm = TRUE))
}

normalize_panel <- function(v) {
  s <- asinh(v)
  lo <- quant_safely(s, 0.02, default = 0)
  hi <- quant_safely(s, 0.998, default = 1)
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
    return(rep(0, length(v)))
  }
  pmin(pmax((s - lo) / (hi - lo), 0), 1)
}

downsample_matrix <- function(mat, step = 2L) {
  mat[seq(1, nrow(mat), by = step), seq(1, ncol(mat), by = step), drop = FALSE]
}

matrix_to_df <- function(mat, panel, step = 2L) {
  m <- downsample_matrix(mat, step = step)
  df <- reshape2::melt(m, varnames = c("Row", "Col"), value.name = "value")
  df$Row <- (df$Row - 1) * step + 1
  df$Col <- (df$Col - 1) * step + 1
  df$panel <- panel
  df$value_norm <- normalize_panel(df$value)
  df
}

fill_invalid <- function(img) {
  mask <- is.finite(img)
  fill <- stats::median(img[mask], na.rm = TRUE)
  out <- img
  out[!mask] <- fill
  list(image = out, mask = mask, fill = fill)
}

cimg_to_matrix <- function(ci) {
  arr <- as.array(ci)
  if (length(dim(arr)) == 2L) {
    return(arr)
  }
  if (length(dim(arr)) == 4L) {
    return(arr[, , 1, 1, drop = TRUE])
  }
  stop("Unsupported cimg dimensionality.")
}

opening_square <- function(ci, size) {
  imager::dilate_square(imager::erode_square(ci, size), size)
}

mask_like <- function(mat, mask) {
  out <- mat
  out[!mask] <- NA_real_
  out
}

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

top_peak_table <- function(method, mat, snr, peak_mask, n_top = 200L) {
  idx <- which(peak_mask, arr.ind = TRUE)
  if (!nrow(idx)) {
    return(data.frame(
      method = character(),
      row = integer(),
      col = integer(),
      value = double(),
      snr = double()
    ))
  }

  tab <- data.frame(
    method = method,
    row = idx[, 1],
    col = idx[, 2],
    value = mat[peak_mask],
    snr = snr[peak_mask]
  )
  tab <- tab[order(-tab$snr, -tab$value), , drop = FALSE]
  utils::head(tab, n_top)
}

score_method <- function(name, mat, mask, elapsed_sec) {
  stats <- estimate_center_sigma(mat, mask = mask)
  snr <- (mat - stats$center) / stats$sigma
  peaks5 <- local_maxima_mask(snr, mask = mask, threshold = 5)
  peak_values <- sort(as.numeric(snr[peaks5]), decreasing = TRUE)

  peak10 <- if (length(peak_values) >= 10L) peak_values[[10L]] else NA_real_
  peak50 <- if (length(peak_values) >= 50L) peak_values[[50L]] else NA_real_
  peak100 <- if (length(peak_values) >= 100L) peak_values[[100L]] else NA_real_
  peak200 <- if (length(peak_values) >= 200L) peak_values[[200L]] else NA_real_

  metrics <- data.frame(
    method = name,
    bg_center = stats$center,
    bg_sigma = stats$sigma,
    peak10_snr = peak10,
    peak50_snr = peak50,
    peak100_snr = peak100,
    peak200_snr = peak200,
    n_peaks_5 = sum(peak_values > 5),
    n_peaks_7 = sum(peak_values > 7),
    elapsed_sec = elapsed_sec,
    stringsAsFactors = FALSE
  )

  metrics$balanced_score <- with(
    metrics,
    pmin(ifelse(is.finite(peak100_snr), peak100_snr, 0), 100) /
      (log1p(n_peaks_5)^2)
  )

  list(metrics = metrics, peaks = top_peak_table(name, mat, snr, peaks5, n_top = 200L))
}

rank_k_background <- function(mat, k) {
  s <- irlba::irlba(mat, nv = k, nu = k)
  s$u %*% (s$d * t(s$v))
}

soft_thresh <- function(x, tau) {
  sign(x) * pmax(abs(x) - tau, 0)
}

svt_full <- function(X, tau) {
  s <- svd(X)
  d <- pmax(s$d - tau, 0)
  r <- sum(d > 0)
  if (r == 0L) {
    return(list(X = matrix(0, nrow(X), ncol(X)), rank = 0L))
  }
  list(
    X = s$u[, seq_len(r), drop = FALSE] %*% (d[seq_len(r)] * t(s$v[, seq_len(r), drop = FALSE])),
    rank = r
  )
}

rpca_ialm <- function(D,
                      lambda = 1 / sqrt(max(dim(D))),
                      max_iter = 20L,
                      tol = 1e-5) {
  m <- nrow(D)
  n <- ncol(D)
  norm_two <- svd(D, nu = 0, nv = 0)$d[1]
  norm_inf <- max(abs(D)) / lambda
  dual_norm <- max(norm_two, norm_inf)
  Y <- D / dual_norm
  mu <- 1.25 / norm_two
  mu_bar <- mu * 1e7
  rho <- 1.5
  L <- matrix(0, m, n)
  S <- matrix(0, m, n)
  dnorm <- norm(D, "F")
  last_rank <- 0L

  for (iter in seq_len(max_iter)) {
    sv <- svt_full(D - S + (1 / mu) * Y, 1 / mu)
    L <- sv$X
    S <- soft_thresh(D - L + (1 / mu) * Y, lambda / mu)
    Z <- D - L - S
    Y <- Y + mu * Z
    mu <- min(mu * rho, mu_bar)
    err <- norm(Z, "F") / dnorm
    last_rank <- sv$rank
    if (err < tol) {
      break
    }
  }

  list(L = L, S = S, rank = last_rank)
}

x <- FITSio::readFITS(fits_path)
img <- x$imDat
stopifnot(is.matrix(img))

prep <- fill_invalid(img)
img_fill <- prep$image
valid_mask <- prep$mask
img_centered <- img_fill - prep$fill
ci <- imager::as.cimg(img_fill)

components <- list()
results <- list()

run_method <- function(name, mat) {
  t0 <- Sys.time()
  out <- mat()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  scored <- score_method(name, out, valid_mask, elapsed_sec = elapsed)
  list(image = out, metrics = scored$metrics, peaks = scored$peaks, elapsed = elapsed)
}

components$raw <- mask_like(img, valid_mask)

results$raw <- run_method("raw", function() {
  components$raw
})

results$tophat_sq7 <- run_method("tophat_sq7", function() {
  opn <- cimg_to_matrix(opening_square(ci, 7))
  residual <- mask_like(img_fill - opn, valid_mask)
  components$tophat_sq7_background <<- mask_like(opn, valid_mask)
  residual
})

results$tophat_sq11 <- run_method("tophat_sq11", function() {
  opn <- cimg_to_matrix(opening_square(ci, 11))
  residual <- mask_like(img_fill - opn, valid_mask)
  components$tophat_sq11_background <<- mask_like(opn, valid_mask)
  residual
})

results$svd_rank8 <- run_method("svd_rank8", function() {
  bg <- rank_k_background(img_centered, 8) + prep$fill
  residual <- img_fill - bg
  components$svd_rank8_background <<- mask_like(bg, valid_mask)
  mask_like(residual, valid_mask)
})

results$svd_rank20 <- run_method("svd_rank20", function() {
  bg <- rank_k_background(img_centered, 20) + prep$fill
  residual <- img_fill - bg
  components$svd_rank20_background <<- mask_like(bg, valid_mask)
  mask_like(residual, valid_mask)
})

t0_rpca <- Sys.time()
rpca_out <- rpca_ialm(img_centered, max_iter = 20L, tol = 1e-5)
rpca_elapsed <- as.numeric(difftime(Sys.time(), t0_rpca, units = "secs"))
components$rpca_background <- mask_like(rpca_out$L + prep$fill, valid_mask)
results$rpca_sparse <- list(
  image = mask_like(rpca_out$S, valid_mask),
  elapsed = rpca_elapsed
)
rpca_scored <- score_method("rpca_sparse", results$rpca_sparse$image, valid_mask, elapsed_sec = rpca_elapsed)
results$rpca_sparse$metrics <- rpca_scored$metrics
results$rpca_sparse$peaks <- rpca_scored$peaks

method_labels <- c(
  raw = "Raw",
  tophat_sq7 = "Top-hat (sq=7)",
  tophat_sq11 = "Top-hat (sq=11)",
  svd_rank8 = "SVD Residual (rank=8)",
  svd_rank20 = "SVD Residual (rank=20)",
  rpca_sparse = sprintf("RPCA Sparse (rank=%d)", rpca_out$rank)
)

summary_df <- do.call(rbind, lapply(results, `[[`, "metrics"))
summary_df$label <- method_labels[summary_df$method]
summary_df <- summary_df[order(-summary_df$balanced_score), , drop = FALSE]
row.names(summary_df) <- NULL

peak_df <- do.call(rbind, lapply(results, `[[`, "peaks"))
peak_df$label <- method_labels[peak_df$method]

write.csv(summary_df, csv_summary_path, row.names = FALSE)
write.csv(peak_df, csv_peaks_path, row.names = FALSE)

detect_methods <- c("raw", "tophat_sq7", "tophat_sq11", "svd_rank8", "svd_rank20", "rpca_sparse")
detect_df <- do.call(
  rbind,
  lapply(detect_methods, function(nm) {
    df <- matrix_to_df(results[[nm]]$image, panel = method_labels[[nm]], step = 2L)
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    df$facet <- sprintf(
      "%s\nbalanced=%.2f | peaks>5sigma=%d",
      method_labels[[nm]],
      met$balanced_score,
      met$n_peaks_5
    )
    df
  })
)

detect_peaks <- do.call(
  rbind,
  lapply(detect_methods, function(nm) {
    pk <- utils::head(results[[nm]]$peaks, 40L)
    if (!nrow(pk)) {
      return(NULL)
    }
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    pk$facet <- sprintf(
      "%s\nbalanced=%.2f | peaks>5sigma=%d",
      method_labels[[nm]],
      met$balanced_score,
      met$n_peaks_5
    )
    pk
  })
)

detect_plot <- ggplot(detect_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  facet_wrap(~facet, ncol = 3) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold")
  )

if (!is.null(detect_peaks) && nrow(detect_peaks)) {
  detect_plot <- detect_plot +
    geom_point(
      data = detect_peaks,
      aes(x = row, y = col),
      inherit.aes = FALSE,
      shape = 21,
      size = 0.9,
      stroke = 0.15,
      fill = "#7ce0ff",
      color = "white",
      alpha = 0.7
    )
}

comp_df <- do.call(
  rbind,
  list(
    matrix_to_df(components$raw, panel = "Raw", step = 2L),
    matrix_to_df(components$svd_rank20_background, panel = "SVD Background (rank=20)", step = 2L),
    matrix_to_df(components$rpca_background, panel = sprintf("RPCA Background (rank=%d)", rpca_out$rank), step = 2L),
    matrix_to_df(results$rpca_sparse$image, panel = "RPCA Sparse Residual", step = 2L)
  )
)

comp_plot <- ggplot(comp_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  facet_wrap(~panel, ncol = 2) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold")
  )

ggsave(png_detect_path, detect_plot, width = 16, height = 10.5, dpi = 220, bg = "black")
ggsave(png_comp_path, comp_plot, width = 13, height = 10, dpi = 220, bg = "black")

saveRDS(
  list(
    fits_path = fits_path,
    summary = summary_df,
    peaks = peak_df,
    components = components,
    results = results,
    rpca_rank = rpca_out$rank
  ),
  rds_path
)

cat("Saved detection PNG to:", png_detect_path, "\n")
cat("Saved components PNG to:", png_comp_path, "\n")
cat("Saved summary CSV to:", csv_summary_path, "\n")
cat("Saved peaks CSV to:", csv_peaks_path, "\n")
cat("Saved RDS to:", rds_path, "\n")
cat("\nTop methods by balanced_score:\n")
print(summary_df[, c("label", "balanced_score", "n_peaks_5", "peak100_snr", "bg_sigma", "elapsed_sec")], row.names = FALSE)
