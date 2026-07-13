args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/benchmark_white_m83_residual_detectors.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

fits_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "..", "capivara_experimental", "white_m83.fits")
}

out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_residual_detectors")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_bg_path <- file.path(out_dir, "white_m83_background_models.png")
png_det_path <- file.path(out_dir, "white_m83_residual_detectors.png")
csv_summary_path <- file.path(out_dir, "white_m83_residual_detectors_summary.csv")
csv_peaks_path <- file.path(out_dir, "white_m83_residual_detectors_peaks.csv")
rds_path <- file.path(out_dir, "white_m83_residual_detectors_results.rds")

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(imager)
  library(mgcv)
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

fit_gam_background <- function(img_fill,
                               valid_mask,
                               init_bg,
                               max_fit_points = 120000L,
                               k_basis = 180L,
                               mask_sigma = 4) {
  idx <- which(valid_mask, arr.ind = TRUE)
  residual0 <- img_fill - init_bg
  st0 <- estimate_center_sigma(residual0, mask = valid_mask)
  keep <- residual0 <= (st0$center + mask_sigma * st0$sigma)
  keep[!valid_mask] <- FALSE

  idx_keep <- which(keep, arr.ind = TRUE)
  if (nrow(idx_keep) > max_fit_points) {
    set.seed(1)
    idx_keep <- idx_keep[sample.int(nrow(idx_keep), max_fit_points), , drop = FALSE]
  }

  dat <- data.frame(
    x = idx_keep[, 1],
    y = idx_keep[, 2],
    z = img_fill[keep][seq_len(nrow(idx_keep))]
  )
  dat$z <- img_fill[cbind(idx_keep[, 1], idx_keep[, 2])]

  fit <- mgcv::bam(
    z ~ s(x, y, bs = "tp", k = k_basis),
    data = dat,
    method = "fREML",
    discrete = TRUE,
    nthreads = 1
  )

  pred_grid <- expand.grid(
    x = seq_len(nrow(img_fill)),
    y = seq_len(ncol(img_fill))
  )
  pred <- matrix(
    as.numeric(stats::predict(fit, newdata = pred_grid)),
    nrow = nrow(img_fill),
    ncol = ncol(img_fill),
    byrow = FALSE
  )

  list(background = pred, fit = fit, keep_mask = keep)
}

dog_filter <- function(residual_ci, sigma_small, sigma_large) {
  cimg_to_matrix(imager::isoblur(residual_ci, sigma = sigma_small)) -
    cimg_to_matrix(imager::isoblur(residual_ci, sigma = sigma_large))
}

matched_gaussian <- function(residual_ci, sigma) {
  cimg_to_matrix(imager::isoblur(residual_ci, sigma = sigma))
}

log_filter <- function(residual_ci, sigma) {
  blur <- imager::isoblur(residual_ci, sigma = sigma)
  dxx <- cimg_to_matrix(imager::deriche(blur, sigma = sigma, order = 2, axis = "x"))
  dyy <- cimg_to_matrix(imager::deriche(blur, sigma = sigma, order = 2, axis = "y"))
  -(dxx + dyy)
}

x <- FITSio::readFITS(fits_path)
img <- x$imDat
stopifnot(is.matrix(img))

prep <- fill_invalid(img)
img_fill <- prep$image
valid_mask <- prep$mask
ci <- imager::as.cimg(img_fill)

gauss_bg <- cimg_to_matrix(imager::isoblur(ci, sigma = 14))
gauss_resid <- img_fill - gauss_bg

gam_bg_fit <- fit_gam_background(
  img_fill = img_fill,
  valid_mask = valid_mask,
  init_bg = gauss_bg,
  max_fit_points = 120000L,
  k_basis = 180L,
  mask_sigma = 4
)
gam_bg <- gam_bg_fit$background
gam_resid <- img_fill - gam_bg

gauss_resid_ci <- imager::as.cimg(gauss_resid)
gam_resid_ci <- imager::as.cimg(gam_resid)

components <- list(
  raw = mask_like(img, valid_mask),
  gaussian_background = mask_like(gauss_bg, valid_mask),
  gaussian_residual = mask_like(gauss_resid, valid_mask),
  gam_background = mask_like(gam_bg, valid_mask),
  gam_residual = mask_like(gam_resid, valid_mask)
)

methods <- list(
  raw = function() components$raw,
  tophat_sq7 = function() {
    opn <- cimg_to_matrix(opening_square(ci, 7))
    mask_like(img_fill - opn, valid_mask)
  },
  tophat_sq11 = function() {
    opn <- cimg_to_matrix(opening_square(ci, 11))
    mask_like(img_fill - opn, valid_mask)
  },
  gaussbg_match = function() {
    mask_like(matched_gaussian(gauss_resid_ci, sigma = 1.6), valid_mask)
  },
  gaussbg_dog = function() {
    mask_like(dog_filter(gauss_resid_ci, sigma_small = 1.2, sigma_large = 2.8), valid_mask)
  },
  gaussbg_log = function() {
    mask_like(log_filter(gauss_resid_ci, sigma = 1.6), valid_mask)
  },
  gambg_match = function() {
    mask_like(matched_gaussian(gam_resid_ci, sigma = 1.6), valid_mask)
  },
  gambg_dog = function() {
    mask_like(dog_filter(gam_resid_ci, sigma_small = 1.2, sigma_large = 2.8), valid_mask)
  },
  gambg_log = function() {
    mask_like(log_filter(gam_resid_ci, sigma = 1.6), valid_mask)
  }
)

method_labels <- c(
  raw = "Raw",
  tophat_sq7 = "Top-hat (sq=7)",
  tophat_sq11 = "Top-hat (sq=11)",
  gaussbg_match = "Gaussian BG -> Matched",
  gaussbg_dog = "Gaussian BG -> DoG",
  gaussbg_log = "Gaussian BG -> LoG",
  gambg_match = "Spline BG -> Matched",
  gambg_dog = "Spline BG -> DoG",
  gambg_log = "Spline BG -> LoG"
)

results <- vector("list", length(methods))
names(results) <- names(methods)

for (nm in names(methods)) {
  t0 <- Sys.time()
  mat <- methods[[nm]]()
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  scored <- score_method(nm, mat, valid_mask, elapsed_sec = elapsed)
  results[[nm]] <- list(
    image = mat,
    metrics = scored$metrics,
    peaks = scored$peaks
  )
}

summary_df <- do.call(rbind, lapply(results, `[[`, "metrics"))
summary_df$label <- method_labels[summary_df$method]
summary_df <- summary_df[order(-summary_df$balanced_score), , drop = FALSE]
row.names(summary_df) <- NULL

peak_df <- do.call(rbind, lapply(results, `[[`, "peaks"))
peak_df$label <- method_labels[peak_df$method]

write.csv(summary_df, csv_summary_path, row.names = FALSE)
write.csv(peak_df, csv_peaks_path, row.names = FALSE)

bg_df <- do.call(
  rbind,
  list(
    matrix_to_df(components$raw, panel = "Raw", step = 2L),
    matrix_to_df(components$gaussian_background, panel = "Large Gaussian Background", step = 2L),
    matrix_to_df(components$gaussian_residual, panel = "Gaussian Residual", step = 2L),
    matrix_to_df(components$gam_background, panel = "Spline/GAM Background", step = 2L),
    matrix_to_df(components$gam_residual, panel = "Spline/GAM Residual", step = 2L)
  )
)

bg_plot <- ggplot(bg_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  facet_wrap(~panel, ncol = 3) +
  theme_void() +
  theme(
    text = element_text(color = "white"),
    plot.background = element_rect(fill = "black", color = NA),
    strip.background = element_rect(fill = "black", color = NA),
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold", color = "white")
  )

det_methods <- c(
  "tophat_sq7",
  "tophat_sq11",
  "gaussbg_match",
  "gaussbg_dog",
  "gambg_match",
  "gambg_dog",
  "gaussbg_log",
  "gambg_log"
)

det_df <- do.call(
  rbind,
  lapply(det_methods, function(nm) {
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

det_peaks <- do.call(
  rbind,
  lapply(det_methods, function(nm) {
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

det_plot <- ggplot(det_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  facet_wrap(~facet, ncol = 4) +
  theme_void() +
  theme(
    text = element_text(color = "white"),
    plot.background = element_rect(fill = "black", color = NA),
    strip.background = element_rect(fill = "black", color = NA),
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold", color = "white")
  )

if (!is.null(det_peaks) && nrow(det_peaks)) {
  det_plot <- det_plot +
    geom_point(
      data = det_peaks,
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

ggsave(png_bg_path, bg_plot, width = 15, height = 9, dpi = 220, bg = "black")
ggsave(png_det_path, det_plot, width = 18, height = 10, dpi = 220, bg = "black")

saveRDS(
  list(
    fits_path = fits_path,
    summary = summary_df,
    peaks = peak_df,
    components = components,
    results = results
  ),
  rds_path
)

cat("Saved background PNG to:", png_bg_path, "\n")
cat("Saved detector PNG to:", png_det_path, "\n")
cat("Saved summary CSV to:", csv_summary_path, "\n")
cat("Saved peaks CSV to:", csv_peaks_path, "\n")
cat("Saved RDS to:", rds_path, "\n")
cat("\nTop methods by balanced_score:\n")
print(summary_df[, c("label", "balanced_score", "n_peaks_5", "peak100_snr", "bg_sigma", "elapsed_sec")], row.names = FALSE)
