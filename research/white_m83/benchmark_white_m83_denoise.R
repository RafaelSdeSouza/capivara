args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/benchmark_white_m83_denoise.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

fits_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "..", "capivara_experimental", "white_m83.fits")
}

out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_panel_path <- file.path(out_dir, "white_m83_denoise_panels.png")
png_metric_path <- file.path(out_dir, "white_m83_denoise_metrics.png")
csv_summary_path <- file.path(out_dir, "white_m83_denoise_summary.csv")
csv_peaks_path <- file.path(out_dir, "white_m83_top_peaks.csv")
rds_path <- file.path(out_dir, "white_m83_denoise_results.rds")

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(imager)
  library(patchwork)
  library(reshape2)
  library(viridis)
})

files <- list.files(file.path(repo_root, "R"), full.names = TRUE)
for (f in files) {
  sys.source(f, envir = .GlobalEnv)
}

quant_safely <- function(x, p, default = NA_real_) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(default)
  }
  as.numeric(stats::quantile(x, probs = p, names = FALSE, type = 8, na.rm = TRUE))
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

starlet_enhance <- function(img,
                            keep_scales,
                            denoise_k,
                            J = 5,
                            mode = "soft") {
  dec <- starlet_mask(img, J = J)
  starlet_reconstruct(
    dec,
    keep_scales = keep_scales,
    include_coarse = FALSE,
    denoise_k = denoise_k,
    mode = mode
  )
}

score_method <- function(name, mat, mask, elapsed_sec) {
  stats <- estimate_center_sigma(mat, mask = mask)
  snr <- (mat - stats$center) / stats$sigma
  peaks3 <- local_maxima_mask(snr, mask = mask, threshold = 3)
  peak_values <- sort(as.numeric(snr[peaks3]), decreasing = TRUE)
  peak10 <- if (length(peak_values) >= 10L) peak_values[[10L]] else NA_real_
  peak50 <- if (length(peak_values) >= 50L) peak_values[[50L]] else NA_real_
  peak100 <- if (length(peak_values) >= 100L) peak_values[[100L]] else NA_real_
  peak200 <- if (length(peak_values) >= 200L) peak_values[[200L]] else NA_real_

  top20 <- utils::head(peak_values, 20)
  top50 <- utils::head(peak_values, 50)
  peak_table <- top_peak_table(name, mat, snr, peaks3, n_top = 200L)

  metrics <- data.frame(
    method = name,
    bg_center = stats$center,
    bg_sigma = stats$sigma,
    q99_snr = quant_safely(snr[mask], 0.99, default = 0),
    q999_snr = quant_safely(snr[mask], 0.999, default = 0),
    peak10_snr = peak10,
    peak50_snr = peak50,
    peak100_snr = peak100,
    peak200_snr = peak200,
    peak20_median_snr = quant_safely(top20, 0.5, default = 0),
    peak50_median_snr = quant_safely(top50, 0.5, default = 0),
    peak95_snr = quant_safely(peak_values, 0.95, default = 0),
    n_peaks_3 = sum(peak_values > 3),
    n_peaks_5 = sum(peak_values > 5),
    n_peaks_7 = sum(peak_values > 7),
    neg_q01_snr = quant_safely(snr[mask], 0.01, default = 0),
    elapsed_sec = elapsed_sec,
    stringsAsFactors = FALSE
  )

  metrics$compact_score <- with(
    metrics,
    peak20_median_snr +
      0.25 * peak95_snr +
      0.50 * log1p(n_peaks_5) -
      0.15 * log1p(n_peaks_3) +
      0.10 * q999_snr
  )

  metrics$balanced_score <- with(
    metrics,
    pmin(ifelse(is.finite(peak100_snr), peak100_snr, 0), 100) /
      (log1p(n_peaks_5)^2)
  )

  list(metrics = metrics, peaks = peak_table, snr = snr)
}

x <- FITSio::readFITS(fits_path)
img <- x$imDat
stopifnot(is.matrix(img))

prep <- fill_invalid(img)
img_fill <- prep$image
valid_mask <- prep$mask
ci <- imager::as.cimg(img_fill)

methods <- list(
  raw = function() {
    mask_like(img, valid_mask)
  },
  gaussian_hp_s2 = function() {
    blur <- cimg_to_matrix(imager::isoblur(ci, sigma = 2))
    mask_like(img_fill - blur, valid_mask)
  },
  median_hp_n2 = function() {
    med <- cimg_to_matrix(imager::medianblur(ci, 2))
    mask_like(img_fill - med, valid_mask)
  },
  tophat_sq7 = function() {
    opn <- cimg_to_matrix(opening_square(ci, 7))
    mask_like(img_fill - opn, valid_mask)
  },
  tophat_sq11 = function() {
    opn <- cimg_to_matrix(opening_square(ci, 11))
    mask_like(img_fill - opn, valid_mask)
  },
  starlet_s1_3_k2 = function() {
    out <- starlet_enhance(
      img,
      keep_scales = 1:3,
      denoise_k = c(2, 2, 1.5, 0, 0),
      J = 5,
      mode = "soft"
    )
    mask_like(out, valid_mask)
  },
  median_hp_starlet = function() {
    med <- cimg_to_matrix(imager::medianblur(ci, 2))
    base <- mask_like(img_fill - med, valid_mask)
    out <- starlet_enhance(
      base,
      keep_scales = 1:2,
      denoise_k = c(1.5, 1.25, 0, 0, 0),
      J = 5,
      mode = "soft"
    )
    mask_like(out, valid_mask)
  },
  tophat11_starlet = function() {
    opn <- cimg_to_matrix(opening_square(ci, 11))
    base <- mask_like(img_fill - opn, valid_mask)
    out <- starlet_enhance(
      base,
      keep_scales = 1:2,
      denoise_k = c(1.5, 1.25, 0, 0, 0),
      J = 5,
      mode = "soft"
    )
    mask_like(out, valid_mask)
  }
)

method_labels <- c(
  raw = "Raw",
  gaussian_hp_s2 = "Gaussian HP (sigma=2)",
  median_hp_n2 = "Median HP (n=2)",
  tophat_sq7 = "Top-hat (sq=7)",
  tophat_sq11 = "Top-hat (sq=11)",
  starlet_s1_3_k2 = "Starlet (s1:3, k=2)",
  median_hp_starlet = "Median HP + Starlet",
  tophat11_starlet = "Top-hat 11 + Starlet"
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
    peaks = scored$peaks,
    snr = scored$snr
  )
}

summary_df <- do.call(rbind, lapply(results, `[[`, "metrics"))
summary_df$label <- method_labels[summary_df$method]
summary_df <- summary_df[order(-summary_df$balanced_score, -summary_df$compact_score), , drop = FALSE]
row.names(summary_df) <- NULL

peak_df <- do.call(rbind, lapply(results, `[[`, "peaks"))
peak_df$label <- method_labels[peak_df$method]

write.csv(summary_df, csv_summary_path, row.names = FALSE)
write.csv(peak_df, csv_peaks_path, row.names = FALSE)

panel_df <- do.call(
  rbind,
  lapply(names(results), function(nm) {
    df <- matrix_to_df(results[[nm]]$image, panel = method_labels[[nm]], step = 2L)
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    df$subtitle <- sprintf(
      "balanced=%.2f | peaks>5sigma=%d | sigma=%.2f",
      met$balanced_score,
      met$n_peaks_5,
      met$bg_sigma
    )
    df
  })
)

peak_overlay <- do.call(
  rbind,
  lapply(names(results), function(nm) {
    pk <- utils::head(results[[nm]]$peaks, 40L)
    if (!nrow(pk)) {
      return(NULL)
    }
    pk$panel <- method_labels[[nm]]
    pk
  })
)

panel_df$facet <- paste0(panel_df$panel, "\n", panel_df$subtitle)
if (!is.null(peak_overlay) && nrow(peak_overlay)) {
  peak_overlay$facet <- paste0(
    peak_overlay$panel,
    "\n",
    sprintf(
      "balanced=%.2f | peaks>5sigma=%d | sigma=%.2f",
      summary_df$balanced_score[match(peak_overlay$method, summary_df$method)],
      summary_df$n_peaks_5[match(peak_overlay$method, summary_df$method)],
      summary_df$bg_sigma[match(peak_overlay$method, summary_df$method)]
    )
  )
}

panel_plot <- ggplot(panel_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  facet_wrap(~facet, ncol = 4) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold")
  )

if (!is.null(peak_overlay) && nrow(peak_overlay)) {
  panel_plot <- panel_plot +
    geom_point(
      data = peak_overlay,
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

metric_long <- reshape2::melt(
  summary_df[, c("label", "balanced_score", "compact_score", "n_peaks_5", "peak100_snr", "bg_sigma", "elapsed_sec")],
  id.vars = "label",
  variable.name = "metric",
  value.name = "value"
)

metric_long$metric <- factor(
  metric_long$metric,
  levels = c("balanced_score", "compact_score", "n_peaks_5", "peak100_snr", "bg_sigma", "elapsed_sec"),
  labels = c("Balanced score", "Compact score", "Peaks > 5 sigma", "100th peak S/N", "Background sigma", "Elapsed seconds")
)
metric_long <- metric_long[is.finite(metric_long$value), , drop = FALSE]

metric_plot <- ggplot(metric_long, aes(x = reorder(label, value), y = value, fill = label)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~metric, scales = "free_x", ncol = 1) +
  scale_fill_viridis_d(option = "C") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

combined_metric_plot <- panel_plot / metric_plot + plot_layout(heights = c(2.3, 1.2))

ggsave(png_panel_path, panel_plot, width = 16, height = 10, dpi = 200, bg = "black")
ggsave(png_metric_path, combined_metric_plot, width = 16, height = 16, dpi = 200, bg = "white")

saveRDS(
  list(
    fits_path = fits_path,
    out_dir = out_dir,
    fill_value = prep$fill,
    summary = summary_df,
    peaks = peak_df,
    results = results
  ),
  rds_path
)

cat("Saved panel PNG to:", png_panel_path, "\n")
cat("Saved metrics PNG to:", png_metric_path, "\n")
cat("Saved summary CSV to:", csv_summary_path, "\n")
cat("Saved peaks CSV to:", csv_peaks_path, "\n")
cat("Saved RDS to:", rds_path, "\n")
cat("\nTop methods by balanced_score:\n")
print(summary_df[, c("label", "balanced_score", "compact_score", "n_peaks_5", "peak100_snr", "bg_sigma", "elapsed_sec")], row.names = FALSE)
