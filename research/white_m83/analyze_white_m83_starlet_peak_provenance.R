args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/analyze_white_m83_starlet_peak_provenance.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

results_rds <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_denoise_results.rds")
}

out_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

png_full_path <- file.path(out_dir, "white_m83_starlet_peak_provenance_full.png")
png_zoom_path <- file.path(out_dir, "white_m83_starlet_peak_provenance_zoom.png")
csv_summary_path <- file.path(out_dir, "white_m83_starlet_peak_provenance_summary.csv")
csv_peaks_path <- file.path(out_dir, "white_m83_starlet_peak_provenance_peaks.csv")

suppressPackageStartupMessages({
  library(ggplot2)
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
    return(list(match = rep(FALSE, nrow(query_df)), dist = rep(Inf, nrow(query_df)), idx = rep(NA_integer_, nrow(query_df))))
  }

  qmat <- as.matrix(query_df[, c("row", "col")])
  rmat <- as.matrix(ref_df[, c("row", "col")])

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = rmat, query = qmat, k = 1)
    dist <- nn$nn.dists[, 1]
    idx <- nn$nn.idx[, 1]
  } else {
    dist <- rep(Inf, nrow(qmat))
    idx <- rep(NA_integer_, nrow(qmat))
    for (i in seq_len(nrow(qmat))) {
      d2 <- (rmat[, 1] - qmat[i, 1])^2 + (rmat[, 2] - qmat[i, 2])^2
      j <- which.min(d2)
      dist[i] <- sqrt(d2[j])
      idx[i] <- j
    }
  }

  list(match = is.finite(dist) & dist <= radius, dist = dist, idx = idx)
}

x <- readRDS(results_rds)
results <- x$results

simple_methods <- c("raw", "gaussian_hp_s2", "median_hp_n2", "tophat_sq7", "tophat_sq11")
combo_methods <- c("median_hp_starlet", "tophat11_starlet")
method_labels <- c(
  raw = "Raw",
  gaussian_hp_s2 = "Gaussian HP",
  median_hp_n2 = "Median HP",
  tophat_sq7 = "Top-hat (sq=7)",
  tophat_sq11 = "Top-hat (sq=11)",
  median_hp_starlet = "Median HP + Starlet",
  tophat11_starlet = "Top-hat 11 + Starlet"
)

peak_threshold <- 5
match_radius_px <- 3

peak_info <- lapply(c(simple_methods, combo_methods), function(nm) {
  extract_peaks(results[[nm]]$image, nm, threshold = peak_threshold)
})
names(peak_info) <- c(simple_methods, combo_methods)

raw_image <- results$raw$image
raw_snr <- peak_info$raw$snr
simple_union <- do.call(rbind, lapply(simple_methods, function(nm) peak_info[[nm]]$peaks))
simple_nonraw_union <- do.call(rbind, lapply(setdiff(simple_methods, "raw"), function(nm) peak_info[[nm]]$peaks))

peak_rows <- list()
summary_rows <- list()

for (nm in combo_methods) {
  peaks <- peak_info[[nm]]$peaks
  if (!nrow(peaks)) {
    next
  }

  match_raw <- match_radius(peaks, peak_info$raw$peaks, radius = match_radius_px)
  match_simple <- match_radius(peaks, simple_nonraw_union, radius = match_radius_px)
  match_any <- match_radius(peaks, simple_union, radius = match_radius_px)

  peaks$match_raw <- match_raw$match
  peaks$match_simple <- match_simple$match
  peaks$match_any <- match_any$match
  peaks$nearest_simple_dist <- match_any$dist
  peaks$support_class <- ifelse(
    peaks$match_raw,
    "Seen in Raw",
    ifelse(peaks$match_simple, "Seen in Simple Filter", "Combo-only")
  )

  peaks$raw_snr <- raw_snr[cbind(peaks$row, peaks$col)]
  peaks$gaussian_hp_snr <- peak_info$gaussian_hp_s2$snr[cbind(peaks$row, peaks$col)]
  peaks$median_hp_snr <- peak_info$median_hp_n2$snr[cbind(peaks$row, peaks$col)]
  peaks$tophat7_snr <- peak_info$tophat_sq7$snr[cbind(peaks$row, peaks$col)]
  peaks$tophat11_snr <- peak_info$tophat_sq11$snr[cbind(peaks$row, peaks$col)]
  peaks$max_simple_snr <- do.call(
    pmax,
    c(peaks[, c("raw_snr", "gaussian_hp_snr", "median_hp_snr", "tophat7_snr", "tophat11_snr")], list(na.rm = TRUE))
  )

  combo_only <- peaks[peaks$support_class == "Combo-only", , drop = FALSE]
  summary_rows[[nm]] <- data.frame(
    method = nm,
    label = method_labels[[nm]],
    total_peaks_5 = nrow(peaks),
    seen_in_raw = sum(peaks$support_class == "Seen in Raw"),
    seen_in_simple_filter = sum(peaks$support_class == "Seen in Simple Filter"),
    combo_only = sum(peaks$support_class == "Combo-only"),
    combo_only_frac = mean(peaks$support_class == "Combo-only"),
    combo_only_median_raw_snr = quant_safely(combo_only$raw_snr, 0.5, default = NA_real_),
    combo_only_q90_raw_snr = quant_safely(combo_only$raw_snr, 0.9, default = NA_real_),
    combo_only_median_max_simple_snr = quant_safely(combo_only$max_simple_snr, 0.5, default = NA_real_),
    combo_only_q90_max_simple_snr = quant_safely(combo_only$max_simple_snr, 0.9, default = NA_real_),
    match_radius_px = match_radius_px,
    stringsAsFactors = FALSE
  )

  peak_rows[[nm]] <- peaks
}

peak_df <- do.call(rbind, peak_rows)
summary_df <- do.call(rbind, summary_rows)

write.csv(summary_df, csv_summary_path, row.names = FALSE)
write.csv(peak_df, csv_peaks_path, row.names = FALSE)

color_map <- c(
  "Seen in Raw" = "#7ce0ff",
  "Seen in Simple Filter" = "#ffe066",
  "Combo-only" = "#ff4d6d"
)

full_bg_df <- do.call(
  rbind,
  lapply(combo_methods, function(nm) {
    df <- matrix_to_df(results[[nm]]$image, panel = method_labels[[nm]], step = 2L)
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    df$facet <- sprintf(
      "%s\nRaw=%d | Simple=%d | Combo-only=%d",
      method_labels[[nm]],
      met$seen_in_raw,
      met$seen_in_simple_filter,
      met$combo_only
    )
    df
  })
)

full_peak_df <- do.call(
  rbind,
  lapply(combo_methods, function(nm) {
    pk <- peak_df[peak_df$method == nm, , drop = FALSE]
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    pk$facet <- sprintf(
      "%s\nRaw=%d | Simple=%d | Combo-only=%d",
      method_labels[[nm]],
      met$seen_in_raw,
      met$seen_in_simple_filter,
      met$combo_only
    )
    pk
  })
)

full_plot <- ggplot(full_bg_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  geom_point(
    data = full_peak_df,
    aes(x = row, y = col, color = support_class),
    inherit.aes = FALSE,
    size = 0.28,
    alpha = 0.55
  ) +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  scale_color_manual(values = color_map, breaks = names(color_map)) +
  facet_wrap(~facet, ncol = 2) +
  theme_void() +
  theme(
    text = element_text(color = "white"),
    legend.position = "bottom",
    legend.text = element_text(color = "white"),
    legend.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    strip.background = element_rect(fill = "black", color = NA),
    strip.text = element_text(size = 10, face = "bold", color = "white")
  ) +
  guides(
    fill = "none",
    color = guide_legend(title = NULL, override.aes = list(size = 3, alpha = 1))
  )

ggsave(png_full_path, full_plot, width = 14, height = 7.5, dpi = 220, bg = "black")

max_idx <- which(raw_image == max(raw_image, na.rm = TRUE), arr.ind = TRUE)[1, ]
zoom_half <- 110
zoom_x <- c(max(1, max_idx[1] - zoom_half), min(nrow(raw_image), max_idx[1] + zoom_half))
zoom_y <- c(max(1, max_idx[2] - zoom_half), min(ncol(raw_image), max_idx[2] + zoom_half))

zoom_raw_df <- do.call(
  rbind,
  lapply(combo_methods, function(nm) {
    df <- matrix_to_df(raw_image, panel = method_labels[[nm]], step = 1L)
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    df$facet <- sprintf(
      "%s on Raw Zoom\nCombo-only=%d | median raw S/N of combo-only=%.2f",
      method_labels[[nm]],
      met$combo_only,
      met$combo_only_median_raw_snr
    )
    df
  })
)

zoom_peak_df <- do.call(
  rbind,
  lapply(combo_methods, function(nm) {
    pk <- peak_df[peak_df$method == nm, , drop = FALSE]
    pk <- pk[
      pk$row >= zoom_x[1] & pk$row <= zoom_x[2] &
        pk$col >= zoom_y[1] & pk$col <= zoom_y[2],
      ,
      drop = FALSE
    ]
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    pk$facet <- sprintf(
      "%s on Raw Zoom\nCombo-only=%d | median raw S/N of combo-only=%.2f",
      method_labels[[nm]],
      met$combo_only,
      met$combo_only_median_raw_snr
    )
    pk
  })
)

zoom_plot <- ggplot(zoom_raw_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  geom_point(
    data = zoom_peak_df,
    aes(x = row, y = col, color = support_class),
    inherit.aes = FALSE,
    size = 0.8,
    alpha = 0.8
  ) +
  coord_fixed(xlim = zoom_x, ylim = zoom_y) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  scale_color_manual(values = color_map, breaks = names(color_map)) +
  facet_wrap(~facet, ncol = 2) +
  theme_void() +
  theme(
    text = element_text(color = "white"),
    legend.position = "bottom",
    legend.text = element_text(color = "white"),
    legend.background = element_rect(fill = "black", color = NA),
    plot.background = element_rect(fill = "black", color = NA),
    strip.background = element_rect(fill = "black", color = NA),
    strip.text = element_text(size = 10, face = "bold", color = "white")
  ) +
  guides(
    fill = "none",
    color = guide_legend(title = NULL, override.aes = list(size = 3, alpha = 1))
  )

ggsave(png_zoom_path, zoom_plot, width = 14, height = 7.5, dpi = 220, bg = "black")

cat("Saved full provenance PNG to:", png_full_path, "\n")
cat("Saved zoom provenance PNG to:", png_zoom_path, "\n")
cat("Saved summary CSV to:", csv_summary_path, "\n")
cat("Saved peak CSV to:", csv_peaks_path, "\n")
cat("\nSummary:\n")
print(summary_df, row.names = FALSE)
