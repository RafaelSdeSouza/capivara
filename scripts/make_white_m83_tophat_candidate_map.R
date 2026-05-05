args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/make_white_m83_tophat_candidate_map.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

results_rds <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_denoise_results.rds")
}

png_path <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_tophat_gc_candidate_map.png")
}

csv_path <- if (length(args) >= 3) {
  args[[3]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_tophat_gc_candidates.csv")
}

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

matrix_to_df <- function(mat, panel, step = 2L) {
  m <- mat[seq(1, nrow(mat), by = step), seq(1, ncol(mat), by = step), drop = FALSE]
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

patchify_candidates <- function(source_mat, peaks, radius = 2L) {
  out <- matrix(NA_real_, nrow(source_mat), ncol(source_mat))
  if (!nrow(peaks)) {
    return(out)
  }

  for (i in seq_len(nrow(peaks))) {
    rr <- seq.int(max(1L, peaks$row[i] - radius), min(nrow(source_mat), peaks$row[i] + radius))
    cc <- seq.int(max(1L, peaks$col[i] - radius), min(ncol(source_mat), peaks$col[i] + radius))
    patch <- source_mat[rr, cc, drop = FALSE]
    cur <- out[rr, cc, drop = FALSE]
    cur[is.na(cur)] <- -Inf
    patch2 <- patch
    patch2[!is.finite(patch2)] <- -Inf
    merged <- pmax(cur, patch2)
    merged[!is.finite(merged)] <- NA_real_
    out[rr, cc] <- merged
  }

  out
}

x <- readRDS(results_rds)
t7 <- x$results$tophat_sq7$image
t11 <- x$results$tophat_sq11$image

t7_info <- extract_peaks(t7, "tophat_sq7", threshold = 5)
t11_info <- extract_peaks(t11, "tophat_sq11", threshold = 5)

match11 <- match_radius(t7_info$peaks, t11_info$peaks, radius = 3)
cand <- t7_info$peaks[match11$match, , drop = FALSE]
cand$nearest_tophat11_dist <- match11$dist[match11$match]
cand$tophat11_snr <- t11_info$snr[cbind(cand$row, cand$col)]
cand <- cand[order(-cand$snr), , drop = FALSE]

utils::write.csv(cand, csv_path, row.names = FALSE)

candidate_map <- patchify_candidates(t7, cand, radius = 2L)

top50 <- utils::head(cand, 50L)

df_resid <- matrix_to_df(t7, "Top-hat (sq=7) Residual", step = 2L)
df_overlay <- matrix_to_df(t7, "Top-hat Residual + Consensus Peaks", step = 2L)
df_cand <- matrix_to_df(candidate_map, sprintf("Consensus Candidate Map\nn=%d peaks", nrow(cand)), step = 2L)

panel_df <- rbind(df_resid, df_overlay, df_cand)

peak_overlay <- top50
peak_overlay$panel <- "Top-hat Residual + Consensus Peaks"

plot_img <- ggplot(panel_df, aes(x = Row, y = Col, fill = value_norm)) +
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

if (nrow(peak_overlay)) {
  plot_img <- plot_img +
    geom_point(
      data = peak_overlay,
      aes(x = row, y = col),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.0,
      stroke = 0.2,
      fill = "#7ce0ff",
      color = "white",
      alpha = 0.8
    )
}

ggsave(png_path, plot_img, width = 16, height = 5.8, dpi = 220, bg = "black")

cat("Saved PNG to:", png_path, "\n")
cat("Saved candidate CSV to:", csv_path, "\n")
cat(sprintf("Consensus candidates kept: %d\n", nrow(cand)))
