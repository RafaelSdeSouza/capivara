args <- commandArgs(trailingOnly = TRUE)

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/make_white_m83_best3_figure.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

results_rds <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_denoise_results.rds")
}

png_path <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(repo_root, "outputs", "white_m83_benchmark", "white_m83_best3_methods.png")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(patchwork)
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

x <- readRDS(results_rds)
summary_df <- x$summary
peak_df <- x$peaks
results <- x$results

keep <- c("raw", "tophat_sq7", "tophat_sq11")
labels <- c(
  raw = "Raw",
  tophat_sq7 = "Top-hat (sq=7)",
  tophat_sq11 = "Top-hat (sq=11)"
)

panel_df <- do.call(
  rbind,
  lapply(keep, function(nm) {
    df <- matrix_to_df(results[[nm]]$image, panel = labels[[nm]], step = 2L)
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    df$facet <- sprintf(
      "%s\nbalanced=%.2f | peaks>5sigma=%d",
      labels[[nm]],
      met$balanced_score,
      met$n_peaks_5
    )
    df
  })
)

peak_overlay <- do.call(
  rbind,
  lapply(keep, function(nm) {
    sub <- subset(peak_df, method == nm)
    sub <- utils::head(sub, 40L)
    if (!nrow(sub)) {
      return(NULL)
    }
    met <- summary_df[match(nm, summary_df$method), , drop = FALSE]
    sub$facet <- sprintf(
      "%s\nbalanced=%.2f | peaks>5sigma=%d",
      labels[[nm]],
      met$balanced_score,
      met$n_peaks_5
    )
    sub
  })
)

plot_img <- ggplot(panel_df, aes(x = Row, y = Col, fill = value_norm)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), na.value = "black") +
  facet_wrap(~facet, ncol = 3) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 11, face = "bold")
  )

if (!is.null(peak_overlay) && nrow(peak_overlay)) {
  plot_img <- plot_img +
    geom_point(
      data = peak_overlay,
      aes(x = row, y = col),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.1,
      stroke = 0.2,
      fill = "#7ce0ff",
      color = "white",
      alpha = 0.75
    )
}

ggsave(png_path, plot_img, width = 15, height = 5.8, dpi = 220, bg = "black")

cat("Saved best-3 PNG to:", png_path, "\n")
