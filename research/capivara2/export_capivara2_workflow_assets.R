#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(ragg)
})

repo_dir <- normalizePath(
  Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara"),
  mustWork = TRUE
)

paper_dir <- normalizePath(
  Sys.getenv("CAPIVARA2_PAPER_DIR", unset = "/Users/rd23aag/Documents/GitHub/Capivara2_0_paper"),
  mustWork = TRUE
)

manga_dir <- file.path(paper_dir, "IFU_data", "MaNGA")
cube_path <- Sys.getenv(
  "CAPIVARA2_WORKFLOW_CUBE",
  unset = file.path(manga_dir, "manga-7443-12703-LOGCUBE.fits")
)

out_dir <- Sys.getenv(
  "CAPIVARA2_WORKFLOW_OUT",
  unset = file.path(manga_dir, "Workflow_capivara2_assets")
)

seg_rds <- Sys.getenv(
  "CAPIVARA2_WORKFLOW_SEG_RDS",
  unset = file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "manga_7443_12703_capivara_n25.rds")
)

kin_root <- file.path(
  repo_dir, "outputs", "capivara2_kinematics_prototype",
  "manga-7443-12703-LOGCUBE", "ifu_halpha_nii_complex"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "spectral_channels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "spectral_channels_masked"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "spectra"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "clean_spectra"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "line_profiles"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "ideal_line_profiles"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "path_response_manifold"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "path_features"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "kinematic_heatmaps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "ppxf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "wavelength_slices_original"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "wavelength_slices_masked"), recursive = TRUE, showWarnings = FALSE)

old_internal_path_names <- file.path(
  out_dir, "path_features",
  paste0("capivara2_path_feature_", c("curl_t", "curl_f", "twist_tf", "jerk_f", "levy_tf"), ".png")
)
unlink(old_internal_path_names)
unlink(Sys.glob(file.path(out_dir, "kinematic_heatmaps", "*.png")))
unlink(Sys.glob(file.path(out_dir, "ppxf", "*.png")))
unlink(Sys.glob(file.path(out_dir, "ideal_line_profiles", "ideal_halpha_*.png")))
unlink(Sys.glob(file.path(out_dir, "path_response_manifold", "*.png")))

manifest <- list()

add_manifest <- function(role, path, source = NA_character_, note = NA_character_) {
  manifest[[length(manifest) + 1L]] <<- data.frame(
    role = role,
    file = normalizePath(path, mustWork = FALSE),
    source = source,
    note = note,
    stringsAsFactors = FALSE
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

vangogh <- c(
  "#80B7FF", "#547FFF", "#405CFF", "#263C8B",
  "#FFFAA3", "#FFDE38", "#BFA524"
)

vangogh_discrete <- function(n) {
  if (n <= length(vangogh)) return(vangogh[seq_len(n)])
  grDevices::colorRampPalette(vangogh, space = "Lab")(n)
}

seq_palette <- grDevices::colorRampPalette(
  c("#0B132B", "#1C4E80", "#58A4B0", "#F2D06B", "#F2A541", "#B84A28"),
  space = "Lab"
)(256)

div_palette <- grDevices::colorRampPalette(
  c("#16345C", "#2C7FB8", "#D7E8E6", "#F2D06B", "#D46A2C", "#7F1D1D"),
  space = "Lab"
)(256)

velocity_palette <- grDevices::colorRampPalette(
  c("#08306B", "#2171B5", "#DCEEF2", "#F7F7F7", "#FDD0A2", "#EF6548", "#99000D"),
  space = "Lab"
)(256)

kin_heat_palette <- grDevices::colorRampPalette(
  c("#14213D", "#2F5D7C", "#6FB1A9", "#F5E6A1", "#F2A541", "#A43E2A"),
  space = "Lab"
)(256)

kin_heat_alt_palette <- grDevices::colorRampPalette(
  c("#263C8B", "#547FFF", "#80B7FF", "#FFFAA3", "#FFDE38", "#BFA524"),
  space = "Lab"
)(256)

theme_void_transparent <- function() {
  theme_void() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      plot.margin = grid::unit(c(0, 0, 0, 0), "pt")
    )
}

theme_tiny_spectrum <- function() {
  theme_void(base_size = 7) +
    theme(
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      plot.margin = grid::unit(c(1, 1, 1, 1), "pt")
    )
}

theme_clean_workflow <- function() {
  theme_void(base_size = 8) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = grid::unit(c(5, 5, 5, 5), "pt")
    )
}

map_df <- function(mat) {
  data.frame(
    x = rep(seq_len(nrow(mat)), times = ncol(mat)),
    y = rep(seq_len(ncol(mat)), each = nrow(mat)),
    value = as.vector(mat)
  )
}

robust01 <- function(x, probs = c(0.01, 0.995), positive = FALSE, gamma = 0.9) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  if (positive) x <- pmax(x, 0)
  vals <- x[is.finite(x)]
  if (!length(vals)) return(rep(0, length(x)))
  q <- stats::quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[2] <= q[1]) q <- range(vals, na.rm = TRUE)
  y <- pmin(pmax((x - q[1]) / (q[2] - q[1]), 0), 1)
  y[!is.finite(y)] <- 0
  y^gamma
}

asinh01 <- function(x, probs = c(0.01, 0.995), gamma = 0.85) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x <- pmax(x, 0)
  vals <- x[is.finite(x)]
  if (!length(vals)) return(rep(0, length(x)))
  q <- stats::quantile(vals, probs = probs, na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[2] <= q[1]) q <- range(vals, na.rm = TRUE)
  y <- asinh(pmin(pmax(x, q[1]), q[2]))
  robust01(y, probs = c(0, 1), gamma = gamma)
}

save_numeric_map <- function(mat, out_path, palette = seq_palette, limits = NULL,
                             transparent = TRUE, width = 1600, height = 1600) {
  df <- map_df(mat)
  if (is.null(limits)) {
    vals <- df$value[is.finite(df$value)]
    limits <- if (length(vals)) as.numeric(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE)) else c(0, 1)
  }
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_gradientn(
      colours = palette,
      limits = limits,
      oob = scales::squish,
      na.value = if (transparent) NA else "white"
    ) +
    theme_void_transparent()
  ragg::agg_png(out_path, width = width, height = height, res = 300, background = if (transparent) "transparent" else "white")
  print(p)
  grDevices::dev.off()
}

save_binary_mask <- function(mask, out_path, fill = "#0B132B", transparent = TRUE) {
  mat <- matrix(ifelse(mask, 1, NA_real_), nrow = nrow(mask), ncol = ncol(mask))
  df <- map_df(mat)
  p <- ggplot(df, aes(x, y, fill = factor(value))) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(values = fill, na.value = if (transparent) NA else "white", guide = "none") +
    theme_void_transparent()
  ragg::agg_png(out_path, width = 1600, height = 1600, res = 300, background = if (transparent) "transparent" else "white")
  print(p)
  grDevices::dev.off()
}

save_cluster_map <- function(cluster_map, out_path, transparent = TRUE) {
  df <- map_df(cluster_map)
  df$value <- factor(df$value)
  n <- length(stats::na.omit(unique(df$value)))
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(values = vangogh_discrete(max(n, 1)), guide = "none", na.value = if (transparent) NA else "white") +
    theme_void_transparent()
  ragg::agg_png(out_path, width = 1800, height = 1800, res = 300, background = if (transparent) "transparent" else "white")
  print(p)
  grDevices::dev.off()
}

save_factor_map <- function(mat, out_path, palette = NULL, transparent = TRUE) {
  df <- map_df(mat)
  df$value <- factor(df$value)
  lev <- levels(stats::na.omit(df$value))
  if (is.null(palette)) {
    palette <- vangogh_discrete(max(length(lev), 1))
    names(palette) <- lev
  } else {
    missing <- setdiff(lev, names(palette))
    if (length(missing)) {
      extra <- vangogh_discrete(length(missing))
      names(extra) <- missing
      palette <- c(palette, extra)
    }
  }
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(values = palette, guide = "none", na.value = if (transparent) NA else "white") +
    theme_void_transparent()
  ragg::agg_png(out_path, width = 1800, height = 1800, res = 300, background = if (transparent) "transparent" else "white")
  print(p)
  grDevices::dev.off()
}

smooth_vec <- function(y, k = 61L) {
  k <- max(3L, as.integer(k))
  if (k %% 2L == 0L) k <- k + 1L
  if (length(y) < k) return(y)
  stats::filter(y, rep(1 / k, k), sides = 2)
}

downsample_curve <- function(x, y, n = 420L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (!length(x)) return(data.frame(x = numeric(), y = numeric()))
  if (length(x) <= n) return(data.frame(x = x, y = y))
  idx <- unique(round(seq(1, length(x), length.out = n)))
  data.frame(x = x[idx], y = y[idx])
}

normalize_spectrum <- function(y) {
  y <- as.numeric(y)
  y[!is.finite(y)] <- NA_real_
  sm <- smooth_vec(y, 75L)
  if (all(!is.finite(sm))) sm <- y
  q <- stats::quantile(sm, c(0.02, 0.98), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[2] <= q[1]) q <- range(sm, na.rm = TRUE)
  out <- (sm - q[1]) / (q[2] - q[1])
  out <- pmin(pmax(out, 0), 1)
  out[!is.finite(out)] <- NA_real_
  out
}

save_spectrum_png <- function(wave, flux, out_path, colour = "#1C4E80", width = 900, height = 280) {
  ok_range <- wave >= 3700 & wave <= 7800
  wave <- wave[ok_range]
  flux <- flux[ok_range]
  y <- normalize_spectrum(flux)
  df <- downsample_curve(wave, y, 420L)
  p <- ggplot(df, aes(x, y)) +
    geom_line(linewidth = 0.72, colour = colour, lineend = "round") +
    coord_cartesian(ylim = c(-0.02, 1.02), expand = FALSE) +
    theme_tiny_spectrum()
  ragg::agg_png(out_path, width = width, height = height, res = 300, background = "transparent")
  print(p)
  grDevices::dev.off()
}

save_clean_spectrum_png <- function(wave, flux, out_path, colour = "#1C4E80", width = 980, height = 330) {
  ok_range <- wave >= 3850 & wave <= 7600
  wave <- wave[ok_range]
  flux <- flux[ok_range]
  y <- normalize_spectrum(flux)
  y <- smooth_vec(y, 101L)
  df <- downsample_curve(wave, y, 360L)
  df$y <- pmin(pmax(df$y, 0.04), 0.96)

  p <- ggplot(df, aes(x, y)) +
    geom_line(linewidth = 1.45, colour = colour, lineend = "round") +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    theme_clean_workflow()

  ragg::agg_png(out_path, width = width, height = height, res = 300, background = "white")
  print(p)
  grDevices::dev.off()
}

save_ppxf_fit_example_png <- function(fit_csv, out_path, selected_bin = NULL, width = 980, height = 330) {
  if (!file.exists(fit_csv)) return(FALSE)
  fit <- utils::read.csv(fit_csv)
  if (!all(c("bin", "wavelength_rest", "observed", "bestfit") %in% names(fit))) return(FALSE)
  if (is.null(selected_bin)) {
    selected_bin <- as.integer(names(sort(table(fit$bin), decreasing = TRUE))[1])
  }
  df <- fit[fit$bin == selected_bin, , drop = FALSE]
  df <- df[is.finite(df$wavelength_rest) & is.finite(df$observed) & is.finite(df$bestfit), , drop = FALSE]
  if (!nrow(df)) return(FALSE)
  df <- df[df$wavelength_rest >= 4800 & df$wavelength_rest <= 7000, , drop = FALSE]
  keep <- unique(round(seq(1, nrow(df), length.out = min(600L, nrow(df)))))
  df <- df[keep, , drop = FALSE]
  q <- stats::quantile(c(df$observed, df$bestfit), c(0.02, 0.98), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[2] <= q[1]) q <- range(c(df$observed, df$bestfit), na.rm = TRUE)
  df$observed_n <- pmin(pmax((df$observed - q[1]) / (q[2] - q[1]), 0), 1)
  df$bestfit_n <- pmin(pmax((df$bestfit - q[1]) / (q[2] - q[1]), 0), 1)
  p <- ggplot(df, aes(wavelength_rest)) +
    geom_line(aes(y = observed_n), linewidth = 0.95, colour = "#1C4E80", lineend = "round") +
    geom_line(aes(y = bestfit_n), linewidth = 0.9, colour = "#F2A541", lineend = "round") +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    theme_clean_workflow()
  ragg::agg_png(out_path, width = width, height = height, res = 300, background = "white")
  print(p)
  grDevices::dev.off()
  TRUE
}

save_profile_png <- function(vel, profile, out_path, colour = "#B84A28", width = 900, height = 280) {
  y <- normalize_spectrum(profile)
  df <- downsample_curve(vel, y, 300L)
  p <- ggplot(df, aes(x, y)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "black") +
    geom_line(linewidth = 0.8, colour = colour, lineend = "round") +
    coord_cartesian(ylim = c(-0.02, 1.04), xlim = stats::quantile(vel[is.finite(vel)], c(0.02, 0.98), na.rm = TRUE), expand = FALSE) +
    theme_tiny_spectrum()
  ragg::agg_png(out_path, width = width, height = height, res = 300, background = "transparent")
  print(p)
  grDevices::dev.off()
}

save_ideal_halpha_profile <- function(out_path, colour = "#B84A28", width = 980, height = 330) {
  velocity <- seq(-900, 900, length.out = 500)
  sigma <- 150
  profile <- exp(-0.5 * (velocity / sigma)^2)
  profile <- profile + 0.11 * exp(-0.5 * ((velocity + 250) / 250)^2)
  profile <- profile / max(profile)
  df <- data.frame(velocity = velocity, profile = profile)

  p <- ggplot(df, aes(velocity, profile)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "#D7DFED") +
    geom_line(linewidth = 1.65, colour = colour, lineend = "round") +
    coord_cartesian(xlim = c(-750, 750), ylim = c(-0.02, 1.04), expand = FALSE) +
    theme_clean_workflow()

  ragg::agg_png(out_path, width = width, height = height, res = 300, background = "white")
  print(p)
  grDevices::dev.off()
}

save_ideal_halpha_family <- function(out_dir) {
  velocity <- seq(-900, 900, length.out = 500)
  gaussian <- function(mu, sigma, amp = 1) {
    amp * exp(-0.5 * ((velocity - mu) / sigma)^2)
  }
  specs <- list(
    systemic_narrow = list(
      colour = "#B84A28",
      profile = gaussian(0, 105)
    ),
    blue_shifted = list(
      colour = "#2F74B5",
      profile = gaussian(-190, 125)
    ),
    red_shifted = list(
      colour = "#C7771F",
      profile = gaussian(190, 125)
    ),
    broad_outflow = list(
      colour = "#7A3E12",
      profile = gaussian(0, 255) + gaussian(0, 95, 0.72)
    ),
    blue_wing = list(
      colour = "#1C4E80",
      profile = gaussian(0, 105) + gaussian(-285, 230, 0.34)
    ),
    red_wing = list(
      colour = "#B84A28",
      profile = gaussian(0, 105) + gaussian(285, 230, 0.34)
    ),
    split_components = list(
      colour = "#1E8078",
      profile = gaussian(-170, 110, 0.82) + gaussian(170, 110, 0.78)
    ),
    asymmetric_complex = list(
      colour = "#7A5FA8",
      profile = gaussian(-95, 120, 1.0) + gaussian(240, 180, 0.42)
    )
  )
  files <- character()
  for (nm in names(specs)) {
    sp <- specs[[nm]]
    profile <- sp$profile
    profile <- profile / max(profile)
    df <- data.frame(velocity = velocity, profile = profile)
    out_path <- file.path(out_dir, paste0("ideal_halpha_", nm, ".png"))
    p <- ggplot(df, aes(velocity, profile)) +
      geom_hline(yintercept = 0, linewidth = 0.25, colour = "#D7DFED") +
      geom_line(linewidth = 1.65, colour = sp$colour, lineend = "round") +
      coord_cartesian(xlim = c(-750, 750), ylim = c(-0.02, 1.04), expand = FALSE) +
      theme_clean_workflow()
    ragg::agg_png(out_path, width = 980, height = 330, res = 300, background = "white")
    print(p)
    grDevices::dev.off()
    files <- c(files, out_path)
  }
  files
}

profile_moments <- function(v, f) {
  f <- pmax(f, 0)
  s <- sum(f)
  if (!is.finite(s) || s <= 0) {
    return(c(centroid = 0, sigma = 0, skew = 0, complexity = 0))
  }
  centroid <- sum(v * f) / s
  sigma <- sqrt(sum((v - centroid)^2 * f) / s)
  skew <- if (sigma > 0) sum(((v - centroid) / sigma)^3 * f) / s else 0
  sf <- stats::filter(f, rep(1 / 7, 7), sides = 2)
  sf[!is.finite(sf)] <- f[!is.finite(sf)]
  dy <- diff(sf)
  turns <- sum(diff(sign(dy[is.finite(dy)])) < 0)
  complexity <- min(turns / 3, 1)
  c(centroid = centroid, sigma = sigma, skew = skew, complexity = complexity)
}

make_profile_manifold_data <- function(n = 10L) {
  velocity <- seq(-700, 700, length.out = 180)
  gaussian <- function(mu, sigma, amp = 1) amp * exp(-0.5 * ((velocity - mu) / sigma)^2)
  curve_list <- list()
  bar_list <- list()
  summary <- data.frame()
  for (row in seq_len(n)) {
    y <- (row - 1) / (n - 1)
    for (col in seq_len(n)) {
      x <- (col - 1) / (n - 1)
      mu <- -315 + 630 * x
      sigma <- 58 + 168 * y
      broad_amp <- 0.28 * y
      blue_wing <- 0.86 * y * (1 - x)^1.35
      red_wing <- 0.86 * y * x^1.35
      split_amp <- 0.78 * y * exp(-((x - 0.5) / 0.24)^2)
      profile <- (1 - 0.55 * split_amp) * gaussian(mu, sigma)
      profile <- profile + broad_amp * gaussian(mu, sigma * 2.35)
      profile <- profile + blue_wing * gaussian(mu - 320 - 80 * y, sigma * 1.55)
      profile <- profile + red_wing * gaussian(mu + 320 + 80 * y, sigma * 1.55)
      profile <- profile + split_amp * gaussian(mu - 185, sigma * 0.70)
      profile <- profile + split_amp * gaussian(mu + 185, sigma * 0.70)
      profile <- profile / max(profile)

      moments <- profile_moments(velocity, profile)
      p3u <- min(abs(moments[["skew"]]) / 1.05 + 0.70 * y * abs(x - 0.5), 1)
      p3F <- min((moments[["sigma"]] - 58) / 275 + 0.22 * broad_amp, 1)
      p4T <- min(moments[["complexity"]] * 0.42 + split_amp * 1.35 + 0.18 * y * x * (1 - x), 1)
      activations <- c(p3u = p3u, p3F = p3F, p4T = p4T)

      cell_x0 <- col - 1
      cell_y0 <- n - row
      curve_list[[length(curve_list) + 1L]] <- data.frame(
        cell = paste(row, col, sep = "_"),
        x = cell_x0 + 0.08 + 0.84 * (velocity - min(velocity)) / diff(range(velocity)),
        y = cell_y0 + 0.31 + 0.55 * profile
      )
      bar_width <- 0.18
      bar_gap <- 0.045
      for (j in seq_along(activations)) {
        xmin <- cell_x0 + 0.18 + (j - 1) * (bar_width + bar_gap)
        bar_list[[length(bar_list) + 1L]] <- data.frame(
          cell = paste(row, col, sep = "_"),
          feature = names(activations)[j],
          xmin = xmin,
          xmax = xmin + bar_width,
          ymin = cell_y0 + 0.08,
          ymax = cell_y0 + 0.08 + 0.18 * activations[j]
        )
      }
      summary <- rbind(
        summary,
        data.frame(row = row, col = col, x = x, y = y, t(activations), check.names = FALSE)
      )
    }
  }
  list(curves = do.call(rbind, curve_list), bars = do.call(rbind, bar_list), summary = summary)
}

save_profile_manifold_png <- function(out_path, with_bars = TRUE, n = 10L) {
  dat <- make_profile_manifold_data(n)
  bar_cols <- c(p3u = "#2F74B5", p3F = "#F2A541", p4T = "#7A5FA8")
  p <- ggplot() +
    geom_path(data = dat$curves, aes(x, y, group = cell), colour = "#1C4E80", linewidth = 0.56, lineend = "round") +
    coord_fixed(expand = FALSE, xlim = c(0, n), ylim = c(0, n)) +
    theme_clean_workflow() +
    theme(plot.margin = grid::unit(c(5, 5, 5, 5), "pt"))
  if (with_bars) {
    p <- p +
      geom_rect(
        data = dat$bars,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = feature),
        colour = NA,
        alpha = 0.96
      ) +
      scale_fill_manual(values = bar_cols, guide = "none")
  }
  ragg::agg_png(out_path, width = 1800, height = 1800, res = 300, background = "white")
  print(p)
  grDevices::dev.off()
  invisible(dat)
}

save_profile_response_heatmaps <- function(out_dir, n = 10L) {
  dat <- make_profile_manifold_data(n)$summary
  feature_cols <- c(p3u = "#2F74B5", p3F = "#F2A541", p4T = "#7A5FA8")
  files <- character()
  for (feature in names(feature_cols)) {
    df <- data.frame(
      x = dat$col,
      y = n + 1 - dat$row,
      value = dat[[feature]]
    )
    p <- ggplot(df, aes(x, y, fill = value)) +
      geom_tile(colour = scales::alpha("white", 0.35), linewidth = 0.28) +
      coord_fixed(expand = FALSE) +
      scale_fill_gradientn(
        colours = c("#101729", scales::alpha(feature_cols[[feature]], 0.52), feature_cols[[feature]]),
        limits = c(0, 1),
        guide = "none"
      ) +
      theme_void_transparent()
    out_path <- file.path(out_dir, paste0("ideal_halpha_path_response_", feature, ".png"))
    ragg::agg_png(out_path, width = 1200, height = 1200, res = 300, background = "transparent")
    print(p)
    grDevices::dev.off()
    files <- c(files, out_path)
  }
  files
}

save_spectral_slice <- function(cube, wave, target, out_path, mask = NULL, half_width = 8) {
  idx <- which(abs(wave - target) <= half_width)
  if (!length(idx)) idx <- which.min(abs(wave - target))
  img <- apply(cube[, , idx, drop = FALSE], c(1, 2), stats::median, na.rm = TRUE)
  if (!is.null(mask)) img <- ifelse(mask, img, NA_real_)
  z <- matrix(asinh01(img), nrow = nrow(img), ncol = ncol(img))
  save_numeric_map(z, out_path, palette = c("#05070D", "#1C4E80", "#58A4B0", "#F2D06B", "#FFFFFF"), limits = c(0, 1))
}

save_colored_slice <- function(cube, wave, target, out_path, colour, mask = NULL, half_width = 8) {
  idx <- which(abs(wave - target) <= half_width)
  if (!length(idx)) idx <- which.min(abs(wave - target))
  img <- apply(cube[, , idx, drop = FALSE], c(1, 2), stats::median, na.rm = TRUE)
  if (!is.null(mask)) img <- ifelse(mask, img, NA_real_)
  z <- matrix(asinh01(img), nrow = nrow(img), ncol = ncol(img))
  df <- map_df(z)
  df$alpha <- pmin(pmax((df$value - 0.035) / 0.965, 0), 1)^0.72
  df$alpha[!is.finite(df$value)] <- NA_real_
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster(aes(alpha = alpha), na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_gradientn(
      colours = c(scales::alpha(colour, 0.25), scales::alpha(colour, 0.62), colour),
      values = c(0, 0.58, 1),
      limits = c(0, 1),
      guide = "none",
      na.value = NA
    ) +
    scale_alpha_continuous(range = c(0, 1), guide = "none", na.value = 0) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, colour = NA),
      plot.margin = grid::unit(c(0, 0, 0, 0), "pt")
    )
  ragg::agg_png(out_path, width = 1400, height = 1400, res = 300, background = "transparent")
  print(p)
  grDevices::dev.off()
}

save_rgb_cube <- function(cube, wave, out_path, mask = NULL) {
  targets <- c(R = 6828, G = 5600, B = 4300)
  arr <- array(0, dim = c(dim(cube)[1], dim(cube)[2], 3))
  for (i in seq_along(targets)) {
    idx <- which(abs(wave - targets[i]) <= 25)
    if (!length(idx)) idx <- which.min(abs(wave - targets[i]))
    img <- apply(cube[, , idx, drop = FALSE], c(1, 2), stats::median, na.rm = TRUE)
    if (!is.null(mask)) img <- ifelse(mask, img, NA_real_)
    arr[, , i] <- matrix(asinh01(img), nrow = dim(cube)[1], ncol = dim(cube)[2])
  }
  alpha <- if (is.null(mask)) matrix(1, dim(cube)[1], dim(cube)[2]) else ifelse(mask, 1, 0)
  arr[!is.finite(arr)] <- 0
  alpha[!is.finite(alpha)] <- 0
  df <- data.frame(
    x = rep(seq_len(dim(cube)[1]), times = dim(cube)[2]),
    y = rep(seq_len(dim(cube)[2]), each = dim(cube)[1]),
    fill = grDevices::rgb(as.vector(arr[, , 1]), as.vector(arr[, , 2]), as.vector(arr[, , 3]), alpha = as.vector(alpha))
  )
  p <- ggplot(df, aes(x, y, fill = fill)) +
    geom_raster() +
    scale_fill_identity() +
    coord_equal(expand = FALSE) +
    theme_void_transparent()
  ragg::agg_png(out_path, width = 1800, height = 1800, res = 300, background = "transparent")
  print(p)
  grDevices::dev.off()
}

write_stack_tikz <- function(out_path, subdir, files, command = "CapivaraSpectrumStack") {
  lines <- c(
    "% Auto-generated by scripts/export_capivara2_workflow_assets.R",
    sprintf("\\newcommand{\\%s}[3]{%%", command),
    "\\begin{tikzpicture}[frame/.style={rounded corners=1.4pt,inner sep=1pt,outer sep=0pt,line width=0.025pt,draw=gray}]",
    "\\def\\PanelW{3cm}",
    "\\def\\N{#1}",
    "\\def\\ANGLE{#2}",
    "\\def\\STEP{#3}",
    "\\foreach \\i in {1,...,\\N} {",
    "  \\pgfmathsetmacro{\\k}{\\i-1}",
    "  \\pgfmathsetmacro{\\x}{-\\k * \\STEP * cos(\\ANGLE)}",
    "  \\pgfmathsetmacro{\\y}{-\\k * \\STEP * sin(\\ANGLE)}"
  )
  for (i in seq_along(files)) {
    lines <- c(lines, sprintf("  \\ifnum\\i=%d\\def\\imgname{%s/%s}\\fi", i, subdir, basename(files[i])))
  }
  lines <- c(
    lines,
    "  \\node[frame] at (\\x,\\y) {\\includegraphics[width=\\PanelW]{\\imgname}};",
    "}",
    "\\end{tikzpicture}%",
    "}",
    "",
    sprintf("\\%s{%d}{50}{0.5}", command, length(files))
  )
  writeLines(lines, out_path)
}

write_spectra_dendrogram_tikz <- function(hc, region_ids, selected_ids, out_path) {
  n <- length(hc$order)
  max_h <- max(hc$height)
  leaf_y <- numeric(n)
  leaf_y[hc$order] <- rev(seq_len(n))
  nodes <- vector("list", n - 1L)
  segments <- list()
  get_child <- function(ref) {
    if (ref < 0) {
      idx <- -ref
      list(x = 1, y = leaf_y[idx], leaf = idx)
    } else {
      nodes[[ref]]
    }
  }
  for (i in seq_len(n - 1L)) {
    left <- get_child(hc$merge[i, 1])
    right <- get_child(hc$merge[i, 2])
    x <- if (max_h > 0) 1 - hc$height[i] / max_h else 0
    y <- mean(c(left$y, right$y))
    nodes[[i]] <- list(x = x, y = y)
    segments[[length(segments) + 1L]] <- data.frame(x = x, xend = x, y = left$y, yend = right$y)
    segments[[length(segments) + 1L]] <- data.frame(x = x, xend = left$x, y = left$y, yend = left$y)
    segments[[length(segments) + 1L]] <- data.frame(x = x, xend = right$x, y = right$y, yend = right$y)
  }
  seg <- do.call(rbind, segments)
  sx <- 6.8
  sy <- 0.55
  lines <- c(
    "% Auto-generated by scripts/export_capivara2_workflow_assets.R",
    "\\begin{tikzpicture}[x=1cm,y=1cm,line cap=round,line join=round]",
    "\\definecolor{vangoghblue}{HTML}{2F6FB3}",
    "\\def\\SpecWidth{2.35cm}"
  )
  for (i in seq_len(nrow(seg))) {
    lines <- c(
      lines,
      sprintf(
        "\\draw[line width=0.55pt, color=vangoghblue!85] (%.3f,%.3f) -- (%.3f,%.3f);",
        seg$x[i] * sx, seg$y[i] * sy, seg$xend[i] * sx, seg$yend[i] * sy
      )
    )
  }
  leaf_df <- data.frame(leaf = seq_len(n), region = region_ids, y = leaf_y)
  leaf_df <- leaf_df[leaf_df$region %in% selected_ids, , drop = FALSE]
  for (i in seq_len(nrow(leaf_df))) {
    lines <- c(
      lines,
      sprintf(
        "\\node[anchor=west,inner sep=0pt] at (7.35,%.3f) {\\includegraphics[width=\\SpecWidth]{spectra/region_%02d.png}};",
        leaf_df$y[i] * sy, leaf_df$region[i]
      )
    )
  }
  lines <- c(lines, "\\end{tikzpicture}")
  writeLines(lines, out_path)
}

write_clean_spectra_dendrogram_tikz <- function(hc, region_ids, selected_ids, out_path) {
  n <- length(hc$order)
  max_h <- max(hc$height)
  leaf_y <- numeric(n)
  leaf_y[hc$order] <- rev(seq_len(n))
  nodes <- vector("list", n - 1L)
  segments <- list()
  get_child <- function(ref) {
    if (ref < 0) {
      idx <- -ref
      list(x = 1, y = leaf_y[idx], leaf = idx)
    } else {
      nodes[[ref]]
    }
  }
  for (i in seq_len(n - 1L)) {
    left <- get_child(hc$merge[i, 1])
    right <- get_child(hc$merge[i, 2])
    x <- if (max_h > 0) 1 - hc$height[i] / max_h else 0
    y <- mean(c(left$y, right$y))
    nodes[[i]] <- list(x = x, y = y)
    segments[[length(segments) + 1L]] <- data.frame(x = x, xend = x, y = left$y, yend = right$y)
    segments[[length(segments) + 1L]] <- data.frame(x = x, xend = left$x, y = left$y, yend = left$y)
    segments[[length(segments) + 1L]] <- data.frame(x = x, xend = right$x, y = right$y, yend = right$y)
  }
  seg <- do.call(rbind, segments)
  sx <- 6.8
  sy <- 0.55
  lines <- c(
    "% Auto-generated by scripts/export_capivara2_workflow_assets.R",
    "\\begin{tikzpicture}[x=1cm,y=1cm,line cap=round,line join=round]",
    "\\definecolor{vangoghblue}{HTML}{2F6FB3}",
    "\\def\\SpecWidth{2.35cm}"
  )
  for (i in seq_len(nrow(seg))) {
    lines <- c(
      lines,
      sprintf(
        "\\draw[line width=0.55pt, color=vangoghblue!85] (%.3f,%.3f) -- (%.3f,%.3f);",
        seg$x[i] * sx, seg$y[i] * sy, seg$xend[i] * sx, seg$yend[i] * sy
      )
    )
  }
  leaf_df <- data.frame(leaf = seq_len(n), region = region_ids, y = leaf_y)
  leaf_df <- leaf_df[leaf_df$region %in% selected_ids, , drop = FALSE]
  for (i in seq_len(nrow(leaf_df))) {
    lines <- c(
      lines,
      sprintf(
        "\\node[anchor=west,inner sep=0pt] at (7.35,%.3f) {\\includegraphics[width=\\SpecWidth]{clean_spectra/region_%02d.png}};",
        leaf_df$y[i] * sy, leaf_df$region[i]
      )
    )
  }
  lines <- c(lines, "\\end{tikzpicture}")
  writeLines(lines, out_path)
}

write_matrix_tikz <- function(out_path, selected_ids, subdir = "spectra") {
  lines <- c(
    "% Auto-generated by scripts/export_capivara2_workflow_assets.R",
    "\\begin{tikzpicture}[frame/.style={inner sep=1pt, outer sep=0pt}]",
    "\\def\\W{2.5cm}",
    "\\def\\DY{0.82}"
  )
  for (i in seq_along(selected_ids)) {
    lines <- c(
      lines,
      sprintf("\\node[frame] at (0, -%d*\\DY) {\\includegraphics[width=\\W]{%s/region_%02d.png}};", i - 1L, subdir, selected_ids[i])
    )
  }
  lines <- c(lines, "\\end{tikzpicture}")
  writeLines(lines, out_path)
}

robust_scale_column <- function(x) {
  x <- as.numeric(x)
  vals <- x[is.finite(x)]
  if (!length(vals)) return(rep(0.5, length(x)))
  q <- stats::quantile(vals, c(0.05, 0.95), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(q)) || q[2] <= q[1]) q <- range(vals, na.rm = TRUE)
  y <- pmin(pmax((x - q[1]) / (q[2] - q[1]), 0), 1)
  y[!is.finite(y)] <- 0.5
  y
}

save_heatmap_matrix_png <- function(mat, out_path, palette = kin_heat_palette,
                                    transparent = TRUE, width = 1200, height = 1200) {
  df <- expand.grid(x = seq_len(ncol(mat)), y = seq_len(nrow(mat)))
  df$value <- as.vector(mat[nrow(mat):1, , drop = FALSE])
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_tile(colour = scales::alpha("white", 0.28), linewidth = 0.22) +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(colours = palette, limits = c(0, 1), guide = "none") +
    theme_void_transparent()
  ragg::agg_png(out_path, width = width, height = height, res = 300, background = if (transparent) "transparent" else "white")
  print(p)
  grDevices::dev.off()
}

save_ideal_kinematic_heatmap_png <- function(out_path, palette = kin_heat_palette) {
  mat <- matrix(0.12, nrow = 9, ncol = 9)
  mat[1:4, 1:4] <- matrix(
    c(0.85, 0.72, 0.58, 0.38,
      0.78, 0.64, 0.46, 0.30,
      0.52, 0.43, 0.34, 0.26,
      0.38, 0.31, 0.24, 0.20),
    nrow = 4, byrow = TRUE
  )
  mat[5:8, 5:8] <- matrix(
    c(0.22, 0.35, 0.52, 0.70,
      0.31, 0.48, 0.68, 0.84,
      0.42, 0.60, 0.82, 0.93,
      0.58, 0.74, 0.88, 0.98),
    nrow = 4, byrow = TRUE
  )
  mat[1:3, 6:8] <- matrix(c(0.36, 0.28, 0.22, 0.40, 0.30, 0.24, 0.45, 0.34, 0.26), nrow = 3, byrow = TRUE)
  mat[6:8, 1:3] <- matrix(c(0.63, 0.54, 0.45, 0.72, 0.62, 0.48, 0.82, 0.68, 0.52), nrow = 3, byrow = TRUE)
  diag(mat) <- c(0.96, 0.84, 0.74, 0.62, 0.58, 0.72, 0.86, 0.98, 0.80)
  mat <- pmin(pmax(mat, 0), 1)
  save_heatmap_matrix_png(mat, out_path, palette = palette, transparent = TRUE, width = 1250, height = 1250)
}

copy_asset <- function(src, dest, role, note = NA_character_) {
  if (!file.exists(src)) return(FALSE)
  file.copy(src, dest, overwrite = TRUE)
  add_manifest(role, dest, source = normalizePath(src), note = note)
  TRUE
}

message("Reading cube: ", cube_path)
fits <- FITSio::readFITS(cube_path)
cube <- fits$imDat
wave <- FITSio::axVec(3, fits$axDat)

message("Reading segmentation: ", seg_rds)
seg <- readRDS(seg_rds)
cluster_map <- seg$cluster_map
mask <- if (!is.null(seg$starlet_info$mask)) seg$starlet_info$mask else is.finite(cluster_map)
white <- seg$starlet_info$collapsed %||% apply(cube, c(1, 2), sum, na.rm = TRUE)
recon <- seg$starlet_info$reconstruction %||% white

white_png <- file.path(out_dir, "capivara2_white_light.png")
save_numeric_map(matrix(asinh01(white), nrow = nrow(white), ncol = ncol(white)), white_png, limits = c(0, 1))
add_manifest("white_light", white_png, source = cube_path, note = "Sagui-style white-light collapse from manga-7443-12703.")

recon_png <- file.path(out_dir, "capivara2_starlet_reconstruction.png")
save_numeric_map(matrix(asinh01(recon), nrow = nrow(recon), ncol = ncol(recon)), recon_png, limits = c(0, 1))
add_manifest("starlet_reconstruction", recon_png, source = seg_rds, note = "Selected starlet scales used for support.")

mask_png <- file.path(out_dir, "capivara2_support_binary.png")
save_binary_mask(mask, mask_png, fill = "#0B132B", transparent = FALSE)
add_manifest("support_binary", mask_png, source = seg_rds, note = "Binary support mask with white background.")

mask_trans_png <- file.path(out_dir, "capivara2_support_transparent.png")
save_binary_mask(mask, mask_trans_png, fill = "#0B132B", transparent = TRUE)
add_manifest("support_transparent", mask_trans_png, source = seg_rds, note = "Transparent support mask for TikZ overlays.")

rgb_png <- file.path(out_dir, "capivara2_cube_rgb.png")
save_rgb_cube(cube, wave, rgb_png)
add_manifest("cube_rgb", rgb_png, source = cube_path, note = "Pseudo-RGB from observed spectral windows; no rotation or transpose.")

rgb_masked_png <- file.path(out_dir, "capivara2_cube_rgb_masked.png")
save_rgb_cube(cube, wave, rgb_masked_png, mask = mask)
add_manifest("cube_rgb_masked", rgb_masked_png, source = cube_path, note = "Pseudo-RGB restricted to starlet support.")

slice_targets <- c(3900, 4300, 4861 * 1.0403, 5600, 6200, 6564 * 1.0403, 6585 * 1.0403, 7200, 7800)
slice_files <- character(length(slice_targets))
slice_mask_files <- character(length(slice_targets))
for (i in seq_along(slice_targets)) {
  stem <- sprintf("slice_%02d_%04dA.png", i, round(slice_targets[i]))
  slice_files[i] <- file.path(out_dir, "spectral_channels", stem)
  slice_mask_files[i] <- file.path(out_dir, "spectral_channels_masked", stem)
  save_spectral_slice(cube, wave, slice_targets[i], slice_files[i], mask = NULL)
  save_spectral_slice(cube, wave, slice_targets[i], slice_mask_files[i], mask = mask)
  add_manifest("spectral_channel", slice_files[i], source = cube_path, note = sprintf("Median slice near %.1f Angstrom.", slice_targets[i]))
  add_manifest("spectral_channel_masked", slice_mask_files[i], source = cube_path, note = sprintf("Masked median slice near %.1f Angstrom.", slice_targets[i]))
}

wavelength_colours <- grDevices::colorRampPalette(
  c("#1C4E80", "#2F74B5", "#58A4B0", "#1E8078", "#F2D06B", "#F2A541", "#C7771F", "#B84A28"),
  space = "Lab"
)(length(slice_targets))

for (i in seq_along(slice_targets)) {
  stem <- sprintf("wavelength_%02d_%04dA.png", i, round(slice_targets[i]))
  out_original <- file.path(out_dir, "wavelength_slices_original", stem)
  out_masked <- file.path(out_dir, "wavelength_slices_masked", stem)
  save_colored_slice(cube, wave, slice_targets[i], out_original, colour = wavelength_colours[i], mask = NULL)
  save_colored_slice(cube, wave, slice_targets[i], out_masked, colour = wavelength_colours[i], mask = mask)
  add_manifest("wavelength_slice_original", out_original, source = cube_path, note = sprintf("Wavelength-coded original cube slice near %.1f Angstrom.", slice_targets[i]))
  add_manifest("wavelength_slice_masked", out_masked, source = cube_path, note = sprintf("Wavelength-coded masked cube slice near %.1f Angstrom.", slice_targets[i]))
}

halpha_target <- 6564 * 1.0403
halpha_slice_original <- file.path(out_dir, "wavelength_slices_original", "capivara2_halpha_slice_original.png")
halpha_slice_masked <- file.path(out_dir, "wavelength_slices_masked", "capivara2_halpha_slice_masked.png")
save_colored_slice(cube, wave, halpha_target, halpha_slice_original, colour = "#B84A28", mask = NULL)
save_colored_slice(cube, wave, halpha_target, halpha_slice_masked, colour = "#B84A28", mask = mask)
add_manifest("halpha_slice_original", halpha_slice_original, source = cube_path, note = "Standalone wavelength-coded Halpha observed-frame slice, native orientation.")
add_manifest("halpha_slice_masked", halpha_slice_masked, source = cube_path, note = "Standalone wavelength-coded Halpha observed-frame slice after starlet support, native orientation.")

write_stack_tikz(file.path(out_dir, "capivara2_spectral_stack.tikz"), "spectral_channels", slice_files, "CapivaraSpectrumStack")
write_stack_tikz(file.path(out_dir, "capivara2_spectral_stack_masked.tikz"), "spectral_channels_masked", slice_mask_files, "CapivaraMaskedSpectrumStack")
add_manifest("tikz_helper", file.path(out_dir, "capivara2_spectral_stack.tikz"), note = "TikZ stack helper for spectral-channel cube.")
add_manifest("tikz_helper", file.path(out_dir, "capivara2_spectral_stack_masked.tikz"), note = "TikZ stack helper for masked spectral-channel cube.")

seg_png <- file.path(out_dir, "capivara2_segmentation_vangogh.png")
save_cluster_map(cluster_map, seg_png, transparent = TRUE)
add_manifest("segmentation", seg_png, source = seg_rds, note = "Capivara N=25 segmentation, Van Gogh categorical palette.")

flux_mat <- matrix(cube, nrow = dim(cube)[1] * dim(cube)[2], ncol = dim(cube)[3])
cl_vec <- as.vector(cluster_map)
cluster_ids <- sort(unique(cl_vec[is.finite(cl_vec)]))
sum_spectra <- matrix(NA_real_, nrow = length(cluster_ids), ncol = ncol(flux_mat))
rownames(sum_spectra) <- cluster_ids
for (i in seq_along(cluster_ids)) {
  idx <- which(cl_vec == cluster_ids[i])
  sum_spectra[i, ] <- colSums(flux_mat[idx, , drop = FALSE], na.rm = TRUE)
}

for (i in seq_along(cluster_ids)) {
  out <- file.path(out_dir, "spectra", sprintf("region_%02d.png", cluster_ids[i]))
  save_spectrum_png(wave, sum_spectra[i, ], out, colour = vangogh_discrete(length(cluster_ids))[i])
  add_manifest("region_spectrum", out, source = cube_path, note = "Real flux-preserving summed Capivara spectrum, smoothed for tiny schematic display.")
}

spectra_ok <- wave >= 3700 & wave <= 7800
norm_spectra <- t(apply(sum_spectra[, spectra_ok, drop = FALSE], 1, normalize_spectrum))
norm_spectra[!is.finite(norm_spectra)] <- 0
hc <- stats::hclust(stats::dist(norm_spectra), method = "ward.D2")
selected <- cluster_ids[hc$order[unique(round(seq(1, length(hc$order), length.out = min(7, length(hc$order)))))]]

clean_cols <- grDevices::colorRampPalette(c("#1C4E80", "#58A4B0", "#F2A541", "#B84A28"), space = "Lab")(length(cluster_ids))
for (i in seq_along(cluster_ids)) {
  out <- file.path(out_dir, "clean_spectra", sprintf("region_%02d.png", cluster_ids[i]))
  save_clean_spectrum_png(wave, sum_spectra[i, ], out, colour = clean_cols[i])
  add_manifest("clean_region_spectrum", out, source = cube_path, note = "Clean standalone summed spectrum for workflow dendrogram assets.")
}

spectra_matrix <- norm_spectra[match(cluster_ids[hc$order], cluster_ids), , drop = FALSE]
spectra_matrix_png <- file.path(out_dir, "capivara2_spectra_matrix.png")
df_mat <- expand.grid(x = seq_len(ncol(spectra_matrix)), y = seq_len(nrow(spectra_matrix)))
df_mat$value <- as.vector(spectra_matrix[nrow(spectra_matrix):1, ])
p_mat <- ggplot(df_mat, aes(x, y, fill = value)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  scale_fill_gradientn(colours = seq_palette, limits = c(0, 1), guide = "none") +
  theme_void_transparent()
ragg::agg_png(spectra_matrix_png, width = 1800, height = 1200, res = 300, background = "transparent")
print(p_mat)
grDevices::dev.off()
add_manifest("spectra_matrix", spectra_matrix_png, source = cube_path, note = "Rows are Capivara summed spectra, ordered by Ward dendrogram.")

dist_mat <- as.matrix(stats::dist(norm_spectra))
dist_mat <- dist_mat[hc$order, hc$order, drop = FALSE]
dmax <- stats::quantile(dist_mat[upper.tri(dist_mat)], 0.98, na.rm = TRUE)
aff <- 1 - pmin(dist_mat / dmax, 1)
diag(aff) <- 1
dist_png <- file.path(out_dir, "capivara2_spectral_dissimilarity.png")
df_dist <- expand.grid(x = seq_len(ncol(aff)), y = seq_len(nrow(aff)))
df_dist$value <- as.vector(aff[nrow(aff):1, ])
p_dist <- ggplot(df_dist, aes(x, y, fill = value)) +
  geom_tile(colour = scales::alpha("white", 0.25), linewidth = 0.35) +
  coord_fixed(expand = FALSE) +
  scale_fill_gradientn(colours = c("#1E2A53", "#355C7D", "#7AA6A1", "#E3C565", "#F3E6AF"), limits = c(0, 1), guide = "none") +
  theme_void_transparent()
ragg::agg_png(dist_png, width = 1600, height = 1600, res = 300, background = "transparent")
print(p_dist)
grDevices::dev.off()
add_manifest("spectral_dissimilarity", dist_png, source = cube_path, note = "Affinity-style visualization derived from distances between real region spectra.")

write_matrix_tikz(file.path(out_dir, "capivara2_spectra_matrix.tikz"), selected)
write_matrix_tikz(file.path(out_dir, "capivara2_clean_spectra_matrix.tikz"), selected, subdir = "clean_spectra")
write_spectra_dendrogram_tikz(hc, cluster_ids, selected, file.path(out_dir, "capivara2_spectral_dendrogram.tikz"))
write_clean_spectra_dendrogram_tikz(hc, cluster_ids, selected, file.path(out_dir, "capivara2_clean_spectral_dendrogram.tikz"))
add_manifest("tikz_helper", file.path(out_dir, "capivara2_spectra_matrix.tikz"), note = "Tiny real spectra stack for TikZ.")
add_manifest("tikz_helper", file.path(out_dir, "capivara2_clean_spectra_matrix.tikz"), note = "Clean tiny spectra stack for TikZ.")
add_manifest("tikz_helper", file.path(out_dir, "capivara2_spectral_dendrogram.tikz"), note = "Dendrogram helper using tiny real spectra PNGs.")
add_manifest("tikz_helper", file.path(out_dir, "capivara2_clean_spectral_dendrogram.tikz"), note = "Dendrogram helper using clean standalone spectra PNGs.")

ideal_halpha <- file.path(out_dir, "ideal_line_profiles", "ideal_halpha_profile.png")
save_ideal_halpha_profile(ideal_halpha)
add_manifest("ideal_halpha_profile", ideal_halpha, note = "Idealized Halpha line profile for the workflow schematic.")

ideal_family <- save_ideal_halpha_family(file.path(out_dir, "ideal_line_profiles"))
for (ff in ideal_family) {
  add_manifest("ideal_halpha_profile_variant", ff, note = "Idealized Halpha profile variant for workflow schematics.")
}

profile_manifold <- file.path(out_dir, "path_response_manifold", "ideal_halpha_profile_manifold.png")
save_profile_manifold_png(profile_manifold, with_bars = FALSE, n = 10L)
add_manifest("path_response_manifold", profile_manifold, note = "Idealized 10x10 Halpha profile manifold; corners converge to distinct kinematic line shapes.")

profile_activation <- file.path(out_dir, "path_response_manifold", "ideal_halpha_path_activation_manifold.png")
save_profile_manifold_png(profile_activation, with_bars = TRUE, n = 10L)
add_manifest("path_response_manifold", profile_activation, note = "Idealized 10x10 Halpha profile manifold with tiny p3u/p3F/p4T activation bars.")

response_heatmaps <- save_profile_response_heatmaps(file.path(out_dir, "path_response_manifold"), n = 10L)
for (ff in response_heatmaps) {
  add_manifest("path_response_heatmap", ff, note = "Standalone path-feature response heatmap over the idealized profile manifold.")
}

assign_csv <- file.path(kin_root, "capivara2_path_feature_assignments.csv")
kin_rds <- file.path(kin_root, "capivara2_kinematics_prototype.rds")
if (file.exists(assign_csv)) {
  assign <- utils::read.csv(assign_csv)
  make_assignment_map <- function(col) {
    mat <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
    ok <- assign$x >= 1 & assign$x <= nrow(mat) & assign$y >= 1 & assign$y <= ncol(mat)
    mat[cbind(assign$x[ok], assign$y[ok])] <- assign[[col]][ok]
    mat
  }
  feature_specs <- data.frame(
    label = c("p3u", "p3F", "p4T", "line_flux", "np_centroid_median0"),
    column = c("curl_t", "curl_f", "twist_tf", "line_flux", "np_centroid_median0"),
    stringsAsFactors = FALSE
  )
  for (ii in seq_len(nrow(feature_specs))) {
    label <- feature_specs$label[ii]
    col <- feature_specs$column[ii]
    if (!col %in% names(assign)) next
    mat <- make_assignment_map(col)
    vals <- mat[is.finite(mat)]
    lim <- if (length(vals)) {
      if (col %in% c("curl_t", "curl_f", "twist_tf", "np_centroid_median0")) {
        a <- max(abs(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE)))
        c(-a, a)
      } else {
        as.numeric(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE))
      }
    } else c(0, 1)
    out <- file.path(out_dir, "path_features", paste0("capivara2_path_feature_", label, ".png"))
    save_numeric_map(mat, out, palette = if (col %in% c("line_flux")) seq_palette else div_palette, limits = lim)
    add_manifest("path_feature_map", out, source = assign_csv, note = paste("Standalone Halpha+[NII] path feature map:", label, "(internal column:", col, ")"))
  }

  heat_specs <- data.frame(
    label = c("p3u", "p3F", "p4T", "velocity", "sigma", "w80", "asymmetry"),
    column = c("curl_t", "curl_f", "twist_tf", "np_centroid_median0", "np_sigma", "np_w80", "skewness"),
    stringsAsFactors = FALSE
  )
  heat_specs <- heat_specs[heat_specs$column %in% names(assign), , drop = FALSE]
  if (nrow(heat_specs)) {
    keep <- is.finite(assign$cluster)
    by_bin <- stats::aggregate(
      assign[keep, heat_specs$column, drop = FALSE],
      by = list(cluster = assign$cluster[keep]),
      FUN = function(z) stats::median(z, na.rm = TRUE)
    )
    by_bin <- by_bin[order(by_bin$cluster), , drop = FALSE]
    feature_mat <- as.matrix(by_bin[, heat_specs$column, drop = FALSE])
    colnames(feature_mat) <- heat_specs$label
    feature_mat <- apply(feature_mat, 2, robust_scale_column)
    if (is.null(dim(feature_mat))) feature_mat <- matrix(feature_mat, ncol = 1)
    rownames(feature_mat) <- by_bin$cluster

    finite_rows <- apply(feature_mat, 1, function(z) all(is.finite(z)))
    ordered_mat <- feature_mat[finite_rows, , drop = FALSE]
    if (nrow(ordered_mat) > 2L) {
      row_order <- stats::hclust(stats::dist(ordered_mat), method = "ward.D2")$order
      ordered_mat <- ordered_mat[row_order, , drop = FALSE]
    }

    out <- file.path(out_dir, "kinematic_heatmaps", "capivara2_kinematic_feature_heatmap.png")
    save_heatmap_matrix_png(ordered_mat, out, palette = kin_heat_palette, transparent = TRUE, width = 1250, height = 1450)
    add_manifest("kinematic_heatmap", out, source = assign_csv, note = "No-bar workflow heatmap of binned kinematic/path features: p3u, p3F, p4T, velocity, sigma, w80, asymmetry.")

    out_alt <- file.path(out_dir, "kinematic_heatmaps", "capivara2_kinematic_feature_heatmap_sagui_palette.png")
    save_heatmap_matrix_png(ordered_mat, out_alt, palette = kin_heat_alt_palette, transparent = TRUE, width = 1250, height = 1450)
    add_manifest("kinematic_heatmap", out_alt, source = assign_csv, note = "No-bar workflow heatmap of binned kinematic/path features using the Sagui-like blue/yellow palette.")

    if (nrow(ordered_mat) > 2L) {
      d <- as.matrix(stats::dist(ordered_mat))
      dmax <- stats::quantile(d[upper.tri(d)], 0.98, na.rm = TRUE)
      affinity <- 1 - pmin(d / dmax, 1)
      diag(affinity) <- 1
      out <- file.path(out_dir, "kinematic_heatmaps", "capivara2_kinematic_dissimilarity_heatmap.png")
      save_heatmap_matrix_png(affinity, out, palette = kin_heat_palette, transparent = TRUE, width = 1400, height = 1400)
      add_manifest("kinematic_heatmap", out, source = assign_csv, note = "No-bar workflow heatmap of kinematic/path-feature affinity between bins.")
    }
  }

  out <- file.path(out_dir, "kinematic_heatmaps", "capivara2_ideal_kinematic_heatmap.png")
  save_ideal_kinematic_heatmap_png(out, palette = kin_heat_palette)
  add_manifest("kinematic_heatmap_ideal", out, note = "Idealized no-bar kinematic heatmap glyph for workflow schematics.")
}

if (file.exists(kin_rds)) {
  kin <- readRDS(kin_rds)
  if (!is.null(kin$segmentation$cluster_map)) {
    out <- file.path(out_dir, "path_features", "capivara2_halpha_nii_path_bins_vangogh.png")
    save_cluster_map(kin$segmentation$cluster_map, out, transparent = TRUE)
    add_manifest("path_segmentation", out, source = kin_rds, note = "Halpha+[NII] path-feature Capivara segmentation.")
  }
  if (!is.null(kin$stacked_profiles)) {
    prof <- kin$stacked_profiles
    pids <- sort(unique(prof$cluster))
    keep <- pids[unique(round(seq(1, length(pids), length.out = min(7, length(pids)))))]
    for (id in keep) {
      pp <- prof[prof$cluster == id, , drop = FALSE]
      out <- file.path(out_dir, "line_profiles", sprintf("line_profile_region_%02d.png", id))
      save_profile_png(pp$velocity_kms, pp$profile, out)
      add_manifest("line_profile", out, source = kin_rds, note = "Real stacked Halpha+[NII] line profile, smoothed for schematic display.")
    }
  }
}

ppxf_population_rds <- file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "local_ppxf_figures", "manga_7443_12703_local_ppxf_laplace_maps.rds")
ppxf_emission_rds <- file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "maps", "manga_7443_12703_ppxf_emission_laplace_maps.rds")
ppxf_fit_csv <- file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "local_ppxf_figures", "manga_7443_12703_local_ppxf_fit_examples_spectra.csv")

ppxf_bins <- file.path(out_dir, "ppxf", "ppxf_bins.png")
save_cluster_map(cluster_map, ppxf_bins, transparent = TRUE)
add_manifest("ppxf_workflow_map", ppxf_bins, source = seg_rds, note = "No-bar pPXF workflow bin map using the Capivara segmentation.")

if (file.exists(ppxf_population_rds)) {
  pop_maps <- readRDS(ppxf_population_rds)
  lap <- pop_maps$laplace %||% pop_maps$raw
  pop_specs <- list(
    ppxf_stellar_velocity = list(mat = lap$stellar_vel_median0 %||% lap$stellar_vel, palette = div_palette, limits = NULL),
    ppxf_stellar_sigma = list(mat = lap$stellar_sigma, palette = seq_palette, limits = NULL),
    ppxf_log_age = list(mat = lap$mean_log_age, palette = seq_palette, limits = NULL),
    ppxf_metallicity = list(mat = lap$mean_metal, palette = seq_palette, limits = NULL),
    ppxf_ml_r = list(mat = lap$ml_r, palette = seq_palette, limits = NULL)
  )
  for (nm in names(pop_specs)) {
    mat <- pop_specs[[nm]]$mat
    if (is.null(mat)) next
    vals <- mat[is.finite(mat)]
    lim <- if (length(vals)) {
      if (identical(pop_specs[[nm]]$palette, div_palette)) {
        a <- max(abs(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE)))
        c(-a, a)
      } else {
        as.numeric(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE))
      }
    } else c(0, 1)
    out <- file.path(out_dir, "ppxf", paste0(nm, ".png"))
    save_numeric_map(mat, out, palette = pop_specs[[nm]]$palette, limits = lim)
    add_manifest("ppxf_workflow_map", out, source = ppxf_population_rds, note = paste("No-bar pPXF workflow map:", nm))
  }
}

if (file.exists(ppxf_emission_rds)) {
  em_maps <- readRDS(ppxf_emission_rds)
  maps <- em_maps$maps %||% list()
  gas_velocity_centered <- maps$gas_velocity_pattern_laplace
  if (!is.null(gas_velocity_centered)) {
    vals <- gas_velocity_centered[is.finite(gas_velocity_centered)]
    if (length(vals)) gas_velocity_centered <- gas_velocity_centered - stats::median(vals, na.rm = TRUE)
  }
  em_specs <- list(
    ppxf_halpha_log_flux = list(mat = maps$halpha_log_flux_laplace, palette = seq_palette, limits = NULL),
    ppxf_gas_velocity = list(mat = gas_velocity_centered, palette = velocity_palette, limits = NULL),
    ppxf_log_nii_ha = list(mat = maps$log_nii_ha_laplace, palette = div_palette, limits = NULL),
    ppxf_log_oiii_hb = list(mat = maps$log_oiii_hb_laplace, palette = div_palette, limits = NULL)
  )
  for (nm in names(em_specs)) {
    mat <- em_specs[[nm]]$mat
    if (is.null(mat)) next
    vals <- mat[is.finite(mat)]
    lim <- if (length(vals)) {
      if (identical(em_specs[[nm]]$palette, div_palette)) {
        a <- max(abs(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE)))
        c(-a, a)
      } else {
        as.numeric(stats::quantile(vals, c(0.02, 0.98), na.rm = TRUE))
      }
    } else c(0, 1)
    out <- file.path(out_dir, "ppxf", paste0(nm, ".png"))
    save_numeric_map(mat, out, palette = em_specs[[nm]]$palette, limits = lim)
    note <- if (nm == "ppxf_gas_velocity") {
      "No-bar pPXF workflow map: ppxf_gas_velocity, median-centered for schematic diverging colors; not absolute calibrated."
    } else {
      paste("No-bar pPXF workflow map:", nm)
    }
    add_manifest("ppxf_workflow_map", out, source = ppxf_emission_rds, note = note)
  }
  if (!is.null(maps$bpt_class_laplace)) {
    bpt_cols <- c(
      "Star-forming" = "#2F74B5", "Composite" = "#58A4B0", "AGN" = "#B84A28",
      "LINER" = "#7A3E12", "Seyfert" = "#C7771F", "Unclassified" = "#D7DFED"
    )
    out <- file.path(out_dir, "ppxf", "ppxf_bpt_class.png")
    save_factor_map(maps$bpt_class_laplace, out, palette = bpt_cols, transparent = TRUE)
    add_manifest("ppxf_workflow_map", out, source = ppxf_emission_rds, note = "No-bar pPXF BPT class workflow map.")
  }
}

fit_out <- file.path(out_dir, "ppxf", "ppxf_fit_example_bin01.png")
if (save_ppxf_fit_example_png(ppxf_fit_csv, fit_out, selected_bin = 1L)) {
  add_manifest("ppxf_fit_example", fit_out, source = ppxf_fit_csv, note = "No-bar pPXF observed/best-fit spectrum glyph.")
}

readme <- c(
  "Capivara 2.0 workflow PNG assets",
  "",
  "Generated by scripts/export_capivara2_workflow_assets.R.",
  paste("Cube:", cube_path),
  paste("Segmentation:", seg_rds),
  "",
  "Important conventions:",
  "- Images are written in the native Capivara matrix orientation.",
  "- No transpose, image rotation, or vertical flip is applied.",
  "- spectra/ contains real flux-preserving summed spectra from Capivara regions.",
  "- clean_spectra/ contains simplified standalone spectra for workflow dendrograms.",
  "- ideal_line_profiles/ contains idealized Halpha profiles for schematic line-profile nodes.",
  "- path_response_manifold/ contains idealized profile manifolds and p3u/p3F/p4T activation glyphs.",
  "- kinematic_heatmaps/ contains no-bar path/kinematic feature heatmap glyphs.",
  "- wavelength_slices_original/ and wavelength_slices_masked/ contain wavelength-coded slices.",
  "- Tiny spectra are smoothed only for schematic readability.",
  "- ppxf/ contains no-bar standalone fitting maps regenerated from RDS products.",
  "- path_features/ uses spectrograph-style public names such as p3u, p3F, and p4T."
)
writeLines(readme, file.path(out_dir, "README_assets.txt"))
add_manifest("readme", file.path(out_dir, "README_assets.txt"), note = "Asset generation notes.")

manifest_df <- if (length(manifest)) do.call(rbind, manifest) else data.frame()
utils::write.csv(manifest_df, file.path(out_dir, "capivara2_workflow_assets_manifest.csv"), row.names = FALSE)

message("Wrote Capivara 2 workflow assets to:")
message(out_dir)
message("Manifest:")
message(file.path(out_dir, "capivara2_workflow_assets_manifest.csv"))
