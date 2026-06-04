#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

repo_dir <- normalizePath(
  Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara"),
  mustWork = TRUE
)
paper_dir <- normalizePath(
  Sys.getenv(
    "CAPIVARA2_PATH_RUN",
    unset = "/Users/rd23aag/Documents/GitHub/spectropath_paper_workspace_2026-05-08/outputs/ifu_gas_path_maps/manga-7443-12703-LOGCUBE/ifu_halpha_nii_complex"
  ),
  mustWork = TRUE
)
out_dir <- Sys.getenv(
  "CAPIVARA2_OUT",
  unset = file.path(repo_dir, "outputs", "capivara2_kinematics_prototype", basename(dirname(paper_dir)), basename(paper_dir))
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

map_file <- Sys.getenv(
  "CAPIVARA2_MAP_VALUES",
  unset = file.path(paper_dir, paste0(basename(paper_dir), "_path_map_values.csv"))
)
profile_file <- Sys.getenv(
  "CAPIVARA2_PROFILE_LONG",
  unset = file.path(paper_dir, paste0(basename(paper_dir), "_profiles_long.csv.gz"))
)
reference_assignments_file <- Sys.getenv(
  "CAPIVARA2_REFERENCE_ASSIGNMENTS",
  unset = file.path(
    paper_dir,
    "path_cluster_velocity_segments_paper_p3u_p3F_p4T_n50_laplace",
    paste0(basename(paper_dir), "_path_cluster_assignments.csv")
  )
)

ncomp <- as.integer(Sys.getenv("CAPIVARA2_NCOMP", unset = "50"))
knn_k <- as.integer(Sys.getenv("CAPIVARA2_KNN", unset = "40"))
spatial_weight <- as.numeric(Sys.getenv("CAPIVARA2_SPATIAL_WEIGHT", unset = "0.15"))
support_mode <- match.arg(
  Sys.getenv("CAPIVARA2_SUPPORT_MODE", unset = "spatial"),
  c("spatial", "strict")
)
feature_impute <- match.arg(
  Sys.getenv("CAPIVARA2_FEATURE_IMPUTE", unset = "neutral"),
  c("neutral", "median")
)
feature_names <- trimws(strsplit(
  Sys.getenv("CAPIVARA2_FEATURES", unset = "p3u,p3F,p4T"),
  ",",
  fixed = TRUE
)[[1]])

pkgload::load_all(repo_dir, quiet = TRUE)

robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  sc <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(sc) || sc <= 0) sc <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(sc) || sc <= 0) sc <- 1
  (x - med) / sc
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  stats::weighted.mean(x[ok], w[ok])
}

baseline_subtract <- function(v, y, n_edge = 2L) {
  edge <- c(seq_len(min(n_edge, length(y))), seq.int(max(1L, length(y) - n_edge + 1L), length(y)))
  y - stats::median(y[edge], na.rm = TRUE)
}

nonparam_velocity <- function(v, y) {
  line <- baseline_subtract(v, y)
  pos <- pmax(line, 0)
  flux <- sum(pos, na.rm = TRUE)
  if (!is.finite(flux) || flux <= 0) {
    return(data.frame(
      np_flux = NA_real_,
      np_peak_v = NA_real_,
      np_centroid = NA_real_,
      np_sigma = NA_real_,
      np_w80 = NA_real_
    ))
  }
  centroid <- sum(v * pos, na.rm = TRUE) / flux
  sigma <- sqrt(sum((v - centroid)^2 * pos, na.rm = TRUE) / flux)
  ord <- order(v)
  cum <- cumsum(pos[ord]) / flux
  qv <- stats::approx(cum, v[ord], xout = c(0.1, 0.9), ties = "ordered", rule = 2)$y
  data.frame(
    np_flux = flux,
    np_peak_v = v[which.max(line)],
    np_centroid = centroid,
    np_sigma = sigma,
    np_w80 = qv[2] - qv[1]
  )
}

safe_nls <- function(expr, data, start, lower, upper) {
  tryCatch(
    suppressWarnings(
      stats::nls(
        expr,
        data = data,
        start = start,
        algorithm = "port",
        lower = lower,
        upper = upper,
        control = stats::nls.control(maxiter = 500, warnOnly = TRUE, minFactor = 1 / 4096)
      )
    ),
    error = function(e) NULL
  )
}

single_gaussian_velocity <- function(v, y) {
  yy <- as.numeric(y)
  line <- baseline_subtract(v, yy)
  base <- stats::median(yy[c(1, 2, length(yy) - 1, length(yy))], na.rm = TRUE)
  amp <- max(yy, na.rm = TRUE) - base
  if (!is.finite(amp) || amp <= 0) amp <- max(line, na.rm = TRUE)
  mu <- v[which.max(line)]
  np <- nonparam_velocity(v, yy)
  sig <- np$np_sigma
  if (!is.finite(sig) || sig < 40) sig <- 120

  fit <- safe_nls(
    y ~ c0 + amp * exp(-0.5 * ((v - mu) / sig)^2),
    data = data.frame(v = v, y = yy),
    start = list(c0 = base, amp = amp, mu = mu, sig = sig),
    lower = c(c0 = -Inf, amp = 0, mu = min(v), sig = 20),
    upper = c(c0 = Inf, amp = Inf, mu = max(v), sig = 900)
  )

  if (is.null(fit)) {
    return(data.frame(
      gauss_mu = NA_real_,
      gauss_sigma = NA_real_,
      gauss_amp = NA_real_,
      gauss_bic = NA_real_
    ))
  }

  co <- stats::coef(fit)
  rss <- sum(stats::resid(fit)^2, na.rm = TRUE)
  data.frame(
    gauss_mu = unname(co[["mu"]]),
    gauss_sigma = abs(unname(co[["sig"]])),
    gauss_amp = unname(co[["amp"]]),
    gauss_bic = length(v) * log(rss / length(v)) + length(co) * log(length(v))
  )
}

resolve_feature_names <- function(features, available_names) {
  aliases <- c(
    p2 = "levy_tf",
    p3u = "curl_t",
    p3F = "curl_f",
    p4F = "jerk_f",
    p4T = "twist_tf"
  )
  labels <- features
  internal <- unname(ifelse(features %in% names(aliases), aliases[features], features))
  missing <- setdiff(internal, available_names)
  if (length(missing)) {
    stop(
      "Missing requested feature columns after alias resolution: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  data.frame(label = labels, internal = internal, stringsAsFactors = FALSE)
}

make_feature_cube <- function(map_values, features) {
  nx <- max(map_values$x, na.rm = TRUE)
  ny <- max(map_values$y, na.rm = TRUE)
  support <- map_values$inside_spatial_mask %in% TRUE
  feature_complete <- support
  for (feature in features$internal) {
    feature_complete <- feature_complete & is.finite(map_values[[feature]])
  }

  feature_valid <- if (support_mode == "strict") feature_complete else support

  if (sum(feature_valid) < ncomp) {
    stop("Too few finite feature pixels for requested Ncomp.", call. = FALSE)
  }

  cube <- array(
    NA_real_,
    dim = c(nx, ny, nrow(features)),
    dimnames = list(NULL, NULL, features$label)
  )
  for (j in seq_len(nrow(features))) {
    raw <- map_values[[features$internal[j]]]
    finite_feature <- feature_valid & is.finite(raw)

    z <- robust_z(raw[finite_feature])
    plane <- matrix(NA_real_, nrow = nx, ncol = ny)

    if (support_mode == "spatial") {
      fill_value <- if (feature_impute == "median") stats::median(z, na.rm = TRUE) else 0
      plane[cbind(map_values$x[feature_valid], map_values$y[feature_valid])] <- fill_value
    }
    plane[cbind(map_values$x[finite_feature], map_values$y[finite_feature])] <- z
    cube[, , j] <- plane
  }

  mask <- matrix(FALSE, nrow = nx, ncol = ny)
  mask[cbind(map_values$x[feature_valid], map_values$y[feature_valid])] <- TRUE

  list(
    cube = cube,
    mask = mask,
    feature_valid = feature_valid,
    feature_complete = feature_complete
  )
}

map_df <- function(mat) {
  out <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  out$value <- as.vector(mat)
  out
}

cluster_palette <- function(n) {
  if (requireNamespace("Polychrome", quietly = TRUE)) {
    return(Polychrome::glasbey.colors(n))
  }
  grDevices::hcl.colors(n, "Dark 3")
}

zero_limit <- function(x, probs = c(0.02, 0.98), fallback = 100) {
  finite <- x[is.finite(x)]
  if (length(finite) < 2L) return(fallback)
  lim <- max(abs(stats::quantile(finite, probs, na.rm = TRUE)), na.rm = TRUE)
  if (is.finite(lim) && lim > 0) lim else fallback
}

map_values <- utils::read.csv(map_file)
profiles <- utils::read.csv(if (grepl("[.]gz$", profile_file)) gzfile(profile_file) else profile_file)
feature_spec <- resolve_feature_names(feature_names, names(map_values))
feature_info <- make_feature_cube(map_values, feature_spec)

message("Running capivara::segment_large() on path-feature cube")
seg <- segment_large(
  input = list(imDat = feature_info$cube),
  Ncomp = ncomp,
  scale_fn = NULL,
  feature_scale = "none",
  mask = feature_info$mask,
  valid_mode = "finite",
  knn_k = knn_k,
  auto_k = TRUE,
  spatial_weight = spatial_weight,
  return_details = TRUE,
  verbose = TRUE
)

cluster_map <- seg$cluster_map
assignments <- map_values
assignments$feature_valid <- feature_info$feature_valid
assignments$feature_complete <- feature_info$feature_complete
assignments$cluster <- cluster_map[cbind(assignments$x, assignments$y)]
assignments$cluster <- ifelse(is.finite(assignments$cluster), as.integer(assignments$cluster), NA_integer_)

profiles_clustered <- merge(
  profiles,
  assignments[, c("x", "y", "cluster", "line_flux")],
  by = c("x", "y"),
  all.x = FALSE
)
profiles_clustered <- profiles_clustered[is.finite(profiles_clustered$cluster), , drop = FALSE]
profiles_clustered$profile_weight <- ifelse(
  is.finite(profiles_clustered$line_flux.y) & profiles_clustered$line_flux.y > 0,
  profiles_clustered$line_flux.y,
  1
)

stacked <- aggregate(
  cbind(profile, profile_weight) ~ cluster + velocity_kms + wavelength,
  data = profiles_clustered,
  FUN = function(x) c(sum = sum(x, na.rm = TRUE), n = length(x))
)

stacked_profiles <- do.call(
  rbind,
  lapply(split(profiles_clustered, list(profiles_clustered$cluster, profiles_clustered$velocity_kms), drop = TRUE), function(df) {
    data.frame(
      cluster = df$cluster[1],
      velocity_kms = df$velocity_kms[1],
      wavelength = df$wavelength[1],
      profile = weighted_mean_safe(df$profile, df$profile_weight),
      n_spaxels_at_velocity = nrow(df)
    )
  })
)
stacked_profiles <- stacked_profiles[order(stacked_profiles$cluster, stacked_profiles$velocity_kms), ]

velocity_rows <- lapply(split(stacked_profiles, stacked_profiles$cluster), function(df) {
  v <- df$velocity_kms
  y <- df$profile
  cbind(
    data.frame(cluster = df$cluster[1]),
    nonparam_velocity(v, y),
    single_gaussian_velocity(v, y)
  )
})
velocity_summary <- do.call(rbind, velocity_rows)

cluster_summary <- aggregate(
  cbind(n_spaxels = feature_valid, total_line_flux = line_flux) ~ cluster,
  data = assignments[is.finite(assignments$cluster), ],
  FUN = function(x) c(n = length(x), sum = sum(x, na.rm = TRUE), median = stats::median(x, na.rm = TRUE))
)
cluster_stats <- aggregate(
  cbind(line_flux, centroid_kms, median_velocity_kms, sigma_kms, skewness, w80_kms) ~ cluster,
  data = assignments[is.finite(assignments$cluster), ],
  FUN = function(x) stats::median(x, na.rm = TRUE)
)
cluster_traditional_weighted <- do.call(
  rbind,
  lapply(split(assignments[is.finite(assignments$cluster), ], assignments$cluster[is.finite(assignments$cluster)]), function(df) {
    data.frame(
      cluster = df$cluster[1],
      trad_centroid_weighted = weighted_mean_safe(df$centroid_kms, df$line_flux),
      trad_median_velocity_weighted = weighted_mean_safe(df$median_velocity_kms, df$line_flux),
      trad_sigma_weighted = weighted_mean_safe(df$sigma_kms, df$line_flux),
      trad_w80_weighted = weighted_mean_safe(df$w80_kms, df$line_flux)
    )
  })
)
cluster_counts <- as.data.frame(table(assignments$cluster, useNA = "no"))
names(cluster_counts) <- c("cluster", "n_spaxels")
cluster_counts$cluster <- as.integer(as.character(cluster_counts$cluster))
velocity_summary <- merge(cluster_counts, velocity_summary, by = "cluster", all.y = TRUE)
velocity_summary <- merge(velocity_summary, cluster_stats, by = "cluster", all.x = TRUE)
velocity_summary <- merge(velocity_summary, cluster_traditional_weighted, by = "cluster", all.x = TRUE)
velocity_summary$np_centroid_median0 <- velocity_summary$np_centroid -
  stats::median(velocity_summary$np_centroid, na.rm = TRUE)
velocity_summary$gauss_mu_median0 <- velocity_summary$gauss_mu -
  stats::median(velocity_summary$gauss_mu, na.rm = TRUE)
velocity_summary$trad_centroid_weighted_median0 <- velocity_summary$trad_centroid_weighted -
  stats::median(velocity_summary$trad_centroid_weighted, na.rm = TRUE)
velocity_summary$trad_median_velocity_weighted_median0 <- velocity_summary$trad_median_velocity_weighted -
  stats::median(velocity_summary$trad_median_velocity_weighted, na.rm = TRUE)

assignments <- merge(
  assignments,
  velocity_summary[, c(
    "cluster",
    "np_centroid", "np_centroid_median0", "np_sigma", "np_w80",
    "gauss_mu", "gauss_mu_median0", "gauss_sigma",
    "trad_centroid_weighted", "trad_centroid_weighted_median0",
    "trad_median_velocity_weighted", "trad_median_velocity_weighted_median0"
  )],
  by = "cluster",
  all.x = TRUE
)

qc <- data.frame(
  map_file = map_file,
  profile_file = profile_file,
  features = paste(feature_spec$label, collapse = ","),
  internal_feature_columns = paste(feature_spec$internal, collapse = ","),
  Ncomp_requested = ncomp,
  Ncomp_returned = length(unique(stats::na.omit(assignments$cluster))),
  knn_k_used = seg$backend_info$knn_k,
  spatial_weight = spatial_weight,
  support_mode = support_mode,
  feature_impute = feature_impute,
  n_feature_pixels = sum(feature_info$feature_valid),
  n_feature_complete_pixels = sum(feature_info$feature_complete),
  reference_file = reference_assignments_file,
  reference_velocity_column = NA_character_,
  reference_np_centroid_cor = NA_real_,
  reference_gauss_mu_cor = NA_real_,
  reference_trad_centroid_cor = NA_real_
)

if (file.exists(reference_assignments_file)) {
  ref <- utils::read.csv(reference_assignments_file)
  ref_col <- intersect(c("tied_hanii_mu_median0", "single_mu_median0", "np_centroid_median0"), names(ref))
  if (length(ref_col)) {
    cmp <- merge(
      assignments[, c("x", "y", "np_centroid_median0", "gauss_mu_median0", "trad_centroid_weighted_median0")],
      ref[, c("x", "y", ref_col[1])],
      by = c("x", "y")
    )
    ref_y <- cmp[[ref_col[1]]]
    cor_safe <- function(x, y) {
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) <= 5L) return(NA_real_)
      stats::cor(x[ok], y[ok])
    }
    qc$reference_velocity_column <- ref_col[1]
    qc$reference_np_centroid_cor <- cor_safe(cmp$np_centroid_median0, ref_y)
    qc$reference_gauss_mu_cor <- cor_safe(cmp$gauss_mu_median0, ref_y)
    qc$reference_trad_centroid_cor <- cor_safe(cmp$trad_centroid_weighted_median0, ref_y)
  }
}

utils::write.csv(assignments, file.path(out_dir, "capivara2_path_feature_assignments.csv"), row.names = FALSE)
utils::write.csv(stacked_profiles, file.path(out_dir, "capivara2_stacked_line_profiles.csv"), row.names = FALSE)
utils::write.csv(velocity_summary, file.path(out_dir, "capivara2_cluster_velocity_summary.csv"), row.names = FALSE)
utils::write.csv(qc, file.path(out_dir, "capivara2_kinematics_qc.csv"), row.names = FALSE)
saveRDS(
  list(
    segmentation = seg,
    assignments = assignments,
    stacked_profiles = stacked_profiles,
    velocity_summary = velocity_summary,
    qc = qc
  ),
  file.path(out_dir, "capivara2_kinematics_prototype.rds")
)

cluster_df <- map_df(cluster_map)
cluster_df$value <- factor(cluster_df$value)
vel_map <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
vel_map[cbind(assignments$x, assignments$y)] <- assignments$np_centroid_median0
vel_df <- map_df(vel_map)
gauss_vel_map <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
gauss_vel_map[cbind(assignments$x, assignments$y)] <- assignments$gauss_mu_median0
gauss_vel_df <- map_df(gauss_vel_map)
trad_vel_map <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
trad_vel_map[cbind(assignments$x, assignments$y)] <- assignments$trad_centroid_weighted_median0
trad_vel_df <- map_df(trad_vel_map)
flux_map <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
flux_map[cbind(assignments$x, assignments$y)] <- assignments$line_flux
flux_df <- map_df(flux_map)

theme_map <- function() {
  theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

p_cluster <- ggplot(cluster_df, aes(x, y, fill = value)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  scale_fill_manual(values = cluster_palette(max(ncomp, 3L)), na.value = "white", guide = "none") +
  labs(
    title = "Capivara path-feature bins",
    subtitle = paste(feature_spec$label, collapse = ", ")
  ) +
  theme_map()

velocity_palette <- c("#173A68", "#2F6F9F", "#9FC7C4", "#F8F1D0", "#F2B84F", "#D5742F", "#8E2418")

plot_velocity_map <- function(df, title, subtitle) {
  lim <- zero_limit(df$value)
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours = velocity_palette,
      limits = c(-lim, lim),
      oob = scales::squish,
      na.value = "white",
      name = expression(Delta*v)
    ) +
    labs(title = title, subtitle = subtitle) +
    theme_map()
}

p_vel <- plot_velocity_map(
  vel_df,
  "Stacked-profile velocity",
  "non-parametric centroid, median-subtracted"
)

p_gauss <- plot_velocity_map(
  gauss_vel_df,
  "Traditional Gaussian velocity",
  "single Gaussian on stacked profile"
)

p_trad <- plot_velocity_map(
  trad_vel_df,
  "Traditional moment velocity",
  "flux-weighted median of per-spaxel centroids"
)

p_flux <- ggplot(flux_df, aes(x, y, fill = asinh(value))) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  scale_fill_gradientn(
    colours = c("#10162F", "#1E3F6F", "#246B87", "#2F9A9A", "#8DBA72", "#E4C450", "#E8892F", "#9B2F2F"),
    na.value = "white",
    name = "asinh flux"
  ) +
  labs(title = "Line support", subtitle = basename(paper_dir)) +
  theme_map()

panel <- (p_flux | p_cluster | p_vel) / (p_trad | p_gauss | patchwork::plot_spacer()) +
  plot_annotation(
    title = "Capivara 2.0 kinematics prototype",
    subtitle = sprintf(
      "Line-feature cube (%s) -> segment_large(N=%d, k=%d); ref cor: np=%s, gauss=%s, trad=%s",
      paste(feature_spec$label, collapse = ","),
      ncomp,
      qc$knn_k_used,
      ifelse(is.finite(qc$reference_np_centroid_cor), sprintf("%.3f", qc$reference_np_centroid_cor), "NA"),
      ifelse(is.finite(qc$reference_gauss_mu_cor), sprintf("%.3f", qc$reference_gauss_mu_cor), "NA"),
      ifelse(is.finite(qc$reference_trad_centroid_cor), sprintf("%.3f", qc$reference_trad_centroid_cor), "NA")
    ),
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  )

png_path <- file.path(out_dir, "capivara2_kinematics_prototype_panel.png")
pdf_path <- file.path(out_dir, "capivara2_kinematics_prototype_panel.pdf")
ggsave(png_path, panel, width = 13.5, height = 8.8, dpi = 320, bg = "white")
ggsave(pdf_path, panel, width = 13.5, height = 8.8, bg = "white")

message("Wrote:")
message(file.path(out_dir, "capivara2_path_feature_assignments.csv"))
message(file.path(out_dir, "capivara2_cluster_velocity_summary.csv"))
message(file.path(out_dir, "capivara2_stacked_line_profiles.csv"))
message(file.path(out_dir, "capivara2_kinematics_qc.csv"))
message(png_path)
