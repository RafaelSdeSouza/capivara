args <- commandArgs(trailingOnly = TRUE)
use_env_inputs <- tolower(Sys.getenv("CAPIVARA_USE_ENV_INPUTS", unset = "false")) %in% c("1", "true", "yes", "y", "on")

# Usage:
#   Rscript scripts/run_manga_10218_bar_merger_maps.R [cube_path] [out_dir]
#
# Common knobs:
#   CAPIVARA_REDSHIFT=0.0461
#   CAPIVARA_LINE=halpha            # aliases include hbeta, oiii5007, nii6583, sii6716
#   CAPIVARA_LINE_REST=6562.8       # optional override, Angstrom
#   CAPIVARA_OUTPUT_PREFIX=manga10218
#   CAPIVARA_KNN=100 CAPIVARA_NCOMP=25
#   CAPIVARA_PATH_KNN=100 CAPIVARA_PATH_NCOMP=45 CAPIVARA_PATH_SPATIAL_WEIGHT=0.10
#
# If [out_dir] is omitted, products are written beside the cube in:
#   dirname(cube_path)/capivara_outputs

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts/run_manga_10218_bar_merger_maps.R")
repo_root <- Sys.getenv("CAPIVARA_REPO_ROOT", unset = "")
if (!nzchar(repo_root)) {
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
}

cube_path <- if (!use_env_inputs && length(args) >= 1) {
  args[[1]]
} else {
  env_cube_path <- Sys.getenv("CAPIVARA_CUBE_PATH", unset = "")
  if (nzchar(env_cube_path)) env_cube_path else "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/bar_merger/manga-10218-12703-LOGCUBE.fits"
}
out_dir <- if (!use_env_inputs && length(args) >= 2) {
  args[[2]]
} else {
  env_out_dir <- Sys.getenv("CAPIVARA_OUTPUT_DIR", unset = "")
  if (nzchar(env_out_dir)) env_out_dir else file.path(dirname(cube_path), "capivara_outputs")
}

env_num <- function(keys, unset) {
  for (key in keys) {
    value <- Sys.getenv(key, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      return(as.numeric(value))
    }
  }
  as.numeric(unset)
}

env_int <- function(keys, unset) {
  as.integer(env_num(keys, unset))
}

env_chr <- function(keys, unset) {
  for (key in keys) {
    value <- Sys.getenv(key, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      return(value)
    }
  }
  unset
}

env_bool <- function(keys, unset = "false") {
  tolower(env_chr(keys, unset)) %in% c("1", "true", "yes", "y", "on")
}

parse_int_seq <- function(x, default) {
  if (is.null(x) || !length(x) || is.na(x) || !nzchar(x)) {
    x <- default
  }
  x <- gsub("\\s+", "", x)
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    z <- strsplit(x, ":", fixed = TRUE)[[1]]
    return(seq.int(as.integer(z[1]), as.integer(z[2])))
  }
  vals <- as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
  vals[is.finite(vals)]
}

line_catalog <- data.frame(
  alias = c(
    "halpha", "ha", "halpha6563",
    "hbeta", "hb", "hbeta4861",
    "oiii5007", "o3", "o3_5007",
    "oiii4959", "o3_4959",
    "nii6548", "n2_6548",
    "nii6583", "n2", "n2_6583",
    "sii6716", "s2_6716",
    "sii6731", "s2_6731",
    "oi6300", "o1", "o1_6300"
  ),
  line_name = c(
    rep("Halpha", 3),
    rep("Hbeta", 3),
    rep("[OIII] 5007", 3),
    rep("[OIII] 4959", 2),
    rep("[NII] 6548", 2),
    rep("[NII] 6583", 3),
    rep("[SII] 6716", 2),
    rep("[SII] 6731", 2),
    rep("[OI] 6300", 3)
  ),
  slug = c(
    rep("halpha", 3),
    rep("hbeta", 3),
    rep("oiii5007", 3),
    rep("oiii4959", 2),
    rep("nii6548", 2),
    rep("nii6583", 3),
    rep("sii6716", 2),
    rep("sii6731", 2),
    rep("oi6300", 3)
  ),
  rest_wave = c(
    rep(6562.8, 3),
    rep(4861.33, 3),
    rep(5006.84, 3),
    rep(4958.91, 2),
    rep(6548.05, 2),
    rep(6583.45, 3),
    rep(6716.44, 2),
    rep(6730.82, 2),
    rep(6300.30, 3)
  ),
  stringsAsFactors = FALSE
)

sanitize_slug <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "line" else x
}

line_spec <- function(line_key, rest_wave_override = NA_real_) {
  key <- sanitize_slug(line_key)
  rec <- line_catalog[match(key, line_catalog$alias), ]
  if (nrow(rec) != 1L || is.na(rec$alias)) {
    label <- line_key
    slug <- sanitize_slug(line_key)
    rest <- rest_wave_override
  } else {
    label <- rec$line_name
    slug <- rec$slug
    rest <- if (is.finite(rest_wave_override)) rest_wave_override else rec$rest_wave
  }
  if (!is.finite(rest) || rest <= 0) {
    stop("Set a known CAPIVARA_LINE alias or provide CAPIVARA_LINE_REST in Angstrom.")
  }
  list(name = label, slug = slug, rest_wave = rest)
}

redshift <- env_num(c("CAPIVARA_REDSHIFT", "CAPIVARA_10218_REDSHIFT"), "0.0461")
ncomp <- env_int(c("CAPIVARA_NCOMP", "CAPIVARA_10218_NCOMP"), "25")
knn_k <- env_int(c("CAPIVARA_KNN", "CAPIVARA_10218_KNN"), "100")
path_ncomp <- env_int(c("CAPIVARA_PATH_NCOMP", "CAPIVARA_10218_PATH_NCOMP"), "45")
path_knn_k <- env_int(c("CAPIVARA_PATH_KNN", "CAPIVARA_10218_PATH_KNN"), as.character(knn_k))
path_spatial_weight <- env_num(c("CAPIVARA_PATH_SPATIAL_WEIGHT", "CAPIVARA_10218_PATH_SPATIAL_WEIGHT"), "0.10")
run_spectral_segmentation <- env_bool(c("CAPIVARA_RUN_SPECTRAL_SEGMENTATION", "CAPIVARA_10218_RUN_SPECTRAL_SEGMENTATION"), "true")
run_path_signatures <- env_bool(c("CAPIVARA_RUN_PATH_SIGNATURES", "CAPIVARA_10218_RUN_PATH_SIGNATURES"), "true")
line_key <- env_chr(c("CAPIVARA_LINE", "CAPIVARA_EMISSION_LINE", "CAPIVARA_10218_LINE"), "halpha")
line_rest_override <- env_num(c("CAPIVARA_LINE_REST", "CAPIVARA_10218_LINE_REST", "CAPIVARA_10218_HALPHA_REST"), NA_real_)
line <- line_spec(line_key, line_rest_override)
line_window_kms <- env_num(c("CAPIVARA_LINE_WINDOW_KMS", "CAPIVARA_10218_LINE_WINDOW_KMS"), "600")
cont_inner_kms <- env_num(c("CAPIVARA_LINE_CONT_INNER_KMS", "CAPIVARA_10218_LINE_CONT_INNER_KMS"), "800")
cont_outer_kms <- env_num(c("CAPIVARA_LINE_CONT_OUTER_KMS", "CAPIVARA_10218_LINE_CONT_OUTER_KMS"), "1400")
path_window_kms <- env_num(c("CAPIVARA_PATH_WINDOW_KMS", "CAPIVARA_10218_PATH_WINDOW_KMS"), as.character(line_window_kms))
centroid_window_kms <- env_num(c("CAPIVARA_CENTROID_WINDOW_KMS", "CAPIVARA_10218_CENTROID_WINDOW_KMS"), "260")
peak_search_kms <- env_num(c("CAPIVARA_PEAK_SEARCH_KMS", "CAPIVARA_10218_PEAK_SEARCH_KMS"), "350")
starlet_scales <- parse_int_seq(env_chr(c("CAPIVARA_STARLET_SCALES", "CAPIVARA_10218_STARLET_SCALES"), "2:5"), "2:5")
starlet_include_coarse <- env_bool(c("CAPIVARA_STARLET_INCLUDE_COARSE", "CAPIVARA_10218_STARLET_INCLUDE_COARSE"), "false")
object_prefix <- sanitize_slug(env_chr(c("CAPIVARA_OUTPUT_PREFIX", "CAPIVARA_10218_OUTPUT_PREFIX"), "manga10218"))
file_prefix <- paste(object_prefix, line$slug, sep = "_")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(gridExtra)
})
if (run_path_signatures && !requireNamespace("spectropath", quietly = TRUE)) {
  stop("Install spectropath or set CAPIVARA_RUN_PATH_SIGNATURES=false for native quick-mode segmentation.")
}

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_root, quiet = TRUE)
} else {
  library(capivara)
}
galaxy_mask_helpers <- file.path(repo_root, "extensions", "capivaraKinematics", "R", "galaxy_mask.R")
if (file.exists(galaxy_mask_helpers)) {
  source(galaxy_mask_helpers)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

palette_van_gogh <- function(n = 256) {
  grDevices::colorRampPalette(
    c("#80B7FF", "#547FFF", "#405CFF", "#263C8B", "#FFFAA3", "#FFDE38", "#BFA524"),
    space = "Lab"
  )(n)
}

div_palette <- grDevices::colorRampPalette(c("#2447A3", "#F7F7F7", "#BFA524"), space = "Lab")(256)
vik_palette <- c("#17335C", "#3F78A8", "#9CC7D8", "#F7F4ED", "#F1B37F", "#C85B3C", "#6F1D1B")

cluster_palette <- function(n) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    return(viridisLite::viridis(n, option = "D", end = 0.92))
  }
  grDevices::colorRampPalette(
    c("#17335C", "#3F78A8", "#4DAF8E", "#F1D66B", "#C85B3C"),
    space = "Lab"
  )(n)
}

matrix_df <- function(mat) {
  df <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  df$value <- as.vector(mat)
  df
}

robust_limits <- function(x, probs = c(0.02, 0.98), symmetric = FALSE) {
  vals <- x[is.finite(x)]
  if (!length(vals)) return(NULL)
  if (isTRUE(symmetric)) {
    lim <- stats::quantile(abs(vals), probs[2], na.rm = TRUE)
    return(c(-lim, lim))
  }
  as.numeric(stats::quantile(vals, probs, na.rm = TRUE))
}

plot_cont <- function(mat, path, title = NULL, palette = palette_van_gogh(256), limits = NULL, midpoint = NULL) {
  df <- matrix_df(mat)
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed() +
    theme_void(base_size = 11) +
    labs(title = title, fill = NULL) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  if (!is.null(midpoint)) {
    p <- p + scale_fill_gradient2(
      low = palette[[1]],
      mid = "#F7F7F7",
      high = palette[[length(palette)]],
      midpoint = midpoint,
      limits = limits,
      na.value = "black"
    )
  } else {
    p <- p + scale_fill_gradientn(colours = palette, limits = limits, na.value = "black")
  }

  ggsave(path, p, width = 5.4, height = 4.8, dpi = 320, bg = "white")
  p
}

plot_seg <- function(mat, path, title = NULL, n = NULL) {
  df <- matrix_df(mat)
  df$value <- factor(df$value)
  if (is.null(n)) {
    n <- length(unique(df$value[!is.na(df$value)]))
  }
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed() +
    scale_fill_manual(values = palette_van_gogh(max(n, 3)), na.value = "black") +
    theme_void(base_size = 11) +
    labs(title = title, fill = "bin") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  ggsave(path, p, width = 5.4, height = 4.8, dpi = 320, bg = "white")
  p
}

read_wave <- function(path, fits) {
  wave <- tryCatch(as.numeric(FITSio::readFITS(path, hdu = 6)$imDat), error = function(e) NULL)
  if (!is.null(wave) && length(wave) == dim(fits$imDat)[3]) {
    return(wave)
  }
  wave <- tryCatch(FITSio::axVec(3, fits$axDat), error = function(e) NULL)
  if (!is.null(wave) && length(wave) == dim(fits$imDat)[3]) {
    return(as.numeric(wave))
  }
  stop("Could not recover wavelength axis from FITS file.")
}

compute_line_maps <- function(cube,
                              wave,
                              mask,
                              z,
                              rest_wave,
                              line_name,
                              line_window_kms = 600,
                              peak_search_kms = 350,
                              centroid_window_kms = 260,
                              cont_inner_kms = 800,
                              cont_outer_kms = 1400) {
  lambda0 <- rest_wave * (1 + z)
  vel <- 299792.458 * (wave / lambda0 - 1)
  line_idx <- which(abs(vel) <= line_window_kms)
  cont_idx <- which(abs(vel) > cont_inner_kms & abs(vel) <= cont_outer_kms)
  if (length(line_idx) < 5L || length(cont_idx) < 5L) {
    stop("Insufficient wavelength channels for ", line_name, " line/continuum windows.")
  }

  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  flux <- matrix(NA_real_, nx, ny)
  velocity <- matrix(NA_real_, nx, ny)
  sigma <- matrix(NA_real_, nx, ny)
  asymmetry <- matrix(NA_real_, nx, ny)
  h3_proxy <- matrix(NA_real_, nx, ny)
  h4_proxy <- matrix(NA_real_, nx, ny)

  for (idx in which(mask)) {
    ij <- arrayInd(idx, .dim = c(nx, ny))
    i <- ij[1]
    j <- ij[2]
    spec <- cube[i, j, ]
    if (!any(is.finite(spec[line_idx]))) next

    cont <- stats::median(spec[cont_idx], na.rm = TRUE)
    if (!is.finite(cont)) next
    line <- spec[line_idx] - cont
    line[!is.finite(line)] <- 0
    v <- vel[line_idx]
    search <- abs(v) <= peak_search_kms
    if (sum(search) < 3L) {
      search <- rep(TRUE, length(v))
    }
    search_idx <- which(search)
    peak <- search_idx[which.max(line[search])]
    if (!length(peak) || !is.finite(line[peak]) || line[peak] <= 0) next

    local <- abs(v - v[peak]) <= centroid_window_kms
    if (sum(local) < 3L) {
      local <- rep(TRUE, length(v))
    }
    pos <- pmax(line, 0)
    pos[!local] <- 0
    sum_pos <- sum(pos)
    if (!is.finite(sum_pos) || sum_pos <= 0) next

    mu <- sum(v * pos) / sum_pos
    sig <- sqrt(sum(pos * (v - mu)^2) / sum_pos)
    blue <- sum(pos[v < 0])
    red <- sum(pos[v > 0])

    flux[i, j] <- sum_pos
    velocity[i, j] <- mu
    sigma[i, j] <- sig
    asymmetry[i, j] <- (red - blue) / (red + blue + .Machine$double.eps)
    if (is.finite(sig) && sig > 0) {
      zvel <- (v - mu) / sig
      h3_proxy[i, j] <- sum(pos * zvel^3) / sum_pos
      h4_proxy[i, j] <- sum(pos * zvel^4) / sum_pos - 3
    }
  }

  ok <- mask & is.finite(flux) & flux > 0 & is.finite(velocity) & is.finite(sigma)
  if (sum(ok) > 5L) {
    velocity[ok] <- velocity[ok] - stats::median(velocity[ok], na.rm = TRUE)
  }

  list(
    lambda0 = lambda0,
    line_name = line_name,
    rest_wave = rest_wave,
    line_window_kms = line_window_kms,
    cont_inner_kms = cont_inner_kms,
    cont_outer_kms = cont_outer_kms,
    velocity_grid = vel,
    line_idx = line_idx,
    cont_idx = cont_idx,
    flux = flux,
    velocity = velocity,
    sigma = sigma,
    asymmetry = asymmetry,
    h3_proxy = h3_proxy,
    h4_proxy = h4_proxy,
    valid = ok
  )
}

baseline_subtract <- function(v, y, n_edge = 2L) {
  edge <- c(seq_len(min(n_edge, length(y))), seq.int(max(1L, length(y) - n_edge + 1L), length(y)))
  base <- stats::median(y[edge], na.rm = TRUE)
  y - base
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
  qv <- stats::approx(cum, v[ord], xout = c(0.1, 0.5, 0.9), ties = "ordered", rule = 2)$y

  data.frame(
    np_flux = flux,
    np_peak_v = v[which.max(line)],
    np_centroid = centroid,
    np_sigma = sigma,
    np_w80 = qv[3] - qv[1]
  )
}

nearest_fill_map <- function(mat, support, source_mask) {
  out <- mat
  source <- support & source_mask & is.finite(mat)
  missing <- support & !source
  if (!any(missing) || !any(source)) {
    return(out)
  }
  source_idx <- which(source, arr.ind = TRUE)
  missing_idx <- which(missing, arr.ind = TRUE)
  for (k in seq_len(nrow(missing_idx))) {
    d2 <- (source_idx[, 1] - missing_idx[k, 1])^2 + (source_idx[, 2] - missing_idx[k, 2])^2
    hit <- which.min(d2)
    out[missing_idx[k, 1], missing_idx[k, 2]] <- mat[source_idx[hit, 1], source_idx[hit, 2]]
  }
  out
}

fill_kinematic_holes <- function(kin, support) {
  measured_valid <- kin$valid
  kin$flux <- nearest_fill_map(kin$flux, support, measured_valid)
  kin$velocity <- nearest_fill_map(kin$velocity, support, measured_valid)
  kin$sigma <- nearest_fill_map(kin$sigma, support, measured_valid)
  kin$asymmetry <- nearest_fill_map(kin$asymmetry, support, measured_valid)
  kin$h3_proxy <- nearest_fill_map(kin$h3_proxy, support, measured_valid)
  kin$h4_proxy <- nearest_fill_map(kin$h4_proxy, support, measured_valid)
  kin$measured_valid <- measured_valid
  kin$valid <- support &
    is.finite(kin$flux) &
    is.finite(kin$velocity) &
    is.finite(kin$sigma) &
    is.finite(kin$h3_proxy) &
    is.finite(kin$h4_proxy)
  kin$imputed <- kin$valid & !measured_valid
  kin
}

nearest_impute_feature_cube <- function(feature_cube, mask) {
  dims <- dim(feature_cube)
  mat <- matrix(feature_cube, nrow = dims[1] * dims[2], ncol = dims[3])
  xy <- expand.grid(x = seq_len(dims[1]), y = seq_len(dims[2]))
  support <- as.vector(mask)
  complete <- support & apply(mat, 1, function(z) all(is.finite(z)))
  missing <- support & !complete

  if (!any(complete)) {
    stop("No complete path-feature spaxels were available for clustering.")
  }

  if (any(missing)) {
    source_rows <- which(complete)
    source_xy <- as.matrix(xy[source_rows, c("x", "y")])
    for (row in which(missing)) {
      d2 <- (source_xy[, 1] - xy$x[row])^2 + (source_xy[, 2] - xy$y[row])^2
      mat[row, ] <- mat[source_rows[which.min(d2)], ]
    }
  }

  feature_cube[] <- mat
  feature_cube
}

build_path_feature_cube <- function(cube,
                                    wave,
                                    observed_wave,
                                    mask,
                                    max_abs_velocity = 600,
                                    feature_names = c("p2", "p3u", "p3F", "p4F", "p4T", "p_pm")) {
  vel <- 299792.458 * (wave / observed_wave - 1)
  keep <- abs(vel) <= max_abs_velocity
  vel <- vel[keep]
  line_cube <- cube[, , keep, drop = FALSE]
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  feature_cube <- array(NA_real_, dim = c(nx, ny, length(feature_names)), dimnames = list(NULL, NULL, feature_names))
  rows <- vector("list", sum(mask, na.rm = TRUE))
  n <- 0L

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      if (!isTRUE(mask[i, j])) next
      flux <- as.numeric(line_cube[i, j, ])
      if (sum(is.finite(flux)) < 8L) next
      flux[!is.finite(flux)] <- 0
      line <- baseline_subtract(vel, flux)
      amp <- max(abs(line), na.rm = TRUE)
      if (!is.finite(amp) || amp <= 0) next
      profile <- line / amp
      pf <- tryCatch(
        spectropath::path_features(cbind(vel, profile), depth = 4, normalize = TRUE, notation = "paper"),
        error = function(e) NULL
      )
      if (is.null(pf)) next
      vals <- as.numeric(pf[1, feature_names, drop = TRUE])
      feature_cube[i, j, ] <- vals

      np <- nonparam_velocity(vel, line)
      n <- n + 1L
      rows[[n]] <- data.frame(
        x = i,
        y = j,
        line_flux = np$np_flux,
        centroid_kms = np$np_centroid,
        sigma_kms = np$np_sigma,
        w80_kms = np$np_w80,
        pf[, feature_names, drop = FALSE],
        check.names = FALSE
      )
    }
  }

  feature_cube <- nearest_impute_feature_cube(feature_cube, mask)
  table <- if (n) do.call(rbind, rows[seq_len(n)]) else data.frame()

  list(
    feature_cube = feature_cube,
    table = table,
    velocity = vel,
    line_cube = line_cube,
    features = feature_names
  )
}

segment_median_map <- function(seg_map, value_map) {
  out <- matrix(NA_real_, nrow = nrow(seg_map), ncol = ncol(seg_map))
  ids <- sort(unique(as.integer(seg_map[is.finite(seg_map)])))
  for (id in ids) {
    pix <- seg_map == id
    vals <- value_map[pix]
    out[pix] <- stats::median(vals[is.finite(vals)], na.rm = TRUE)
  }
  out
}

plot_path_segment <- function(seg_map, path, title = NULL) {
  df <- matrix_df(seg_map)
  ids <- sort(unique(as.integer(df$value[is.finite(df$value)])))
  df$value <- factor(df$value, levels = ids)
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_manual(values = cluster_palette(max(length(ids), 3)), na.value = "white", guide = "none") +
    theme_void(base_size = 11) +
    labs(title = title) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.02),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  ggsave(path, p, width = 5.4, height = 4.8, dpi = 320, bg = "white")
  p
}

plot_velocity_colored_segments <- function(seg_map, velocity_map, path, title = NULL) {
  med_map <- segment_median_map(seg_map, velocity_map)
  df <- matrix_df(med_map)
  lim <- robust_limits(df$value, symmetric = TRUE)
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    scale_fill_gradient2(
      low = vik_palette[1],
      mid = vik_palette[4],
      high = vik_palette[7],
      midpoint = 0,
      limits = lim,
      oob = scales::squish,
      na.value = "white"
    ) +
    theme_void(base_size = 11) +
    labs(title = title, fill = "km/s") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.02),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position = "right"
    )
  ggsave(path, p, width = 5.4, height = 4.8, dpi = 320, bg = "white")
  p
}

message("Reading cube: ", cube_path)
message(sprintf(
  "Using %s rest=%.3f A, z=%.5f, observed=%.3f A",
  line$name,
  line$rest_wave,
  redshift,
  line$rest_wave * (1 + redshift)
))
fits <- FITSio::readFITS(cube_path, hdu = 1)
wave <- read_wave(cube_path, fits)
cube <- fits$imDat

message("Building full-frame starlet support...")
star <- build_starlet_mask(
  fits,
  starlet_J = 5,
  starlet_scales = starlet_scales,
  include_coarse = starlet_include_coarse,
  denoise_k = 0,
  positive_only = TRUE
)
support_raw <- star$mask
support <- support_raw
if (exists("better_galaxy_mask", mode = "function")) {
  support <- better_galaxy_mask(support_raw, close_iterations = 1L, fill_holes = TRUE, preserve_input = TRUE, connectivity = 8L)
}

seg <- list(
  cluster_map = matrix(NA_integer_, nrow = dim(cube)[1], ncol = dim(cube)[2]),
  backend_info = list(actual_Ncomp = NA_integer_)
)
if (run_spectral_segmentation) {
  message("Running full-spectrum Capivara segmentation...")
  seg <- segment_large(
    fits,
    Ncomp = ncomp,
    redshift = redshift,
    use_starlet_mask = FALSE,
    mask = support,
    knn_k = knn_k,
    auto_k = FALSE,
    max_k = knn_k,
    verbose = TRUE
  )
} else {
  message("Skipping full-spectrum Capivara segmentation; using native kinematic maps only.")
}

message("Computing ", line$name, " kinematic feature maps...")
kin <- compute_line_maps(
  cube,
  wave,
  support,
  redshift,
  rest_wave = line$rest_wave,
  line_name = line$name,
  line_window_kms = line_window_kms,
  peak_search_kms = peak_search_kms,
  centroid_window_kms = centroid_window_kms,
  cont_inner_kms = cont_inner_kms,
  cont_outer_kms = cont_outer_kms
)
kin <- fill_kinematic_holes(kin, support)

kin_cube <- array(NA_real_, dim = c(dim(cube)[1], dim(cube)[2], 6L))
kin_cube[, , 1] <- log10(kin$flux)
kin_cube[, , 2] <- kin$velocity
kin_cube[, , 3] <- kin$sigma
kin_cube[, , 4] <- kin$asymmetry
kin_cube[, , 5] <- kin$h3_proxy
kin_cube[, , 6] <- kin$h4_proxy
kin_input <- list(imDat = kin_cube, hdr = fits$hdr, axDat = NULL)

message("Running kinematic-aware Capivara segmentation...")
kin_seg <- segment_large(
  kin_input,
  Ncomp = ncomp,
  scale_fn = identity,
  knn_k = knn_k,
  auto_k = FALSE,
  max_k = knn_k,
  feature_scale = "robust_col",
  mask = kin$valid,
  valid_mode = "finite",
  verbose = TRUE
)

path_features <- NULL
path_seg <- NULL
if (run_path_signatures) {
  message("Computing ", line$name, " path-signature features...")
  path_features <- build_path_feature_cube(
    cube = cube,
    wave = wave,
    observed_wave = kin$lambda0,
    mask = kin$valid,
    max_abs_velocity = path_window_kms
  )
  path_input <- list(imDat = path_features$feature_cube, hdr = fits$hdr, axDat = NULL)

  message("Running path-signature kinematic-aware Capivara segmentation...")
  path_seg <- segment_large(
    path_input,
    Ncomp = path_ncomp,
    scale_fn = identity,
    knn_k = path_knn_k,
    auto_k = FALSE,
    max_k = path_knn_k,
    feature_scale = "robust_col",
    spatial_weight = path_spatial_weight,
    mask = kin$valid,
    valid_mode = "finite",
    verbose = TRUE
  )
  if (nrow(path_features$table)) {
    path_features$table$path_signature_segment <- path_seg$cluster_map[cbind(path_features$table$x, path_features$table$y)]
  }
}

message("Saving products...")
plot_cont(star$collapsed, file.path(out_dir, paste0(object_prefix, "_white_light.png")), "white light")
plot_cont(ifelse(support_raw, 1, NA_real_), file.path(out_dir, paste0(object_prefix, "_starlet_support_raw.png")), sprintf("raw starlet support: %d spaxels", sum(support_raw)), palette = c("#F2D06B", "#F2D06B"), limits = c(0, 1))
plot_cont(ifelse(support, 1, NA_real_), file.path(out_dir, paste0(object_prefix, "_starlet_support.png")), sprintf("better starlet mask: %d spaxels", sum(support)), palette = c("#F2D06B", "#F2D06B"), limits = c(0, 1))
if (run_spectral_segmentation) {
  plot_seg(seg$cluster_map, file.path(out_dir, sprintf("%s_capivara_segments_n%d.png", object_prefix, ncomp)), sprintf("Capivara full-spectrum segments (N=%d)", ncomp), n = ncomp)
}
plot_cont(log10(kin$flux), file.path(out_dir, paste0(file_prefix, "_flux_log.png")), paste(line$name, "log flux"), limits = robust_limits(log10(kin$flux)))
plot_cont(kin$velocity, file.path(out_dir, paste0(file_prefix, "_velocity_centered.png")), paste(line$name, "velocity, median centered"), palette = div_palette, limits = robust_limits(kin$velocity, symmetric = TRUE), midpoint = 0)
plot_cont(kin$sigma, file.path(out_dir, paste0(file_prefix, "_sigma.png")), paste(line$name, "sigma"), limits = robust_limits(kin$sigma))
plot_cont(kin$asymmetry, file.path(out_dir, paste0(file_prefix, "_asymmetry.png")), paste(line$name, "red-blue asymmetry"), palette = div_palette, limits = c(-1, 1), midpoint = 0)
plot_cont(kin$h3_proxy, file.path(out_dir, paste0(file_prefix, "_h3_proxy.png")), paste(line$name, "h3 proxy"), palette = div_palette, limits = robust_limits(kin$h3_proxy, symmetric = TRUE), midpoint = 0)
plot_cont(kin$h4_proxy, file.path(out_dir, paste0(file_prefix, "_h4_proxy.png")), paste(line$name, "h4 proxy"), palette = div_palette, limits = robust_limits(kin$h4_proxy, symmetric = TRUE), midpoint = 0)
plot_cont(ifelse(kin$imputed, 1, NA_real_), file.path(out_dir, paste0(file_prefix, "_imputed_line_holes.png")), sprintf("%s nearest-filled holes: %d spaxels", line$name, sum(kin$imputed)), palette = c("#C85B3C", "#C85B3C"), limits = c(0, 1))
plot_seg(kin_seg$cluster_map, file.path(out_dir, sprintf("%s_%s_kinematic_aware_segments_n%d.png", object_prefix, line$slug, ncomp)), sprintf("%s kinematic-aware segments (N=%d)", line$name, ncomp), n = ncomp)
if (run_path_signatures) {
  plot_path_segment(
    path_seg$cluster_map,
    file.path(out_dir, sprintf("%s_path_signature_segments_n%d.png", file_prefix, path_ncomp)),
    sprintf("%s path-signature segments (N=%d)", line$name, path_ncomp)
  )
  plot_velocity_colored_segments(
    path_seg$cluster_map,
    kin$velocity,
    file.path(out_dir, sprintf("%s_path_signature_velocity_segments_n%d.png", file_prefix, path_ncomp)),
    sprintf("path-signature segments, median %s velocity (N=%d)", line$name, path_ncomp)
  )
}

segment_panel_plot <- if (run_path_signatures) {
  plot_velocity_colored_segments(path_seg$cluster_map, kin$velocity, file.path(out_dir, "_tmp_path_velocity_segments.png"), "path-aware segments")
} else {
  plot_velocity_colored_segments(kin_seg$cluster_map, kin$velocity, file.path(out_dir, "_tmp_kin_velocity_segments.png"), "kinematic segments")
}
spectral_panel_plot <- if (run_spectral_segmentation) {
  plot_seg(seg$cluster_map, file.path(out_dir, "_tmp_segments.png"), "Capivara spectral", n = ncomp)
} else {
  plot_cont(log10(kin$flux), file.path(out_dir, "_tmp_line_flux_panel.png"), paste(line$name, "flux"), limits = robust_limits(log10(kin$flux)))
}

panel <- gridExtra::arrangeGrob(
  plot_cont(star$collapsed, file.path(out_dir, "_tmp_white_light.png"), "white light"),
  spectral_panel_plot,
  plot_cont(kin$velocity, file.path(out_dir, "_tmp_velocity.png"), paste(line$name, "velocity"), palette = vik_palette, limits = robust_limits(kin$velocity, symmetric = TRUE), midpoint = 0),
  segment_panel_plot,
  ncol = 4
)
ggsave(file.path(out_dir, paste0(file_prefix, "_capivara_kinematic_panel.png")), panel, width = 14, height = 3.5, dpi = 320, bg = "white")

if (run_path_signatures) {
  path_panel <- gridExtra::arrangeGrob(
    plot_cont(log10(kin$flux), file.path(out_dir, "_tmp_line_flux.png"), paste(line$name, "flux"), limits = robust_limits(log10(kin$flux))),
    plot_cont(kin$velocity, file.path(out_dir, "_tmp_line_velocity.png"), paste(line$name, "velocity"), palette = vik_palette, limits = robust_limits(kin$velocity, symmetric = TRUE), midpoint = 0),
    plot_path_segment(path_seg$cluster_map, file.path(out_dir, "_tmp_path_segments.png"), "path signatures"),
    plot_velocity_colored_segments(path_seg$cluster_map, kin$velocity, file.path(out_dir, "_tmp_path_velocity.png"), "velocity-colored path groups"),
    ncol = 4
  )
  ggsave(file.path(out_dir, sprintf("%s_path_signature_pretty_panel_n%d.png", file_prefix, path_ncomp)), path_panel, width = 14, height = 3.5, dpi = 320, bg = "white")
}

tab <- expand.grid(x = seq_len(dim(cube)[1]), y = seq_len(dim(cube)[2]))
tab$starlet_support <- as.vector(support)
tab$capivara_segment <- as.vector(seg$cluster_map)
tab[[paste0(line$slug, "_flux")]] <- as.vector(kin$flux)
tab[[paste0(line$slug, "_velocity_centered")]] <- as.vector(kin$velocity)
tab[[paste0(line$slug, "_sigma")]] <- as.vector(kin$sigma)
tab[[paste0(line$slug, "_asymmetry")]] <- as.vector(kin$asymmetry)
tab[[paste0(line$slug, "_h3_proxy")]] <- as.vector(kin$h3_proxy)
tab[[paste0(line$slug, "_h4_proxy")]] <- as.vector(kin$h4_proxy)
tab[[paste0(line$slug, "_measured_valid")]] <- as.vector(kin$measured_valid)
tab[[paste0(line$slug, "_imputed")]] <- as.vector(kin$imputed)
tab$kinematic_aware_segment <- as.vector(kin_seg$cluster_map)
if (run_path_signatures) {
  tab$path_signature_segment <- as.vector(path_seg$cluster_map)
}
utils::write.csv(tab, file.path(out_dir, paste0(file_prefix, "_capivara_kinematic_spaxel_table.csv")), row.names = FALSE)
if (run_path_signatures) {
  utils::write.csv(path_features$table, file.path(out_dir, paste0(file_prefix, "_path_signature_spaxel_features.csv")), row.names = FALSE)
}

maps <- list(
  starlet_support = support + 0,
  capivara_segment = seg$cluster_map,
  setNames(list(log10(kin$flux)), paste0(line$slug, "_log_flux"))[[1]],
  setNames(list(kin$velocity), paste0(line$slug, "_velocity_centered"))[[1]],
  setNames(list(kin$sigma), paste0(line$slug, "_sigma"))[[1]],
  setNames(list(kin$asymmetry), paste0(line$slug, "_asymmetry"))[[1]],
  setNames(list(kin$h3_proxy), paste0(line$slug, "_h3_proxy"))[[1]],
  setNames(list(kin$h4_proxy), paste0(line$slug, "_h4_proxy"))[[1]],
  kinematic_aware_segment = kin_seg$cluster_map
)
names(maps)[3:8] <- c(
  paste0(line$slug, "_log_flux"),
  paste0(line$slug, "_velocity_centered"),
  paste0(line$slug, "_sigma"),
  paste0(line$slug, "_asymmetry"),
  paste0(line$slug, "_h3_proxy"),
  paste0(line$slug, "_h4_proxy")
)
if (run_path_signatures) {
  maps$path_signature_segment <- path_seg$cluster_map
  for (k in seq_along(path_features$features)) {
    maps[[paste0("path_", path_features$features[[k]])]] <- path_features$feature_cube[, , k]
  }
}
map_stack <- array(NA_real_, dim = c(dim(cube)[1], dim(cube)[2], length(maps)))
for (k in seq_along(maps)) {
  map_stack[, , k] <- maps[[k]]
}
try(FITSio::writeFITSim(map_stack, file.path(out_dir, paste0(file_prefix, "_capivara_kinematic_maps.fits")), type = "double"), silent = TRUE)
utils::write.csv(
  data.frame(
    channel = seq_len(dim(map_stack)[3]),
    name = names(maps)
  ),
  file.path(out_dir, paste0(file_prefix, "_fits_channels.csv")),
  row.names = FALSE
)

saveRDS(
  list(
    cube_path = cube_path,
    redshift = redshift,
    ncomp = ncomp,
    knn_k = knn_k,
    run_spectral_segmentation = run_spectral_segmentation,
    run_path_signatures = run_path_signatures,
    path_ncomp = path_ncomp,
    path_knn_k = path_knn_k,
    path_spatial_weight = path_spatial_weight,
    line = line,
    line_observed = kin$lambda0,
    line_window_kms = line_window_kms,
    cont_inner_kms = cont_inner_kms,
    cont_outer_kms = cont_outer_kms,
    peak_search_kms = peak_search_kms,
    centroid_window_kms = centroid_window_kms,
    path_window_kms = path_window_kms,
    starlet_scales = starlet_scales,
    starlet_include_coarse = starlet_include_coarse,
    starlet = star,
    support_raw = support_raw,
    support = support,
    capivara = seg,
    kinematics = kin,
    kinematic_features = kin_cube,
    kinematic_aware = kin_seg,
    path_features = path_features,
    path_signature = path_seg
  ),
  file.path(out_dir, paste0(file_prefix, "_capivara_kinematic_results.rds"))
)

unlink(file.path(out_dir, c(
  "_tmp_white_light.png",
  "_tmp_segments.png",
  "_tmp_line_flux_panel.png",
  "_tmp_velocity.png",
  "_tmp_kin_velocity_segments.png",
  "_tmp_path_velocity_segments.png",
  "_tmp_line_flux.png",
  "_tmp_line_velocity.png",
  "_tmp_path_segments.png",
  "_tmp_path_velocity.png"
)))

message("Wrote outputs to: ", out_dir)
if (run_path_signatures) {
  message(sprintf(
    "Summary: support=%d, spectral valid=%d, kin valid=%d, path valid=%d, spectral actual=%d, kin actual=%d, path actual=%d",
    sum(support),
    sum(!is.na(seg$cluster_map)),
    sum(kin$valid),
    sum(!is.na(path_seg$cluster_map)),
    seg$backend_info$actual_Ncomp,
    kin_seg$backend_info$actual_Ncomp,
    path_seg$backend_info$actual_Ncomp
  ))
} else {
  message(sprintf(
    "Summary: support=%d, spectral valid=%d, kin valid=%d, spectral actual=%d, kin actual=%d, path signatures skipped",
    sum(support),
    sum(!is.na(seg$cluster_map)),
    sum(kin$valid),
    seg$backend_info$actual_Ncomp,
    kin_seg$backend_info$actual_Ncomp
  ))
}
