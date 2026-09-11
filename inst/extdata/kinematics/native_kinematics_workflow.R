args <- commandArgs(trailingOnly = TRUE)
use_env_inputs <- tolower(Sys.getenv("CAPIVARA_USE_ENV_INPUTS", unset = "false")) %in% c("1", "true", "yes", "y", "on")

# Internal workflow called by run_manga_bar_model().
# It is bundled with capivara so installed users do not need a source checkout.
#
# Common knobs:
#   CAPIVARA_REDSHIFT=0.0461
#   CAPIVARA_LINE=halpha            # aliases include hbeta, oiii5007, nii6583, sii6716
#   CAPIVARA_LINE_REST=6564.608     # optional registry-consistent vacuum Angstrom
#   CAPIVARA_OUTPUT_PREFIX=manga10218
#   CAPIVARA_KNN=100 CAPIVARA_NCOMP=25
#   CAPIVARA_PATH_KNN=100 CAPIVARA_PATH_NCOMP=45 CAPIVARA_PATH_SPATIAL_WEIGHT=0.10
#
# If [out_dir] is omitted, products are written beside the cube in:
#   dirname(cube_path)/capivara_outputs

cube_path <- if (!use_env_inputs && length(args) >= 1) {
  args[[1]]
} else {
  env_cube_path <- Sys.getenv("CAPIVARA_CUBE_PATH", unset = "")
  if (nzchar(env_cube_path)) env_cube_path else stop("Supply a cube path explicitly.", call. = FALSE)
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

sanitize_slug <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "line" else x
}

line_spec <- function(line_key, rest_wave_override = NA_real_, medium = "vacuum") {
  rec <- .capivara_match_emission_lines(line_key, medium)
  if (nrow(rec) != 1L) stop("Choose one registered emission line.")
  if (is.finite(rest_wave_override) && abs(rest_wave_override-rec$rest_wavelength)>1e-6) {
    stop("CAPIVARA_LINE_REST disagrees with the referenced line registry/medium.")
  }
  list(name = rec$label, slug = rec$name, rest_wave = rec$rest_wavelength)
}

redshift <- env_num(c("CAPIVARA_REDSHIFT", "CAPIVARA_10218_REDSHIFT"), NA_real_)
ncomp <- env_int(c("CAPIVARA_NCOMP", "CAPIVARA_10218_NCOMP"), "25")
knn_k <- env_int(c("CAPIVARA_KNN", "CAPIVARA_10218_KNN"), "100")
path_ncomp <- env_int(c("CAPIVARA_PATH_NCOMP", "CAPIVARA_10218_PATH_NCOMP"), "45")
path_knn_k <- env_int(c("CAPIVARA_PATH_KNN", "CAPIVARA_10218_PATH_KNN"), as.character(knn_k))
path_spatial_weight <- env_num(c("CAPIVARA_PATH_SPATIAL_WEIGHT", "CAPIVARA_10218_PATH_SPATIAL_WEIGHT"), "0.10")
run_spectral_segmentation <- env_bool(c("CAPIVARA_RUN_SPECTRAL_SEGMENTATION", "CAPIVARA_10218_RUN_SPECTRAL_SEGMENTATION"), "true")
run_path_signatures <- env_bool(c("CAPIVARA_RUN_PATH_SIGNATURES", "CAPIVARA_10218_RUN_PATH_SIGNATURES"), "true")
line_key <- env_chr(c("CAPIVARA_LINE", "CAPIVARA_EMISSION_LINE", "CAPIVARA_10218_LINE"), "halpha")
line_rest_override <- env_num(c("CAPIVARA_LINE_REST", "CAPIVARA_10218_LINE_REST", "CAPIVARA_10218_HALPHA_REST"), NA_real_)
wavelength_frame <- env_chr("CAPIVARA_WAVELENGTH_FRAME", "observed")
wavelength_medium <- env_chr("CAPIVARA_WAVELENGTH_MEDIUM", "vacuum")
systemic_redshift_source <- env_chr("CAPIVARA_REDSHIFT_SOURCE", "explicit workflow environment")
profile_centering_mode <- env_chr("CAPIVARA_PROFILE_CENTERING_MODE", "systemic")
line <- line_spec(line_key, line_rest_override, wavelength_medium)
line_window_kms <- env_num(c("CAPIVARA_LINE_WINDOW_KMS", "CAPIVARA_10218_LINE_WINDOW_KMS"), "600")
cont_inner_kms <- env_num(c("CAPIVARA_LINE_CONT_INNER_KMS", "CAPIVARA_10218_LINE_CONT_INNER_KMS"), "800")
cont_outer_kms <- env_num(c("CAPIVARA_LINE_CONT_OUTER_KMS", "CAPIVARA_10218_LINE_CONT_OUTER_KMS"), "1400")
path_window_kms <- env_num(c("CAPIVARA_PATH_WINDOW_KMS", "CAPIVARA_10218_PATH_WINDOW_KMS"), as.character(line_window_kms))
centroid_window_kms <- env_num(c("CAPIVARA_CENTROID_WINDOW_KMS", "CAPIVARA_10218_CENTROID_WINDOW_KMS"), "260")
peak_search_kms <- env_num(c("CAPIVARA_PEAK_SEARCH_KMS", "CAPIVARA_10218_PEAK_SEARCH_KMS"), "350")
starlet_scales <- parse_int_seq(env_chr(c("CAPIVARA_STARLET_SCALES", "CAPIVARA_10218_STARLET_SCALES"), "2:5"), "2:5")
starlet_include_coarse <- env_bool(c("CAPIVARA_STARLET_INCLUDE_COARSE", "CAPIVARA_10218_STARLET_INCLUDE_COARSE"), "false")
kinematic_support <- tolower(env_chr(c("CAPIVARA_KINEMATIC_SUPPORT", "CAPIVARA_10218_KINEMATIC_SUPPORT"), "starlet"))
if (!kinematic_support %in% c("starlet", "line_flux")) {
  stop("CAPIVARA_KINEMATIC_SUPPORT must be 'starlet' or 'line_flux'.", call. = FALSE)
}
line_flux_sigma <- env_num(c("CAPIVARA_KINEMATIC_LINE_FLUX_SIGMA", "CAPIVARA_10218_KINEMATIC_LINE_FLUX_SIGMA"), "3")
if (!is.finite(line_flux_sigma) || line_flux_sigma <= 0) line_flux_sigma <- 3
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

if (!requireNamespace("capivara", quietly = TRUE)) {
  stop("Install the `capivara` package before running the kinematics workflow.", call. = FALSE)
}
library(capivara)

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
  wave <- tryCatch(as.numeric(.capivara_read_fits(path, hdu = 6)$imDat), error = function(e) NULL)
  if (!is.null(wave) && length(wave) == dim(fits$imDat)[3]) {
    return(wave)
  }
  wave <- tryCatch(FITSio::axVec(3, fits$axDat), error = function(e) NULL)
  if (!is.null(wave) && length(wave) == dim(fits$imDat)[3]) {
    return(as.numeric(wave))
  }
  stop("Could not recover wavelength axis from FITS file.")
}

clean_support <- function(mask) {
  if (!exists("clean_galaxy_support", mode = "function", inherits = TRUE)) {
    return(mask)
  }
  clean_galaxy_support(
    mask,
    close_iterations = 1L,
    fill_holes = TRUE,
    preserve_input = FALSE,
    connectivity = 8L
  )
}

line_flux_support <- function(flux, candidate_mask, z_threshold = 3) {
  candidate_mask <- is.finite(candidate_mask) & candidate_mask
  nr <- nrow(flux)
  nc <- ncol(flux)
  border <- matrix(FALSE, nr, nc)
  br <- max(1L, floor(0.10 * nr))
  bc <- max(1L, floor(0.10 * nc))
  border[seq_len(br), ] <- TRUE
  border[(nr - br + 1L):nr, ] <- TRUE
  border[, seq_len(bc)] <- TRUE
  border[, (nc - bc + 1L):nc] <- TRUE

  sky <- flux[border & candidate_mask & is.finite(flux)]
  if (length(sky) < 20L) {
    values <- flux[candidate_mask & is.finite(flux)]
    if (!length(values)) {
      stop("No finite line flux is available for line-flux support.", call. = FALSE)
    }
    sky <- values[values <= stats::quantile(values, 0.25, na.rm = TRUE)]
  }
  center <- stats::median(sky, na.rm = TRUE)
  scale <- stats::mad(sky, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    q <- stats::quantile(sky, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    scale <- diff(q) / 1.349
  }
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(sky, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    stop("Could not estimate border noise for line-flux support.", call. = FALSE)
  }

  threshold <- center + z_threshold * scale
  raw_mask <- candidate_mask & is.finite(flux) & flux > threshold
  mask <- clean_support(raw_mask)
  if (sum(mask, na.rm = TRUE) < 6L) {
    stop("Line-flux support is too small; lower `line_flux_sigma` or use starlet support.", call. = FALSE)
  }
  list(
    mask = mask,
    raw_mask = raw_mask,
    threshold = threshold,
    center = center,
    scale = scale,
    z_threshold = z_threshold
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
fits <- .capivara_read_fits(cube_path, hdu = 1)
wavelength_frame <- .input_wavelength_frame(fits, wavelength_frame, required = TRUE)
if (!is.null(fits$wavelength_medium) && fits$wavelength_medium != wavelength_medium) {
  stop("The requested wavelength_medium conflicts with the native MaNGA metadata.")
}
wave <- read_wave(cube_path, fits)
lsf <- read_manga_lsf(cube_path, fits)
lsf_provenance <- c(lsf$provenance, list(
  fitting_method = "native profile moments and SpectroPath descriptors",
  selected_lsf = "neither: observed profile representation",
  correction_applied = FALSE, width_interpretation = "observed; includes instrumental and pixel broadening"))
# Retain both native LSFs at every extracted line sample in the returned product.
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
support_starlet <- clean_support(support_raw)
support <- support_starlet
support_info <- list(method = "starlet")
if (identical(kinematic_support, "line_flux")) {
  message("Building line-flux kinematic support...")
  preliminary_kin <- compute_line_maps(
    cube,
    wave,
    support_starlet,
    redshift,
    rest_wave = line$rest_wave,
    line_name = line$name,
    line_window_kms = line_window_kms,
    peak_search_kms = peak_search_kms,
    centroid_window_kms = centroid_window_kms,
    cont_inner_kms = cont_inner_kms,
    cont_outer_kms = cont_outer_kms,
    wavelength_frame = wavelength_frame, wavelength_medium = wavelength_medium,
    systemic_redshift_source = systemic_redshift_source
  )
  support_info <- line_flux_support(
    preliminary_kin$flux,
    support_starlet,
    z_threshold = line_flux_sigma
  )
  support_info$method <- "line_flux"
  support <- support_info$mask
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
  cont_outer_kms = cont_outer_kms,
  wavelength_frame = wavelength_frame, wavelength_medium = wavelength_medium,
  systemic_redshift_source = systemic_redshift_source
)
# Preserve measured systemic velocities. The model runner can explicitly flag
# preview imputation; segmentation here uses only measured line profiles.
kin$measured_valid <- kin$valid
kin$imputed <- matrix(FALSE, nrow(support), ncol(support))

kin_cube <- array(NA_real_, dim = c(dim(cube)[1], dim(cube)[2], 6L))
kin_cube[, , 1] <- log10(kin$flux)
kin_cube[, , 2] <- kin$velocity
kin_cube[, , 3] <- kin$sigma
kin_cube[, , 4] <- kin$asymmetry
kin_cube[, , 5] <- kin$h3_proxy
kin_cube[, , 6] <- kin$h4_proxy
kin$frame_provenance$lsf <- lsf_provenance
kin$lsf_sigma_angstrom_pre <- lsf$lsf_sigma_angstrom_pre[,,kin$line_idx,drop=FALSE]
kin$lsf_sigma_angstrom_post <- lsf$lsf_sigma_angstrom_post[,,kin$line_idx,drop=FALSE]
rm(lsf)
gc()
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

kin_seg$feature_axis_provenance <- kin_seg$wavelength_provenance
kin_seg$wavelength_provenance <- kin$wavelength_provenance
kin_seg$kinematic_provenance <- kin$frame_provenance

path_features <- NULL
path_seg <- NULL
if (run_path_signatures) {
  message("Computing ", line$name, " path-signature features...")
  path_features <- build_path_feature_cube(
    cube = cube,
    wave = wave,
    observed_wave = kin$lambda0,
    mask = kin$valid,
    max_abs_velocity = path_window_kms,
    rest_wave = line$rest_wave, redshift = redshift, line_name = line$name,
    wavelength_frame = wavelength_frame, wavelength_medium = wavelength_medium,
    systemic_redshift_source = systemic_redshift_source,
    profile_centering_mode = profile_centering_mode
  )
  path_features$frame_provenance$lsf <- lsf_provenance
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
  path_seg$feature_axis_provenance <- path_seg$wavelength_provenance
  path_seg$wavelength_provenance <- path_features$wavelength_provenance
  path_seg$kinematic_provenance <- path_features$frame_provenance
  if (nrow(path_features$table)) {
    path_features$table$path_signature_segment <- path_seg$cluster_map[cbind(path_features$table$x, path_features$table$y)]
  }
}

message("Saving products...")
plot_cont(star$collapsed, file.path(out_dir, paste0(object_prefix, "_white_light.png")), "white light")
plot_cont(ifelse(support_raw, 1, NA_real_), file.path(out_dir, paste0(object_prefix, "_starlet_support_raw.png")), sprintf("raw starlet support: %d spaxels", sum(support_raw)), palette = c("#F2D06B", "#F2D06B"), limits = c(0, 1))
plot_cont(ifelse(support_starlet, 1, NA_real_), file.path(out_dir, paste0(object_prefix, "_starlet_support.png")), sprintf("connected starlet support: %d spaxels", sum(support_starlet)), palette = c("#F2D06B", "#F2D06B"), limits = c(0, 1))
if (identical(kinematic_support, "line_flux")) {
  plot_cont(ifelse(support, 1, NA_real_), file.path(out_dir, paste0(file_prefix, "_line_flux_support.png")), sprintf("%s flux support (%.1f sigma): %d spaxels", line$name, line_flux_sigma, sum(support)), palette = c("#F2D06B", "#F2D06B"), limits = c(0, 1))
}
if (run_spectral_segmentation) {
  plot_seg(seg$cluster_map, file.path(out_dir, sprintf("%s_capivara_segments_n%d.png", object_prefix, ncomp)), sprintf("Capivara full-spectrum segments (N=%d)", ncomp), n = ncomp)
}
plot_cont(log10(kin$flux), file.path(out_dir, paste0(file_prefix, "_flux_log.png")), paste(line$name, "log flux"), limits = robust_limits(log10(kin$flux)))
plot_cont(kin$velocity, file.path(out_dir, paste0(file_prefix, "_velocity_systemic.png")), paste(line$name, "systemic-relative velocity"), palette = div_palette, limits = robust_limits(kin$velocity, symmetric = TRUE), midpoint = 0)
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
tab$starlet_support <- as.vector(support_starlet)
tab$kinematic_support <- as.vector(support)
tab$capivara_segment <- as.vector(seg$cluster_map)
tab[[paste0(line$slug, "_flux")]] <- as.vector(kin$flux)
tab[[paste0(line$slug, "_velocity_systemic")]] <- as.vector(kin$velocity)
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
  starlet_support = support_starlet + 0,
  kinematic_support = support + 0,
  capivara_segment = seg$cluster_map,
  setNames(list(log10(kin$flux)), paste0(line$slug, "_log_flux"))[[1]],
  setNames(list(kin$velocity), paste0(line$slug, "_velocity_systemic"))[[1]],
  setNames(list(kin$sigma), paste0(line$slug, "_sigma"))[[1]],
  setNames(list(kin$asymmetry), paste0(line$slug, "_asymmetry"))[[1]],
  setNames(list(kin$h3_proxy), paste0(line$slug, "_h3_proxy"))[[1]],
  setNames(list(kin$h4_proxy), paste0(line$slug, "_h4_proxy"))[[1]],
  kinematic_aware_segment = kin_seg$cluster_map
)
names(maps)[4:9] <- c(
  paste0(line$slug, "_log_flux"),
  paste0(line$slug, "_velocity_systemic"),
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
stopifnot(length(names(maps)) == length(maps), !anyNA(names(maps)),
          all(nzchar(names(maps))), !anyDuplicated(names(maps)))
FITSio::writeFITSim(map_stack, file.path(out_dir, paste0(file_prefix, "_capivara_kinematic_maps.fits")), type = "double")
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
    frame_provenance = kin$frame_provenance,
    wavelength_provenance = kin$wavelength_provenance,
    profile_centering_mode = profile_centering_mode,
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
    support_starlet = support_starlet,
    support_method = kinematic_support,
    support_info = support_info,
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
