#!/usr/bin/env Rscript

# Capivara-native bisymmetric kinematic model for one MaNGA LOGCUBE.
# Edit the small block below and run this file. It creates all images and RDS
# products in `output_dir` without requiring DAP products.

cube_path <- "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/normal_bar/manga-8078-12703-LOGCUBE.fits"
redshift <- 0.0281
output_dir <- "/private/tmp/capivara_native_8078_example"

plateifu <- "8078-12703"
line <- "halpha"
n_segments <- 25
n_path_segments <- 45
knn_k <- 100

starlet_scales <- "1:5"
include_coarse_starlet <- TRUE

line_window_kms <- 600
peak_search_kms <- 350
centroid_window_kms <- 220

disc_pa_deg <- 14.7
disc_inc_deg <- 36.4
bar_phi_deg <- 40.9

# No edits are usually needed below this line.

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "extensions/capivaraKinematics/examples/native_manga_bisymmetric_8078.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)

native_runner <- file.path(repo_root, "scripts", "run_manga_10218_bar_merger_maps.R")
model_runner <- file.path(repo_root, "extensions", "capivaraKinematics", "scripts", "run_8078_native_bisymmetric_model.R")

stopifnot(file.exists(cube_path))
stopifnot(file.exists(native_runner))
stopifnot(file.exists(model_runner))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- paste0("manga", gsub("-", "_", plateifu), "_native")
native_env <- c(
  CAPIVARA_REDSHIFT = as.character(redshift),
  CAPIVARA_LINE = line,
  CAPIVARA_OUTPUT_PREFIX = prefix,
  CAPIVARA_KNN = as.character(knn_k),
  CAPIVARA_NCOMP = as.character(n_segments),
  CAPIVARA_PATH_KNN = as.character(knn_k),
  CAPIVARA_PATH_NCOMP = as.character(n_path_segments),
  CAPIVARA_PATH_SPATIAL_WEIGHT = "0.10",
  CAPIVARA_STARLET_SCALES = starlet_scales,
  CAPIVARA_STARLET_INCLUDE_COARSE = if (include_coarse_starlet) "true" else "false",
  CAPIVARA_LINE_WINDOW_KMS = as.character(line_window_kms),
  CAPIVARA_PEAK_SEARCH_KMS = as.character(peak_search_kms),
  CAPIVARA_CENTROID_WINDOW_KMS = as.character(centroid_window_kms)
)

status <- system2(
  "Rscript",
  c(native_runner, cube_path, output_dir),
  env = paste(names(native_env), native_env, sep = "="),
  stdout = file.path(output_dir, paste0(prefix, "_native_capivara.log")),
  stderr = file.path(output_dir, paste0(prefix, "_native_capivara.log"))
)
if (!identical(status, 0L)) {
  stop("Native Capivara kinematic extraction failed.", call. = FALSE)
}

native_rds <- file.path(output_dir, paste0(prefix, "_", line, "_capivara_kinematic_results.rds"))
model_prefix <- paste0(prefix, "_", line, "_bisymmetric")

model_env <- c(
  CAPIVARA_MODEL_PLATEIFU = plateifu,
  CAPIVARA_MODEL_PA_DEG = as.character(disc_pa_deg),
  CAPIVARA_MODEL_INC_DEG = as.character(disc_inc_deg),
  CAPIVARA_MODEL_BAR_PHI_DEG = as.character(bar_phi_deg)
)
status <- system2(
  "Rscript",
  c(model_runner, native_rds, output_dir, model_prefix),
  env = paste(names(model_env), model_env, sep = "="),
  stdout = file.path(output_dir, paste0(prefix, "_bisymmetric_model.log")),
  stderr = file.path(output_dir, paste0(prefix, "_bisymmetric_model.log"))
)
if (!identical(status, 0L)) {
  stop("Native bisymmetric model failed.", call. = FALSE)
}

manifest <- data.frame(
  product = c("native_results", "model_results", "model_figure", "component_figure"),
  path = c(
    native_rds,
    file.path(output_dir, paste0(model_prefix, ".rds")),
    file.path(output_dir, paste0(model_prefix, "_model.png")),
    file.path(output_dir, paste0(model_prefix, "_components.png"))
  )
)
utils::write.csv(manifest, file.path(output_dir, paste0(prefix, "_outputs_manifest.csv")), row.names = FALSE)

invisible(manifest)
