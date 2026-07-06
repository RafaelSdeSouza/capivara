#!/usr/bin/env Rscript

# Capivara-native bisymmetric kinematic modelling for any MaNGA LOGCUBE.
# Either edit this block and run:
#   Rscript extensions/capivaraKinematics/examples/native_manga_bisymmetric.R
# or pass the main inputs:
#   Rscript extensions/capivaraKinematics/examples/native_manga_bisymmetric.R cube.fits 0.035 /tmp/out manga8333

cube_path <- "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/normal_bar/manga-8078-12703-LOGCUBE.fits"
redshift <- 0.0281
output_dir <- NA_character_
object_id <- "manga8078_12703"




emission_line <- "halpha"

n_segments <- 25
n_path_segments <- 45
knn_k <- 50
run_spectral_segmentation <- FALSE
run_path_signatures <- FALSE

starlet_scales <- "2:5"
include_coarse_starlet <- TRUE

line_window_kms <- 600
peak_search_kms <- 350
centroid_window_kms <- 220

# Supply geometry when known. If `disc_pa_deg` is NA, it is estimated from the
# native velocity gradient. If `disc_inc_deg` is NA, 60 deg is used as an
# exploratory placeholder.
disc_pa_deg <- NA_real_
disc_inc_deg <- NA_real_
bar_phi_deg <- NA_real_
display_orientation <- "rot90_cw"

robust_fit <- TRUE
smooth_lambda <- 10
second_order_lambda <- 25
max_v2_fraction <- 0.8
max_mean_v2 <- 350

# No edits are usually needed below this line. The code below runs in the same
# R session, so RStudio keeps the intermediate objects available for inspection:
# `star`, `support`, `seg`, `kin`, `kin_seg`, `result`, and `diagnostics`.
# Set `run_spectral_segmentation <- TRUE` for the slower full-spectrum
# Capivara segmentation panel; quick model checks do not need it.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L) cube_path <- args[[1]]
if (length(args) >= 2L) redshift <- as.numeric(args[[2]])
if (length(args) >= 3L) output_dir <- args[[3]]
if (length(args) >= 4L) object_id <- args[[4]]
if (length(args) >= 5L) emission_line <- args[[5]]

script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "extensions/capivaraKinematics/examples/native_manga_bisymmetric.R")
}

is_capivara_repo <- function(path) {
  if (!dir.exists(path)) {
    return(FALSE)
  }
  file.exists(file.path(path, "DESCRIPTION")) &&
    file.exists(file.path(path, "scripts", "run_manga_10218_bar_merger_maps.R"))
}

find_capivara_repo <- function(script_path) {
  candidates <- unique(c(
    Sys.getenv("CAPIVARA_REPO_ROOT", unset = NA_character_),
    dirname(script_path),
    getwd(),
    "/Users/rd23aag/Documents/GitHub/capivara"
  ))

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(active_path)) {
      candidates <- unique(c(dirname(active_path), candidates))
    }
  }

  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (candidate in candidates) {
    here <- normalizePath(candidate, mustWork = FALSE)
    repeat {
      if (is_capivara_repo(here)) {
        return(normalizePath(here, mustWork = TRUE))
      }
      parent <- dirname(here)
      if (identical(parent, here)) {
        break
      }
      here <- parent
    }
  }

  stop(
    "Could not find the Capivara repo. In RStudio, open/run this file from ",
    "/Users/rd23aag/Documents/GitHub/capivara or set CAPIVARA_REPO_ROOT.",
    call. = FALSE
  )
}

repo_root <- find_capivara_repo(script_path)

native_runner <- file.path(repo_root, "scripts", "run_manga_10218_bar_merger_maps.R")
model_runner <- file.path(repo_root, "extensions", "capivaraKinematics", "scripts", "run_8078_native_bisymmetric_model.R")

if (!file.exists(cube_path)) {
  stop("Set `cube_path` at the top of this script or pass it as the first command-line argument.", call. = FALSE)
}
if (is.na(output_dir) || !nzchar(output_dir)) {
  output_dir <- file.path(dirname(cube_path), "capivara_outputs", object_id)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(object_id))

native_settings <- c(
  CAPIVARA_REPO_ROOT = repo_root,
  CAPIVARA_USE_ENV_INPUTS = "true",
  CAPIVARA_CUBE_PATH = cube_path,
  CAPIVARA_OUTPUT_DIR = output_dir,
  CAPIVARA_REDSHIFT = as.character(redshift),
  CAPIVARA_LINE = emission_line,
  CAPIVARA_OUTPUT_PREFIX = prefix,
  CAPIVARA_KNN = as.character(knn_k),
  CAPIVARA_NCOMP = as.character(n_segments),
  CAPIVARA_RUN_SPECTRAL_SEGMENTATION = if (run_spectral_segmentation) "true" else "false",
  CAPIVARA_RUN_PATH_SIGNATURES = if (run_path_signatures) "true" else "false",
  CAPIVARA_PATH_KNN = as.character(knn_k),
  CAPIVARA_PATH_NCOMP = as.character(n_path_segments),
  CAPIVARA_PATH_SPATIAL_WEIGHT = "0.10",
  CAPIVARA_STARLET_SCALES = starlet_scales,
  CAPIVARA_STARLET_INCLUDE_COARSE = if (include_coarse_starlet) "true" else "false",
  CAPIVARA_LINE_WINDOW_KMS = as.character(line_window_kms),
  CAPIVARA_PEAK_SEARCH_KMS = as.character(peak_search_kms),
  CAPIVARA_CENTROID_WINDOW_KMS = as.character(centroid_window_kms)
)
restore_env <- function(keys, old_values) {
  for (i in seq_along(keys)) {
    if (is.na(old_values[[i]])) {
      Sys.unsetenv(keys[[i]])
    } else {
      do.call(Sys.setenv, stats::setNames(list(old_values[[i]]), keys[[i]]))
    }
  }
}
old_native_env <- Sys.getenv(names(native_settings), unset = NA_character_)
on.exit(restore_env(names(native_settings), old_native_env), add = TRUE)
do.call(Sys.setenv, as.list(native_settings))

source(native_runner, local = environment())
if (interactive() && exists("panel", inherits = FALSE)) {
  print(panel)
}

native_rds <- file.path(output_dir, paste0(prefix, "_", emission_line, "_capivara_kinematic_results.rds"))
model_prefix <- paste0(prefix, "_", emission_line, "_bisymmetric")

model_settings <- c(
  CAPIVARA_REPO_ROOT = repo_root,
  CAPIVARA_MODEL_USE_ENV_INPUTS = "true",
  CAPIVARA_MODEL_NATIVE_RDS = native_rds,
  CAPIVARA_MODEL_OUTPUT_DIR = output_dir,
  CAPIVARA_MODEL_OUTPUT_PREFIX = model_prefix,
  CAPIVARA_MODEL_PLATEIFU = object_id,
  CAPIVARA_MODEL_DISPLAY_ORIENTATION = display_orientation,
  CAPIVARA_MODEL_ROBUST = if (robust_fit) "true" else "false",
  CAPIVARA_MODEL_SMOOTH_LAMBDA = as.character(smooth_lambda),
  CAPIVARA_MODEL_SECOND_ORDER_LAMBDA = as.character(second_order_lambda),
  CAPIVARA_MODEL_MAX_V2_FRACTION = as.character(max_v2_fraction),
  CAPIVARA_MODEL_MAX_MEAN_V2 = as.character(max_mean_v2)
)
if (is.finite(disc_pa_deg)) {
  model_settings <- c(model_settings, CAPIVARA_MODEL_PA_DEG = as.character(disc_pa_deg))
}
if (is.finite(disc_inc_deg)) {
  model_settings <- c(model_settings, CAPIVARA_MODEL_INC_DEG = as.character(disc_inc_deg))
}
if (is.finite(bar_phi_deg)) {
  model_settings <- c(model_settings, CAPIVARA_MODEL_BAR_PHI_DEG = as.character(bar_phi_deg))
}
old_model_env <- Sys.getenv(names(model_settings), unset = NA_character_)
on.exit(restore_env(names(model_settings), old_model_env), add = TRUE)
do.call(Sys.setenv, as.list(model_settings))

source(model_runner, local = environment())

manifest <- data.frame(
  product = c(
    "native_results",
    "model_results",
    "model_figure",
    "component_figure",
    "native_panel",
    "spaxel_table"
  ),
  path = c(
    native_rds,
    file.path(output_dir, paste0(model_prefix, ".rds")),
    file.path(output_dir, paste0(model_prefix, "_model.png")),
    file.path(output_dir, paste0(model_prefix, "_components.png")),
    file.path(output_dir, paste0(prefix, "_", emission_line, "_capivara_kinematic_panel.png")),
    file.path(output_dir, paste0(prefix, "_", emission_line, "_capivara_kinematic_spaxel_table.csv"))
  )
)
utils::write.csv(manifest, file.path(output_dir, paste0(prefix, "_outputs_manifest.csv")), row.names = FALSE)

if (interactive()) {
  print(manifest)
}

invisible(list(
  manifest = manifest,
  native = if (exists("kin", inherits = FALSE)) list(star = star, support = support, seg = seg, kin = kin, kin_seg = kin_seg) else NULL,
  model = if (exists("result", inherits = FALSE)) result else NULL
))
