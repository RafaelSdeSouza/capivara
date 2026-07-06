.capivara_find_repo_root <- function(repo_root = NULL) {
  candidates <- character()
  if (!is.null(repo_root) && length(repo_root) && !is.na(repo_root[[1]]) && nzchar(repo_root[[1]])) {
    candidates <- c(candidates, repo_root[[1]])
  }
  env_root <- Sys.getenv("CAPIVARA_REPO_ROOT", unset = "")
  if (nzchar(env_root)) {
    candidates <- c(candidates, env_root)
  }
  candidates <- c(candidates, getwd(), normalizePath(file.path(getwd(), ".."), mustWork = FALSE))
  for (path in candidates) {
    if (dir.exists(path) && file.exists(file.path(path, "DESCRIPTION")) && dir.exists(file.path(path, "R"))) {
      return(normalizePath(path, mustWork = TRUE))
    }
  }
  stop("Could not find the Capivara repository root. Pass repo_root explicitly.", call. = FALSE)
}

.capivara_restore_env <- function(keys, old_values) {
  for (i in seq_along(keys)) {
    if (is.na(old_values[[i]])) {
      Sys.unsetenv(keys[[i]])
    } else {
      do.call(Sys.setenv, stats::setNames(list(old_values[[i]]), keys[[i]]))
    }
  }
}

.capivara_set_env <- function(values) {
  old <- Sys.getenv(names(values), unset = NA_character_)
  do.call(Sys.setenv, as.list(values))
  old
}

#' Run the prototype native MaNGA bar-modelling workflow
#'
#' This is the friendly Capivara-2.0-style entry point for the current
#' prototype. It reads a MaNGA LOGCUBE, resolves the redshift locally when
#' needed, makes native Capivara kinematic products, and fits the bisymmetric
#' model.
#'
#' @param cube_path Path to a MaNGA LOGCUBE.
#' @param redshift Optional redshift. If `NA`, uses FITS header then local
#'   DRPall metadata.
#' @param emission_line Emission line alias, currently typically `"halpha"`.
#' @param segmentation_mode One of `"kinematic"`, `"path_signature"`,
#'   `"spectral"`, or `"all"`.
#' @param output_dir Directory for products.
#' @param object_id Optional object identifier. Defaults to inferred plate-IFU.
#' @param repo_root Capivara repository root for this prototype.
#' @return A list with plots, paths, native result, model result, and metadata.
#' @export
run_manga_bar_model <- function(cube_path,
                                redshift = NA_real_,
                                emission_line = "halpha",
                                segmentation_mode = c("kinematic", "path_signature", "spectral", "all"),
                                output_dir = NULL,
                                object_id = NULL,
                                repo_root = NULL,
                                knn_k = 50,
                                n_segments = 25,
                                n_path_segments = 45,
                                starlet_scales = "2:5",
                                include_coarse_starlet = FALSE,
                                display_orientation = "rot90_cw",
                                disc_pa_deg = NA_real_,
                                disc_inc_deg = NA_real_,
                                bar_phi_deg = NA_real_,
                                use_bar_support_mask = is.finite(bar_phi_deg),
                                bar_support_width_deg = 25,
                                robust_fit = TRUE,
                                smooth_lambda = 10,
                                second_order_lambda = 25,
                                max_v2_fraction = 0.8,
                                max_mean_v2 = 350,
                                show_plots = interactive()) {
  segmentation_mode <- match.arg(segmentation_mode)
  repo_root <- .capivara_find_repo_root(repo_root)
  cube_path <- normalizePath(cube_path, mustWork = TRUE)

  z_info <- resolve_manga_redshift(cube_path, redshift = redshift)
  object_id <- if (!is.null(object_id) && nzchar(object_id)) object_id else z_info$plateifu
  if (is.na(object_id) || !nzchar(object_id)) {
    object_id <- tools::file_path_sans_ext(basename(cube_path))
  }
  if (is.null(output_dir) || !nzchar(output_dir)) {
    output_dir <- file.path(dirname(cube_path), "capivara_outputs", object_id)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  run_spectral_segmentation <- segmentation_mode %in% c("spectral", "all")
  run_path_signatures <- segmentation_mode %in% c("path_signature", "all")
  prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(object_id))
  native_runner <- file.path(repo_root, "scripts", "run_manga_10218_bar_merger_maps.R")
  model_runner <- file.path(repo_root, "extensions", "capivaraKinematics", "scripts", "run_8078_native_bisymmetric_model.R")
  if (!file.exists(native_runner)) {
    stop("Could not find native runner: ", native_runner, call. = FALSE)
  }
  if (!file.exists(model_runner)) {
    stop("Could not find model runner: ", model_runner, call. = FALSE)
  }

  native_env <- c(
    CAPIVARA_REPO_ROOT = repo_root,
    CAPIVARA_USE_ENV_INPUTS = "true",
    CAPIVARA_CUBE_PATH = cube_path,
    CAPIVARA_OUTPUT_DIR = output_dir,
    CAPIVARA_REDSHIFT = as.character(z_info$redshift),
    CAPIVARA_LINE = emission_line,
    CAPIVARA_OUTPUT_PREFIX = prefix,
    CAPIVARA_KNN = as.character(knn_k),
    CAPIVARA_NCOMP = as.character(n_segments),
    CAPIVARA_RUN_SPECTRAL_SEGMENTATION = if (run_spectral_segmentation) "true" else "false",
    CAPIVARA_RUN_PATH_SIGNATURES = if (run_path_signatures) "true" else "false",
    CAPIVARA_PATH_KNN = as.character(knn_k),
    CAPIVARA_PATH_NCOMP = as.character(n_path_segments),
    CAPIVARA_STARLET_SCALES = starlet_scales,
    CAPIVARA_STARLET_INCLUDE_COARSE = if (include_coarse_starlet) "true" else "false",
    CAPIVARA_LINE_WINDOW_KMS = "600",
    CAPIVARA_PEAK_SEARCH_KMS = "350",
    CAPIVARA_CENTROID_WINDOW_KMS = "220"
  )
  old_native_env <- .capivara_set_env(native_env)
  on.exit(.capivara_restore_env(names(native_env), old_native_env), add = TRUE)
  native_scope <- new.env(parent = globalenv())
  source(native_runner, local = native_scope)

  native_rds <- file.path(output_dir, paste0(prefix, "_", emission_line, "_capivara_kinematic_results.rds"))
  model_prefix <- paste0(prefix, "_", emission_line, "_bisymmetric")
  model_env <- c(
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
    CAPIVARA_MODEL_MAX_MEAN_V2 = as.character(max_mean_v2),
    CAPIVARA_MODEL_USE_BAR_MASK = if (use_bar_support_mask) "true" else "false",
    CAPIVARA_MODEL_BAR_MASK_WIDTH_DEG = as.character(bar_support_width_deg)
  )
  if (is.finite(disc_pa_deg)) {
    model_env <- c(model_env, CAPIVARA_MODEL_PA_DEG = as.character(disc_pa_deg))
  }
  if (is.finite(disc_inc_deg)) {
    model_env <- c(model_env, CAPIVARA_MODEL_INC_DEG = as.character(disc_inc_deg))
  }
  if (is.finite(bar_phi_deg)) {
    model_env <- c(model_env, CAPIVARA_MODEL_BAR_PHI_DEG = as.character(bar_phi_deg))
  }
  old_model_env <- .capivara_set_env(model_env)
  on.exit(.capivara_restore_env(names(model_env), old_model_env), add = TRUE)
  model_scope <- new.env(parent = globalenv())
  source(model_runner, local = model_scope)

  model_rds <- file.path(output_dir, paste0(model_prefix, ".rds"))
  result <- list(
    object_id = object_id,
    redshift = z_info$redshift,
    redshift_source = z_info$source,
    output_dir = output_dir,
    native_rds = native_rds,
    model_rds = model_rds,
    model_png = file.path(output_dir, paste0(model_prefix, "_model.png")),
    components_png = file.path(output_dir, paste0(model_prefix, "_components.png")),
    native = readRDS(native_rds),
    model = readRDS(model_rds),
    panel = native_scope$panel,
    model_plot = model_scope$model_plot,
    component_plot = model_scope$component_plot
  )
  class(result) <- c("capivara_manga_bar_result", class(result))
  if (isTRUE(show_plots)) {
    print(result$panel)
    print(result$model_plot)
    print(result$component_plot)
  }
  result
}

#' @export
print.capivara_manga_bar_result <- function(x, ...) {
  cat("Capivara MaNGA bar result\n")
  cat("  object: ", x$object_id, "\n", sep = "")
  cat("  redshift: ", x$redshift, " (", x$redshift_source, ")\n", sep = "")
  cat("  output_dir: ", x$output_dir, "\n", sep = "")
  cat("  model_png: ", x$model_png, "\n", sep = "")
  cat("  components_png: ", x$components_png, "\n", sep = "")
  invisible(x)
}

#' @export
plot.capivara_manga_bar_result <- function(x, which = c("all", "summary", "model", "components"), ...) {
  which <- match.arg(which)
  if (which %in% c("all", "summary")) {
    print(x$panel)
  }
  if (which %in% c("all", "model")) {
    print(x$model_plot)
  }
  if (which %in% c("all", "components")) {
    print(x$component_plot)
  }
  invisible(x)
}
