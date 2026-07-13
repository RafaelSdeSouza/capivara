.capivara_workflow_file <- function(name, repo_root = NULL) {
  installed <- system.file("extdata", "kinematics", name, package = "capivara")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }

  # This fallback is only for running directly from a source checkout.
  if (!is.null(repo_root) && length(repo_root) && !is.na(repo_root[[1]]) && nzchar(repo_root[[1]])) {
    candidate <- file.path(
      normalizePath(repo_root[[1]], mustWork = TRUE),
      "inst", "extdata", "kinematics", name
    )
    if (file.exists(candidate)) {
      return(candidate)
    }
  }

  stop(
    "Could not locate the bundled Capivara kinematics workflow. Reinstall capivara.",
    call. = FALSE
  )
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

#' Run the native MaNGA kinematics and bisymmetric-bar workflow
#'
#' This is the Capivara 2.0 entry point for native emission-line kinematics.
#' It reads a MaNGA LOGCUBE, resolves the redshift locally when needed, creates
#' standard, kinematic-aware, and path-signature segmentations, and fits a
#' local NIRVANA-style bisymmetric velocity model.
#'
#' @param cube_path Path to a MaNGA LOGCUBE.
#' @param redshift Optional redshift. If `NA`, uses FITS header then local
#'   DRPall metadata.
#' @param emission_line Emission line alias, currently typically `"halpha"`.
#' @param segmentation_mode One of `"kinematic"`, `"path_signature"`,
#'   `"spectral"`, or `"all"`.
#' @param output_dir Directory for products.
#' @param object_id Optional object identifier. Defaults to inferred plate-IFU.
#' @param repo_root Optional Capivara repository root used only when running
#'   directly from a development checkout. Installed packages do not need it.
#' @param knn_k Nearest-neighbour graph size used by the native Capivara
#'   segmentation runners.
#' @param n_segments Number of standard kinematic/spectral segments.
#' @param n_path_segments Number of path-signature-aware segments.
#' @param starlet_scales Starlet scales used for the support mask.
#' @param include_coarse_starlet If `TRUE`, include the coarse starlet scale in
#'   the support mask.
#' @param display_orientation Plot display orientation passed to the bar plots.
#' @param disc_pa_deg Optional fixed disc position angle in degrees.
#' @param disc_inc_deg Optional fixed disc inclination in degrees.
#' @param bar_phi_deg Optional fixed bar angle in degrees.
#' @param use_bar_support_mask If `TRUE`, restrict the bar term support around
#'   the bar angle.
#' @param bar_support_width_deg Angular half-width for the optional bar support.
#' @param robust_fit If `TRUE`, use robust reweighting in the model fit.
#' @param smooth_lambda First-order smoothing penalty for velocity profiles.
#' @param second_order_lambda Curvature penalty for second-order bar terms.
#' @param max_v2_fraction Maximum allowed mean second-order amplitude relative
#'   to the circular term before shrinkage.
#' @param max_mean_v2 Absolute cap on mean second-order velocity amplitude.
#' @param show_plots If `TRUE`, print the native, model, and component plots.
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
  native_runner <- .capivara_workflow_file("native_kinematics_workflow.R", repo_root)
  model_runner <- .capivara_workflow_file("native_bisymmetric_workflow.R", repo_root)

  native_env <- c(
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
  native_scope <- new.env(parent = environment(run_manga_bar_model))
  source(native_runner, local = native_scope)

  native_rds <- file.path(output_dir, paste0(prefix, "_", emission_line, "_capivara_kinematic_results.rds"))
  model_prefix <- paste0(prefix, "_", emission_line, "_bisymmetric")
  model_env <- c(
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
  model_scope <- new.env(parent = environment(run_manga_bar_model))
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

#' Run the Capivara kinematic analysis module
#'
#' Friendly alias for the native MaNGA workflow. `segmentation_mode = "all"`
#' creates traditional spectral, kinematic-aware, and path-signature
#' segmentations before fitting the disc and bisymmetric models.
#'
#' @inheritParams run_manga_bar_model
#' @param ... Additional geometry or robust-fitting controls accepted by
#'   [run_manga_bar_model()].
#' @return A `capivara_manga_bar_result` containing the maps, segmentations,
#'   fitted models, plots, and paths to saved products.
#' @export
run_kinematic_analysis <- function(cube_path,
                                   redshift = NA_real_,
                                   emission_line = "halpha",
                                   segmentation_mode = "all",
                                   output_dir = NULL,
                                   object_id = NULL,
                                   knn_k = 100,
                                   n_segments = 25,
                                   n_path_segments = 45,
                                   show_plots = interactive(),
                                   ...) {
  run_manga_bar_model(
    cube_path = cube_path,
    redshift = redshift,
    emission_line = emission_line,
    segmentation_mode = segmentation_mode,
    output_dir = output_dir,
    object_id = object_id,
    knn_k = knn_k,
    n_segments = n_segments,
    n_path_segments = n_path_segments,
    show_plots = show_plots,
    ...
  )
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
