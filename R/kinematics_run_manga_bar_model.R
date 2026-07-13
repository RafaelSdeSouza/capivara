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

.capivara_run_native_kinematics <- function(cube_path,
                                             redshift,
                                             emission_line,
                                             output_dir,
                                             object_id,
                                             repo_root,
                                             knn_k,
                                             n_segments,
                                             n_path_segments,
                                             run_spectral_segmentation,
                                             run_path_signatures,
                                             starlet_scales,
                                             include_coarse_starlet) {
  prefix <- gsub("[^A-Za-z0-9]+", "_", tolower(object_id))
  native_runner <- .capivara_workflow_file("native_kinematics_workflow.R", repo_root)
  native_env <- c(
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
    CAPIVARA_STARLET_SCALES = starlet_scales,
    CAPIVARA_STARLET_INCLUDE_COARSE = if (include_coarse_starlet) "true" else "false",
    CAPIVARA_LINE_WINDOW_KMS = "600",
    CAPIVARA_PEAK_SEARCH_KMS = "350",
    CAPIVARA_CENTROID_WINDOW_KMS = "220"
  )
  old_native_env <- .capivara_set_env(native_env)
  on.exit(.capivara_restore_env(names(native_env), old_native_env), add = TRUE)

  native_scope <- new.env(parent = environment(.capivara_run_native_kinematics))
  source(native_runner, local = native_scope)

  list(
    prefix = prefix,
    native_rds = file.path(
      output_dir,
      paste0(prefix, "_", emission_line, "_capivara_kinematic_results.rds")
    ),
    native = readRDS(file.path(
      output_dir,
      paste0(prefix, "_", emission_line, "_capivara_kinematic_results.rds")
    )),
    panel = native_scope$panel
  )
}

#' Segment an IFU cube using emission-line kinematics
#'
#' Builds native emission-line maps, then clusters either the kinematic feature
#' cube or the path-signature feature cube. It deliberately does not run the
#' ordinary full-spectrum Capivara segmentation: use [segment()] or
#' [segment_large()] when spectral segmentation is the scientific objective.
#'
#' @param cube_path Path to a MaNGA LOGCUBE.
#' @param redshift Optional redshift. If `NA`, uses FITS header then local
#'   MaNGA metadata.
#' @param emission_line Emission line alias, such as `"halpha"` or
#'   `"oiii5007"`.
#' @param segmentation_mode `"kinematic"` clusters flux, velocity, dispersion,
#'   and profile-shape maps. `"path_signature"` additionally computes the
#'   path-signature feature segmentation.
#' @param output_dir Directory for saved products. Defaults beside the cube.
#' @param object_id Optional output identifier.
#' @param knn_k kNN graph size for sparse Ward clustering.
#' @param n_segments Number of kinematic-aware segments.
#' @param n_path_segments Number of path-signature segments.
#' @param starlet_scales Starlet support scales.
#' @param include_coarse_starlet Include the coarse starlet plane in support.
#' @param show_plots Print the compact kinematic panel.
#' @return A `capivara_kinematic_segmentation` object containing native maps,
#'   support, kinematic segmentation, optional path segmentation, and paths.
#' @export
segment_kinematics <- function(cube_path,
                               redshift = NA_real_,
                               emission_line = "halpha",
                               segmentation_mode = c("kinematic", "path_signature"),
                               output_dir = NULL,
                               object_id = NULL,
                               knn_k = 50,
                               n_segments = 25,
                               n_path_segments = 45,
                               starlet_scales = "2:5",
                               include_coarse_starlet = FALSE,
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

  native_run <- .capivara_run_native_kinematics(
    cube_path = cube_path,
    redshift = z_info$redshift,
    emission_line = emission_line,
    output_dir = output_dir,
    object_id = object_id,
    repo_root = NULL,
    knn_k = knn_k,
    n_segments = n_segments,
    n_path_segments = n_path_segments,
    run_spectral_segmentation = FALSE,
    run_path_signatures = identical(segmentation_mode, "path_signature"),
    starlet_scales = starlet_scales,
    include_coarse_starlet = include_coarse_starlet
  )

  out <- list(
    object_id = object_id,
    redshift = z_info$redshift,
    redshift_source = z_info$source,
    output_dir = output_dir,
    native_rds = native_run$native_rds,
    native = native_run$native,
    panel = native_run$panel,
    segmentation_mode = segmentation_mode,
    segmentation = if (identical(segmentation_mode, "path_signature")) {
      native_run$native$path_signature
    } else {
      native_run$native$kinematic_aware
    }
  )
  class(out) <- c("capivara_kinematic_segmentation", class(out))
  if (isTRUE(show_plots)) {
    print(out$panel)
  }
  out
}

.capivara_model_control <- function(control = list()) {
  defaults <- list(
    starlet_scales = "2:5",
    include_coarse_starlet = FALSE,
    display_orientation = "rot90_cw",
    disc_pa_deg = NA_real_,
    disc_inc_deg = NA_real_,
    bar_phi_deg = NA_real_,
    use_bar_support_mask = FALSE,
    bar_support_width_deg = 25,
    robust_fit = TRUE,
    smooth_lambda = 10,
    second_order_lambda = 25,
    max_v2_fraction = 0.8,
    max_mean_v2 = 350
  )
  if (is.null(control)) {
    return(defaults)
  }
  if (!is.list(control) || is.null(names(control))) {
    stop("`model_control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown `model_control` entries: ", paste(unknown, collapse = ", "), call. = FALSE)
  }
  utils::modifyList(defaults, control)
}

.capivara_run_kinematic_model <- function(cube_path,
                                           redshift,
                                           emission_line,
                                           segmentation_mode,
                                           model,
                                           output_dir,
                                           object_id,
                                           knn_k,
                                           n_segments,
                                           n_path_segments,
                                           model_control,
                                           show_plots,
                                           repo_root = NULL) {
  segmentation_mode <- match.arg(segmentation_mode, c("kinematic", "path_signature"))
  model <- .capivara_match_kinematic_model(model)
  control <- .capivara_model_control(model_control)
  if (.capivara_is_bar_model(model) && !is.finite(control$bar_phi_deg)) {
    stop(
      "`model = \"bisymmetric_bar\"` requires `model_control = list(bar_phi_deg = ...)`. ",
      "The bar angle is a physical prior and is never assumed from the disc PA.",
      call. = FALSE
    )
  }

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

  native_run <- .capivara_run_native_kinematics(
    cube_path = cube_path,
    redshift = z_info$redshift,
    emission_line = emission_line,
    output_dir = output_dir,
    object_id = object_id,
    repo_root = repo_root,
    knn_k = knn_k,
    n_segments = n_segments,
    n_path_segments = n_path_segments,
    run_spectral_segmentation = FALSE,
    run_path_signatures = identical(segmentation_mode, "path_signature"),
    starlet_scales = control$starlet_scales,
    include_coarse_starlet = control$include_coarse_starlet
  )
  model_prefix <- paste(native_run$prefix, emission_line, model, sep = "_")
  model_runner <- .capivara_workflow_file("native_bisymmetric_workflow.R", repo_root)
  model_env <- c(
    CAPIVARA_MODEL_USE_ENV_INPUTS = "true",
    CAPIVARA_MODEL_KIND = model,
    CAPIVARA_MODEL_NATIVE_RDS = native_run$native_rds,
    CAPIVARA_MODEL_OUTPUT_DIR = output_dir,
    CAPIVARA_MODEL_OUTPUT_PREFIX = model_prefix,
    CAPIVARA_MODEL_PLATEIFU = object_id,
    CAPIVARA_MODEL_DISPLAY_ORIENTATION = control$display_orientation,
    CAPIVARA_MODEL_ROBUST = if (isTRUE(control$robust_fit)) "true" else "false",
    CAPIVARA_MODEL_SMOOTH_LAMBDA = as.character(control$smooth_lambda),
    CAPIVARA_MODEL_SECOND_ORDER_LAMBDA = as.character(control$second_order_lambda),
    CAPIVARA_MODEL_MAX_V2_FRACTION = as.character(control$max_v2_fraction),
    CAPIVARA_MODEL_MAX_MEAN_V2 = as.character(control$max_mean_v2),
    CAPIVARA_MODEL_USE_BAR_MASK = if (isTRUE(control$use_bar_support_mask)) "true" else "false",
    CAPIVARA_MODEL_BAR_MASK_WIDTH_DEG = as.character(control$bar_support_width_deg)
  )
  if (is.finite(control$disc_pa_deg)) {
    model_env <- c(model_env, CAPIVARA_MODEL_PA_DEG = as.character(control$disc_pa_deg))
  }
  if (is.finite(control$disc_inc_deg)) {
    model_env <- c(model_env, CAPIVARA_MODEL_INC_DEG = as.character(control$disc_inc_deg))
  }
  if (.capivara_is_bar_model(model)) {
    model_env <- c(model_env, CAPIVARA_MODEL_BAR_PHI_DEG = as.character(control$bar_phi_deg))
  }
  old_model_env <- .capivara_set_env(model_env)
  on.exit(.capivara_restore_env(names(model_env), old_model_env), add = TRUE)
  model_scope <- new.env(parent = environment(.capivara_run_kinematic_model))
  source(model_runner, local = model_scope)

  result <- list(
    object_id = object_id,
    redshift = z_info$redshift,
    redshift_source = z_info$source,
    model = model,
    model_control = control,
    output_dir = output_dir,
    native_rds = native_run$native_rds,
    model_rds = file.path(output_dir, paste0(model_prefix, ".rds")),
    model_png = file.path(output_dir, paste0(model_prefix, "_model.png")),
    components_png = if (.capivara_is_bar_model(model)) {
      file.path(output_dir, paste0(model_prefix, "_components.png"))
    } else {
      NA_character_
    },
    native = native_run$native,
    model_result = readRDS(file.path(output_dir, paste0(model_prefix, ".rds"))),
    panel = native_run$panel,
    model_plot = model_scope$model_plot,
    component_plot = if (exists("component_plot", model_scope, inherits = FALSE)) {
      model_scope$component_plot
    } else {
      NULL
    }
  )
  class(result) <- c("capivara_kinematic_result", "capivara_manga_bar_result", class(result))
  if (isTRUE(show_plots)) {
    print(result$panel)
    print(result$model_plot)
    if (!is.null(result$component_plot)) {
      print(result$component_plot)
    }
  }
  result
}

#' Run a Capivara kinematic analysis
#'
#' This is the compact entry point for native emission-line kinematics. It
#' makes a kinematic-aware segmentation first, then fits the requested model
#' module. The default model is an axisymmetric disc because it is applicable
#' to both barred and unbarred galaxies. Full-spectrum segmentation is a
#' separate analysis and is run with [segment()] or [segment_large()].
#'
#' @param cube_path Path to an IFU FITS cube.
#' @param redshift Redshift. Leave as `NA` for a MaNGA cube with local metadata.
#' @param emission_line Emission-line alias, such as `"halpha"` or `"oiii5007"`.
#' @param segmentation_mode `"kinematic"` for line-map features or
#'   `"path_signature"` to also compute path-signature regions.
#' @param model Kinematic model module. Use `"axisymmetric"` by default;
#'   `"bisymmetric_bar"` is an explicit bar hypothesis.
#' @param output_dir Directory for saved products. Defaults beside the cube.
#' @param object_id Optional output identifier.
#' @param knn_k kNN graph size for kinematic clustering.
#' @param n_segments Number of kinematic-aware segments.
#' @param n_path_segments Number of path-signature segments.
#' @param model_control Named list of model controls. For a bar model it must
#'   contain `bar_phi_deg`. See [kinematic_models()] for available modules.
#' @param show_plots Print figures as they are generated.
#' @return A `capivara_kinematic_result` with maps, segmentation, model, plots,
#'   and saved-product paths.
#' @export
run_kinematic_analysis <- function(cube_path,
                                   redshift = NA_real_,
                                   emission_line = "halpha",
                                   segmentation_mode = c("kinematic", "path_signature"),
                                   model = c("axisymmetric", "bisymmetric_bar"),
                                   output_dir = NULL,
                                   object_id = NULL,
                                   knn_k = 50,
                                   n_segments = 25,
                                   n_path_segments = 45,
                                   model_control = list(),
                                   show_plots = interactive()) {
  segmentation_mode <- match.arg(segmentation_mode)
  model <- .capivara_match_kinematic_model(model)
  .capivara_run_kinematic_model(
    cube_path = cube_path,
    redshift = redshift,
    emission_line = emission_line,
    segmentation_mode = segmentation_mode,
    model = model,
    output_dir = output_dir,
    object_id = object_id,
    knn_k = knn_k,
    n_segments = n_segments,
    n_path_segments = n_path_segments,
    model_control = model_control,
    show_plots = show_plots
  )
}

#' Run the explicit bisymmetric-bar module for a MaNGA cube
#'
#' This convenience wrapper is intentionally bar-specific. It requires a bar
#' angle measured from imaging or supplied by the user; it never substitutes
#' the disc position angle for a bar angle. For an ordinary rotation model, use
#' [run_kinematic_analysis()] with its axisymmetric default.
#'
#' @inheritParams run_kinematic_analysis
#' @param bar_phi_deg In-plane bar angle in degrees relative to the disc major
#'   axis.
#' @return A `capivara_kinematic_result` with the bisymmetric model and its
#'   component decomposition.
#' @export
run_manga_bar_model <- function(cube_path,
                                redshift = NA_real_,
                                emission_line = "halpha",
                                segmentation_mode = c("kinematic", "path_signature"),
                                bar_phi_deg,
                                output_dir = NULL,
                                object_id = NULL,
                                knn_k = 50,
                                n_segments = 25,
                                n_path_segments = 45,
                                model_control = list(),
                                show_plots = interactive()) {
  if (missing(bar_phi_deg) || !is.finite(bar_phi_deg)) {
    stop("Supply `bar_phi_deg` to run the bisymmetric bar module.", call. = FALSE)
  }
  if (is.null(model_control)) {
    model_control <- list()
  }
  if (!is.list(model_control)) {
    stop("`model_control` must be a named list.", call. = FALSE)
  }
  model_control$bar_phi_deg <- bar_phi_deg
  run_kinematic_analysis(
    cube_path = cube_path,
    redshift = redshift,
    emission_line = emission_line,
    segmentation_mode = segmentation_mode,
    model = "bisymmetric_bar",
    output_dir = output_dir,
    object_id = object_id,
    knn_k = knn_k,
    n_segments = n_segments,
    n_path_segments = n_path_segments,
    model_control = model_control,
    show_plots = show_plots
  )
}

#' @export
print.capivara_kinematic_result <- function(x, ...) {
  cat("Capivara kinematic result\n")
  cat("  object: ", x$object_id, "\n", sep = "")
  cat("  model: ", x$model, "\n", sep = "")
  cat("  redshift: ", x$redshift, " (", x$redshift_source, ")\n", sep = "")
  cat("  output_dir: ", x$output_dir, "\n", sep = "")
  cat("  model_png: ", x$model_png, "\n", sep = "")
  if (!is.na(x$components_png)) {
    cat("  components_png: ", x$components_png, "\n", sep = "")
  }
  invisible(x)
}

#' @export
plot.capivara_kinematic_result <- function(x,
                                           which = c("all", "segmentation", "model", "components"),
                                           ...) {
  which <- match.arg(which)
  if (which %in% c("all", "segmentation")) {
    print(x$panel)
  }
  if (which %in% c("all", "model")) {
    print(x$model_plot)
  }
  if (identical(which, "components")) {
    if (is.null(x$component_plot)) {
      stop("Component decomposition is available only for `model = \"bisymmetric_bar\"`.", call. = FALSE)
    }
    print(x$component_plot)
  } else if (identical(which, "all") && !is.null(x$component_plot)) {
    print(x$component_plot)
  }
  invisible(x)
}
