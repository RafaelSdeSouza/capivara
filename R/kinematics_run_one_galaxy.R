.capivara_read_config <- function(config) {
  if (is.list(config)) {
    return(config)
  }
  .capivara_require("yaml")
  if (!file.exists(config)) {
    stop("Config file does not exist: ", config, call. = FALSE)
  }
  yaml::read_yaml(config)
}

.capivara_cfg <- function(config, name, default = NULL) {
  value <- config[[name]]
  if (is.null(value) || length(value) == 0L) default else value
}

.capivara_nested_cfg <- function(config, section, name, default = NULL) {
  block <- config[[section]]
  if (is.null(block)) {
    return(default)
  }
  value <- block[[name]]
  if (is.null(value) || length(value) == 0L) default else value
}

.capivara_build_spaxels <- function(manga, capivara) {
  velocity <- manga$velocity
  segmentation_map <- capivara$segmentation_map
  if (!identical(dim(velocity), dim(segmentation_map))) {
    stop(
      "Velocity map and Capivara segmentation map dimensions differ: ",
      paste(dim(velocity), collapse = "x"),
      " vs ",
      paste(dim(segmentation_map), collapse = "x"),
      call. = FALSE
    )
  }

  dims <- dim(velocity)
  grid <- expand.grid(x = seq_len(dims[2]), y = seq_len(dims[1]))
  idx <- cbind(grid$y, grid$x)
  grid$velocity <- velocity[idx]
  grid$velocity_error <- if (is.null(manga$velocity_error)) NA_real_ else manga$velocity_error[idx]
  grid$flux <- if (is.null(manga$flux)) NA_real_ else manga$flux[idx]
  grid$snr <- if (is.null(manga$snr)) NA_real_ else manga$snr[idx]
  grid$manga_mask <- if (is.null(manga$mask)) 0 else manga$mask[idx]
  grid$seg_label <- as.integer(segmentation_map[idx])

  tab <- capivara$segment_table
  hit <- match(grid$seg_label, tab$label)
  grid$seg_class <- tab$class[hit]
  grid$use_for_disc_fit <- tab$use_for_disc_fit[hit]
  grid$use_for_bar_diagnostics <- tab$use_for_bar_diagnostics[hit]
  grid$seg_class[is.na(grid$seg_class)] <- "unassigned"
  grid$use_for_disc_fit[is.na(grid$use_for_disc_fit)] <- FALSE
  grid$use_for_bar_diagnostics[is.na(grid$use_for_bar_diagnostics)] <- FALSE
  grid$has_capivara_segment <- is.finite(grid$seg_label) & grid$seg_label > 0L
  grid$valid <- is.finite(grid$velocity) &
    (is.na(grid$manga_mask) | grid$manga_mask == 0)
  grid
}

.capivara_apply_velocity_window <- function(spaxels, velocity_window = NULL) {
  if (is.null(velocity_window) || length(velocity_window) != 2L) {
    return(spaxels)
  }
  velocity_window <- sort(as.numeric(velocity_window))
  if (!all(is.finite(velocity_window))) {
    return(spaxels)
  }
  inside <- is.finite(spaxels$velocity) &
    spaxels$velocity >= velocity_window[1] &
    spaxels$velocity <= velocity_window[2]
  spaxels$valid <- spaxels$valid & inside
  spaxels
}

.capivara_write_log <- function(path, lines) {
  .capivara_dir_create(dirname(path))
  writeLines(as.character(lines), path)
}

.capivara_report_text <- function(summary) {
  paste0(
    "For galaxy ", summary$plateifu,
    ", the axisymmetric disc model was fitted using ", summary$N_disc,
    " Capivara-labelled clean disc spaxels. The bar segment contains ",
    summary$N_bar, " valid spaxels. The fitted arctangent model has v_sys = ",
    round(summary$vsys, 2), " km/s, v_max = ", round(summary$vmax, 2),
    " km/s, and R_t = ", round(summary$Rt, 2),
    " spaxels. The median absolute disc-subtracted residual inside the bar is A_bar = ",
    round(summary$A_bar, 2), " km/s. The normalized bar residual strength is Q_kin = ",
    round(summary$Q_kin, 2),
    ". Values above unity indicate stronger residuals inside the bar than in the clean disc."
  )
}

#' Run one Capivara-native MaNGA barred-galaxy kinematic analysis
#'
#' @param config_file YAML config path, or a config list.
#' @return An analysis result list.
#' @export
run_one_galaxy <- function(config_file) {
  .capivara_require(c("readr", "yaml", "ggplot2", "patchwork"))
  config <- .capivara_read_config(config_file)

  plateifu <- as.character(.capivara_cfg(config, "plateifu", "unknown"))
  output_dir <- .capivara_cfg(config, "output_dir", file.path("outputs", plateifu))
  dirs <- file.path(output_dir, c("figures", "tables", "models", "logs"))
  invisible(lapply(dirs, .capivara_dir_create))

  log_lines <- c(
    paste("plateifu:", plateifu),
    paste("manga_maps_file:", .capivara_cfg(config, "manga_maps_file")),
    paste("capivara_segmentation_file:", .capivara_cfg(config, "capivara_segmentation_file")),
    paste("capivara_segment_table_file:", .capivara_cfg(config, "capivara_segment_table_file"))
  )
  log_file <- file.path(output_dir, "logs", paste0(plateifu, "_run_log.txt"))

  manga <- read_manga_maps(
    .capivara_cfg(config, "manga_maps_file"),
    velocity_component = .capivara_cfg(config, "velocity_component", "stellar")
  )
  cap <- read_capivara_output(
    .capivara_cfg(config, "capivara_segmentation_file"),
    .capivara_cfg(config, "capivara_segment_table_file")
  )
  spaxels <- .capivara_build_spaxels(manga, cap)
  spaxels <- .capivara_apply_velocity_window(
    spaxels,
    velocity_window = .capivara_cfg(config, "velocity_window_kms", NULL)
  )

  geometry_config <- .capivara_cfg(config, "geometry", list())
  geometry <- estimate_disc_geometry(
    spaxels,
    geometry = geometry_config,
    allow_placeholder_inclination = isTRUE(.capivara_cfg(config, "allow_placeholder_inclination", FALSE)),
    placeholder_inc_deg = as.numeric(.capivara_cfg(config, "placeholder_inc_deg", 60))
  )
  dep <- deproject_coordinates(spaxels$x, spaxels$y, geometry)
  spaxels <- cbind(spaxels, dep)

  class_counts <- as.data.frame(table(spaxels$seg_class[spaxels$valid]), stringsAsFactors = FALSE)
  names(class_counts) <- c("seg_class", "N_valid")
  log_lines <- c(
    log_lines,
    paste("geometry:", paste(names(geometry), unlist(geometry), sep = "=", collapse = "; ")),
    paste("geometry_status:", geometry$geometry_status),
    "valid_spaxels_by_class:",
    paste(class_counts$seg_class, class_counts$N_valid, sep = "=", collapse = "; ")
  )

  n_rings <- as.integer(.capivara_cfg(config, "n_rings", 10))
  axisym <- fit_axisymmetric_piecewise_model(
    spaxels,
    geometry,
    use_errors = isTRUE(.capivara_cfg(config, "use_velocity_errors", TRUE)),
    n_rings = n_rings,
    smooth_lambda = as.numeric(.capivara_cfg(config, "smooth_lambda", 10)),
    fixed_vsys = .capivara_cfg(config, "fixed_vsys", NULL)
  )
  spaxels <- axisym$spaxels

  bar_cfg <- .capivara_cfg(config, "bar", list())
  bar_geometry <- estimate_bar_geometry(
    spaxels,
    phi_b_deg = .capivara_nested_cfg(config, "bar", "phi_b_deg", .capivara_nested_cfg(config, "bar", "bar_pa_deg", NULL))
  )
  fit <- fit_bisymmetric_model(
    spaxels,
    geometry,
    bar_geometry,
    use_errors = isTRUE(.capivara_cfg(config, "use_velocity_errors", TRUE)),
    n_rings = n_rings,
    smooth_lambda = as.numeric(.capivara_cfg(config, "smooth_lambda", 10)),
    second_order_lambda = as.numeric(.capivara_cfg(config, "second_order_lambda", 500)),
    fixed_vsys = .capivara_cfg(config, "fixed_vsys", NULL)
  )
  spaxels <- fit$spaxels

  diagnostics <- compute_residual_diagnostics(
    spaxels,
    plateifu = plateifu,
    fit_parameters = fit$parameters,
    fit_status = fit$fit_status,
    geometry_status = geometry$geometry_status
  )

  summary_file <- file.path(output_dir, "tables", paste0(plateifu, "_kinematic_summary.csv"))
  diag_file <- file.path(output_dir, "tables", paste0(plateifu, "_segment_residual_diagnostics.csv"))
  readr::write_csv(diagnostics$summary, summary_file)
  readr::write_csv(diagnostics$by_class, diag_file)

  result <- list(
    plateifu = plateifu,
    config = config,
    manga = manga,
    segmentation_map = cap$segmentation_map,
    segment_table = cap$segment_table,
    spaxels = spaxels,
    geometry = geometry,
    bar_geometry = bar_geometry,
    axisym = axisym,
    fit = fit,
    diagnostics = diagnostics,
    dims = dim(manga$velocity),
    output_dir = output_dir
  )
  saveRDS(result, file.path(output_dir, "models", paste0(plateifu, "_bisymmetric_model.rds")))

  plot_capivara_kinematics(
    result,
    png_file = file.path(output_dir, "figures", paste0(plateifu, "_bisymmetric_model.png")),
    pdf_file = file.path(output_dir, "figures", paste0(plateifu, "_bisymmetric_model.pdf"))
  )
  plot_capivara_component_decomposition(
    result,
    png_file = file.path(output_dir, "figures", paste0(plateifu, "_capivara_components.png")),
    pdf_file = file.path(output_dir, "figures", paste0(plateifu, "_capivara_components.pdf"))
  )

  report <- .capivara_report_text(diagnostics$summary)
  writeLines(report, file.path(output_dir, "report.txt"))
  log_lines <- c(
    log_lines,
    paste("N_axisym_fit:", axisym$N_fit),
    paste("N_bar_fit:", fit$N_fit),
    paste("bar_geometry:", paste(names(bar_geometry), unlist(bar_geometry), sep = "=", collapse = "; ")),
    paste("fitted_parameters:", paste(names(fit$parameters), unlist(fit$parameters), sep = "=", collapse = "; ")),
    paste("fit_status:", fit$fit_status),
    paste("summary_file:", summary_file),
    paste("diagnostics_file:", diag_file)
  )
  .capivara_write_log(log_file, log_lines)

  invisible(result)
}
