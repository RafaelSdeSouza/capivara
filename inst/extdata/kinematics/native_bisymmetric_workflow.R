env_arg <- function(key, unset) {
  value <- Sys.getenv(key, unset = NA_character_)
  if (!is.na(value) && nzchar(value)) value else unset
}

args <- commandArgs(trailingOnly = TRUE)
use_env_inputs <- tolower(Sys.getenv("CAPIVARA_MODEL_USE_ENV_INPUTS", unset = "false")) %in% c("1", "true", "yes", "y", "on")
native_rds <- if (!use_env_inputs && length(args) >= 1L) args[[1]] else env_arg("CAPIVARA_MODEL_NATIVE_RDS", "")
if (!nzchar(native_rds)) stop("Supply a native result path explicitly.", call. = FALSE)
output_dir <- if (!use_env_inputs && length(args) >= 2L) args[[2]] else env_arg("CAPIVARA_MODEL_OUTPUT_DIR", dirname(native_rds))
prefix <- if (!use_env_inputs && length(args) >= 3L) args[[3]] else env_arg("CAPIVARA_MODEL_OUTPUT_PREFIX", "capivara_model")
plateifu <- Sys.getenv("CAPIVARA_MODEL_PLATEIFU", unset = "native cube")

env_num <- function(key, unset = NA_real_) {
  value <- Sys.getenv(key, unset = NA_character_)
  if (!is.na(value) && nzchar(value)) as.numeric(value) else as.numeric(unset)
}
env_chr <- function(key, unset) {
  value <- Sys.getenv(key, unset = NA_character_)
  if (!is.na(value) && nzchar(value)) value else unset
}
env_bool <- function(key, unset = "false") {
  tolower(env_chr(key, unset)) %in% c("1", "true", "yes", "y", "on")
}
model_kind <- .capivara_match_kinematic_model(
  env_chr("CAPIVARA_MODEL_KIND", "axisymmetric")
)
analysis_mode <- match.arg(env_chr("CAPIVARA_MODEL_ANALYSIS_MODE", "science"), c("science", "preview"))
if (!exists("input_bar_geometry", inherits = FALSE)) input_bar_geometry <- NULL

make_bar_support_mask <- function(spaxels, phi_bar_disc_rad, width_deg = 25, r_min_frac = 0, r_max_frac = 1) {
  out <- rep(FALSE, nrow(spaxels))
  ok <- spaxels$valid & is.finite(spaxels$theta) & is.finite(spaxels$R) & is.finite(phi_bar_disc_rad)
  if (!any(ok)) {
    return(out)
  }

  r_max <- stats::quantile(spaxels$R[ok], 0.98, na.rm = TRUE)
  if (!is.finite(r_max) || r_max <= 0) {
    return(out)
  }

  width_rad <- abs(as.numeric(width_deg)) * pi / 180
  r_min <- max(0, as.numeric(r_min_frac)) * r_max
  r_max_use <- max(r_min, as.numeric(r_max_frac) * r_max)
  dphi <- atan2(sin(spaxels$theta - phi_bar_disc_rad), cos(spaxels$theta - phi_bar_disc_rad))
  dphi <- pmin(abs(dphi), abs(abs(dphi) - pi))
  out[ok] <- dphi[ok] <= width_rad & spaxels$R[ok] >= r_min & spaxels$R[ok] <= r_max_use
  out
}

if (!file.exists(native_rds)) {
  stop("Missing native Capivara kinematic RDS: ", native_rds, call. = FALSE)
}

native <- readRDS(native_rds)
kin <- native$kinematics
dims <- dim(kin$velocity)
support <- native$support
valid <- support & kin$valid & is.finite(kin$velocity)

native_white_light_centre <- function(native, dims) {
  white_light <- native$starlet$collapsed
  if (!is.matrix(white_light) || !identical(dim(white_light), dims) || !any(is.finite(white_light))) {
    return(list(x = NA_real_, y = NA_real_, source = "kinematic_flux_fallback"))
  }
  hit <- which.max(ifelse(is.finite(white_light), white_light, -Inf))
  pixel <- arrayInd(hit, dim(white_light))
  list(x = pixel[2], y = pixel[1], source = "native_white_light_peak")
}

white_centre <- native_white_light_centre(native, dims)
model_x0 <- env_num("CAPIVARA_MODEL_X0", NA_real_)
model_y0 <- env_num("CAPIVARA_MODEL_Y0", NA_real_)
if (!is.finite(model_x0) || !is.finite(model_y0)) {
  model_x0 <- white_centre$x
  model_y0 <- white_centre$y
}

manga <- list(
  velocity = kin$velocity,
  velocity_error = NULL,
  velocity_ivar = NULL,
  mask = ifelse(valid, 0, 1),
  flux = kin$flux,
  snr = NULL,
  line_flux = kin$flux,
  line_sigma = kin$sigma,
  header = NULL,
  maps_file = native_rds,
  velocity_component = paste0("native_", native$line$slug)
)
selected_mode <- env_chr("CAPIVARA_MODEL_SEGMENTATION_MODE", "kinematic")
kinematic_segments <- if (selected_mode == "path_signature") {
  native$path_signature$cluster_map
} else native$kinematic_aware$cluster_map
if (is.null(kinematic_segments) || !identical(dim(kinematic_segments), dims)) {
  kinematic_segments <- matrix(NA_integer_, dims[1], dims[2])
  kinematic_segments[valid] <- 1L
}
segment_labels <- sort(unique(as.integer(stats::na.omit(as.vector(kinematic_segments)))))
segment_labels <- segment_labels[segment_labels > 0L]
cap <- list(
  segmentation_map = kinematic_segments,
  segment_table = data.frame(
    label = segment_labels,
    class = "disc",
    use_for_disc_fit = TRUE,
    use_for_bar_diagnostics = FALSE,
    stringsAsFactors = FALSE
  )
)

spaxels <- .capivara_build_spaxels(manga, cap)
pixel_index <- cbind(spaxels$y, spaxels$x)
spaxels <- .capivara_measurement_flags(spaxels, kin, preview = analysis_mode == "preview")

geometry <- estimate_disc_geometry(
  spaxels,
  geometry = list(
    x0 = model_x0,
    y0 = model_y0,
    vsys = env_num("CAPIVARA_MODEL_VSYS", NA_real_),
    pa_image_deg = env_num("CAPIVARA_MODEL_PA_DEG", NA_real_),
    inc_deg = env_num("CAPIVARA_MODEL_INC_DEG", NA_real_),
    coordinate_convention = "nirvana"
  ),
  allow_placeholder_inclination = analysis_mode == "preview",
  placeholder_inc_deg = 60
)
# Fit in the original IFU array coordinates. Display orientation is handled
# only by the plotting layer and must not rotate the physical model.
spaxels <- cbind(spaxels, deproject_coordinates(spaxels$x, spaxels$y, geometry))

use_bar_mask <- identical(model_kind, "bisymmetric_bar") &&
  env_bool("CAPIVARA_MODEL_USE_BAR_MASK", "false")
bar_mask <- rep(FALSE, nrow(spaxels))
bar_geometry <- NULL
if (identical(model_kind, "bisymmetric_bar")) {
  phi_bar_disc_deg <- env_num("CAPIVARA_MODEL_BAR_PHI_DEG", NA_real_)
  if (!is.null(input_bar_geometry)) {
    bar_geometry <- .capivara_vetted_bar_geometry(input_bar_geometry, geometry, dims,
                                               science = analysis_mode == "science")
    if (is.finite(phi_bar_disc_deg) &&
        abs(.capivara_axial_angle(phi_bar_disc_deg - bar_geometry$phi_bar_disc_deg)) > 1e-6) {
      stop("Manual bar angle conflicts with supplied bar geometry.")
    }
  } else {
    if (analysis_mode == "science" && !is.finite(phi_bar_disc_deg)) {
      stop("Scientific mode requires an explicitly supplied bar angle or vetted geometry.")
    }
    bar_geometry <- estimate_bar_geometry(
    spaxels,
    phi_bar_disc_deg = phi_bar_disc_deg,
    white_light = native$starlet$collapsed,
    geometry = geometry
  )
  }
}
bar_support_source <- "not_supplied"
if (!is.null(bar_geometry$bar_mask)) {
  bar_mask <- bar_geometry$bar_mask[pixel_index]
  bar_support_source <- "supplied_native_grid_mask"
} else if (isTRUE(use_bar_mask)) {
  if (analysis_mode != "preview") {
    stop("An angular wedge is preview-only; supply a vetted bar_geometry$bar_mask.")
  }
  bar_mask <- make_bar_support_mask(
    spaxels,
    phi_bar_disc_rad = bar_geometry$phi_bar_disc_rad,
    width_deg = env_num("CAPIVARA_MODEL_BAR_MASK_WIDTH_DEG", 25),
    r_min_frac = env_num("CAPIVARA_MODEL_BAR_MASK_R_MIN_FRAC", 0),
    r_max_frac = env_num("CAPIVARA_MODEL_BAR_MASK_R_MAX_FRAC", 1)
  )
  bar_support_source <- "preview_angular_wedge"
}
bar_support_available <- bar_support_source != "not_supplied"
if (bar_support_available) {
  spaxels$seg_class[bar_mask] <- "bar"
  spaxels$use_for_bar_diagnostics[bar_mask] <- TRUE
}
n_rings <- as.integer(env_num("CAPIVARA_MODEL_N_RINGS", 22))
robust_fit <- env_bool("CAPIVARA_MODEL_ROBUST", "true")
smooth_lambda <- env_num("CAPIVARA_MODEL_SMOOTH_LAMBDA", 10)
second_order_lambda <- env_num("CAPIVARA_MODEL_SECOND_ORDER_LAMBDA", 25)
axisym <- fit_axisymmetric_piecewise_model(
  spaxels,
  geometry,
  n_rings = n_rings,
  smooth_lambda = smooth_lambda,
  fixed_vsys = geometry$vsys,
  robust = robust_fit
)
spaxels <- axisym$spaxels
if (identical(model_kind, "bisymmetric_bar")) {
  fit <- fit_bisymmetric_model(
    spaxels,
    geometry,
    bar_geometry,
    n_rings = n_rings,
    smooth_lambda = smooth_lambda,
    second_order_lambda = second_order_lambda,
    fixed_vsys = geometry$vsys,
    robust = robust_fit,
    max_v2_fraction = env_num("CAPIVARA_MODEL_MAX_V2_FRACTION", 0.8),
    max_mean_v2 = env_num("CAPIVARA_MODEL_MAX_MEAN_V2", 350)
  )
  spaxels <- fit$spaxels
} else {
  fit <- axisym
  spaxels$v_resid <- spaxels$v_axisym_resid
}

diagnostics <- compute_residual_diagnostics(
  spaxels,
  plateifu = plateifu,
  fit_parameters = fit$parameters,
  fit_status = fit$fit_status,
  geometry_status = geometry$geometry_status,
  bar_support_available = bar_support_available
)

result <- list(
  plateifu = plateifu,
  config = list(
    model = model_kind,
    analysis_mode = analysis_mode,
    display_orientation = env_chr("CAPIVARA_MODEL_DISPLAY_ORIENTATION", "rot90_cw"),
    velocity_source = paste("Capivara-native", native$line$slug),
    support = if (identical(native$support_method, "line_flux")) {
      "connected line-flux kinematic support"
    } else {
      paste0(
        "connected starlet support, scales ",
        Sys.getenv("CAPIVARA_STARLET_SCALES", unset = "2:5")
      )
    },
    segmentation_role = "kinematic regions retained for visualization and diagnostics; velocity model fit is spaxel-wise",
    centre_source = white_centre$source,
    bar_angle_source = if (identical(model_kind, "bisymmetric_bar")) bar_geometry$bar_status else NA_character_,
    imputed_fit_weight = if (analysis_mode == "preview") 0.35 else 0,
    n_rings = n_rings,
    robust_fit = robust_fit,
    smooth_lambda = smooth_lambda,
    second_order_lambda = if (identical(model_kind, "bisymmetric_bar")) second_order_lambda else NA_real_,
    bar_support_mask = bar_support_source,
    bar_support_width_deg = if (use_bar_mask) env_num("CAPIVARA_MODEL_BAR_MASK_WIDTH_DEG", 25) else NA_real_,
    bar_support_n_spaxels = if (bar_support_available) sum(bar_mask & spaxels$valid, na.rm = TRUE) else NA_integer_,
    fixed_vsys = geometry$vsys
  ),
  native = native,
  manga = manga,
  segmentation_map = cap$segmentation_map,
  segment_table = cap$segment_table,
  spaxels = spaxels,
  geometry = geometry,
  bar_geometry = bar_geometry,
  qc = .capivara_kinematic_fit_qc(fit, analysis_mode, bar_support_available),
  axisym = axisym,
  fit = fit,
  diagnostics = diagnostics,
  dims = dims,
  output_dir = output_dir
)

.capivara_dir_create(output_dir)
saveRDS(result, file.path(output_dir, paste0(prefix, ".rds")))
readr::write_csv(diagnostics$summary, file.path(output_dir, paste0(prefix, "_summary.csv")))
readr::write_csv(diagnostics$by_class, file.path(output_dir, paste0(prefix, "_residuals_by_class.csv")))

model_plot <- plot_capivara_kinematics(
  result,
  png_file = file.path(output_dir, paste0(prefix, "_model.png")),
  pdf_file = file.path(output_dir, paste0(prefix, "_model.pdf"))
)
if (identical(model_kind, "bisymmetric_bar")) {
  component_plot <- plot_capivara_component_decomposition(
    result,
    png_file = file.path(output_dir, paste0(prefix, "_components.png")),
    pdf_file = file.path(output_dir, paste0(prefix, "_components.pdf"))
  )
}
