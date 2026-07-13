.capivara_display_matrix <- function(mat, orientation = "transpose") {
  if (is.null(orientation) || !length(orientation)) {
    orientation <- "transpose"
  }
  orientation <- tolower(as.character(orientation[[1]]))
  switch(
    orientation,
    identity = mat,
    transpose = t(mat),
    flip_x = mat[, ncol(mat):1, drop = FALSE],
    flip_y = mat[nrow(mat):1, , drop = FALSE],
    rot90_cw = t(mat[nrow(mat):1, , drop = FALSE]),
    rot90_ccw = t(mat)[ncol(mat):1, , drop = FALSE],
    rot180 = mat[nrow(mat):1, ncol(mat):1, drop = FALSE],
    stop("Unsupported display orientation: ", orientation, call. = FALSE)
  )
}

.capivara_squish_oob <- function(x, range) {
  pmin(pmax(x, range[1]), range[2])
}

.capivara_map_plot <- function(mat, title, fill_label = NULL, diverging = FALSE,
                               discrete = FALSE, limits = NULL, legend_position = "bottom",
                               display_orientation = "transpose") {
  path_div_palette <- c("#082A55", "#145A96", "#5EA4C8", "#F7F3E8", "#F2A15F", "#C43D2E", "#6E1419")
  path_flux_palette <- c("#101D3A", "#174A7C", "#1E7A9C", "#36A793", "#B8D66A", "#F0D343", "#F59A2F", "#C94E27", "#7F1D1D")
  path_segment_palette <- c("#082A55", "#145A96", "#1E7A9C", "#36A793", "#B8D66A", "#F0D343", "#F59A2F", "#C94E27", "#7F1D1D", "#6B3FA0", "#157A62", "#C9A227")
  mat <- .capivara_display_matrix(mat, display_orientation)
  df <- matrix_to_long(mat)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.height = ggplot2::unit(0.11, "in"),
      legend.key.width = ggplot2::unit(0.35, "in"),
      plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
      panel.border = ggplot2::element_rect(colour = "grey25", fill = NA, linewidth = 0.35),
      plot.margin = ggplot2::margin(4, 4, 10, 4)
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      title.position = "left",
      title.vjust = 0.5,
      barwidth = ggplot2::unit(1.5, "in"),
      barheight = ggplot2::unit(0.09, "in")
    )) +
    ggplot2::labs(title = title, fill = fill_label)

  if (isTRUE(discrete)) {
    df$value <- factor(df$value)
    n <- max(2L, length(levels(stats::na.omit(df$value))))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::geom_raster() +
      ggplot2::coord_fixed(expand = FALSE) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_fill_manual(
        values = grDevices::colorRampPalette(path_segment_palette, space = "Lab")(n),
        na.value = "white"
      ) +
      ggplot2::theme_void(base_size = 10) +
      ggplot2::theme(
        legend.position = legend_position,
        legend.direction = "horizontal",
        legend.title = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 8),
        legend.key.height = ggplot2::unit(0.11, "in"),
        legend.key.width = ggplot2::unit(0.35, "in"),
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
        panel.border = ggplot2::element_rect(colour = "grey25", fill = NA, linewidth = 0.35),
        plot.margin = ggplot2::margin(4, 4, 10, 4)
      ) +
      ggplot2::labs(title = title, fill = fill_label)
  } else if (isTRUE(diverging)) {
    if (is.null(limits) || !all(is.finite(limits)) || diff(limits) <= 0) {
      limits <- range(df$value, finite = TRUE)
    }
    zero_value <- (0 - limits[1]) / diff(limits)
    zero_value <- min(1, max(0, zero_value))
    neutral <- ceiling(length(path_div_palette) / 2)
    div_values <- c(
      seq(0, zero_value, length.out = neutral),
      seq(zero_value, 1, length.out = length(path_div_palette) - neutral + 1)[-1]
    )
    p <- p + ggplot2::scale_fill_gradientn(
      colours = path_div_palette,
      values = div_values,
      limits = limits,
      oob = .capivara_squish_oob,
      na.value = "grey92"
    )
  } else {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = path_flux_palette,
      limits = limits,
      oob = .capivara_squish_oob,
      na.value = "grey92"
    )
  }
  p
}

#' Plot Capivara segmentation-aware kinematic modelling products
#'
#' @param result Result from `run_one_galaxy()`.
#' @param png_file,pdf_file Output paths.
#' @return The patchwork plot object.
#' @export
plot_capivara_kinematics <- function(result, png_file = NULL, pdf_file = NULL) {
  .capivara_require(c("ggplot2", "patchwork"))

  display_orientation <- result$config$display_orientation
  if (is.null(display_orientation) || !length(display_orientation)) {
    display_orientation <- "transpose"
  }
  dims <- result$dims
  to_mat <- function(value) {
    mat <- matrix(NA_real_, dims[1], dims[2])
    idx <- cbind(result$spaxels$y, result$spaxels$x)
    mat[idx] <- value
    mat
  }

  flux <- if (!is.null(result$manga$flux)) result$manga$flux else result$manga$snr
  if (is.null(flux)) {
    flux <- matrix(as.numeric(result$spaxels$valid), dims[1], dims[2])
  }
  valid_footprint <- to_mat(as.numeric(result$spaxels$valid))
  flux[!is.finite(valid_footprint) | valid_footprint <= 0] <- NA_real_
  observed <- to_mat(result$spaxels$velocity)
  bar_model <- to_mat(result$spaxels$v_bar_model)
  axisym_model <- to_mat(result$spaxels$v_axisym_model)
  bar_resid <- to_mat(result$spaxels$v_bar_resid)
  axisym_resid <- to_mat(result$spaxels$v_axisym_resid)

  vabs <- stats::quantile(abs(c(observed, bar_model, axisym_model)), 0.995, na.rm = TRUE)
  vel_lim <- c(-vabs, vabs)
  res_lim <- stats::quantile(abs(c(bar_resid, axisym_resid)), 0.98, na.rm = TRUE)
  res_lim <- if (is.finite(res_lim)) c(-res_lim, res_lim) else NULL

  prof <- result$fit$profile
  prof$Axy <- approx(result$axisym$profile$R, result$axisym$profile$Vt, xout = prof$R, rule = 2)$y
  profile_sign <- sign(stats::median(prof$Vt, na.rm = TRUE))
  if (!is.finite(profile_sign) || profile_sign == 0) {
    profile_sign <- 1
  }
  prof$Vt <- prof$Vt * profile_sign
  prof$V2t <- prof$V2t * profile_sign
  prof$V2r <- prof$V2r * profile_sign
  prof$Axy <- prof$Axy * profile_sign
  curve <- ggplot2::ggplot(prof, ggplot2::aes(R)) +
    ggplot2::geom_line(ggplot2::aes(y = Vt, colour = "Vt"), linewidth = 0.75) +
    ggplot2::geom_line(ggplot2::aes(y = V2t, colour = "V2t"), linewidth = 0.75, linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = V2r, colour = "V2r"), linewidth = 0.75, linetype = "dashed") +
    ggplot2::geom_line(ggplot2::aes(y = Axy, colour = "Axy."), linewidth = 0.75, linetype = "dotdash") +
    ggplot2::scale_colour_manual(
      values = c(Vt = "#17223B", V2t = "#C85B3C", V2r = "#6E8F28", "Axy." = "#3158A4"),
      breaks = c("Vt", "V2t", "V2r", "Axy.")
    ) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::theme(legend.position = c(0.78, 0.32), legend.title = ggplot2::element_blank()) +
    ggplot2::labs(title = "Velocity Profiles", x = "Radius (spaxels)", y = "amplitude (km/s)")

  summary <- result$diagnostics$summary[1, ]
  param_text <- paste0(
    "i: ", round(result$geometry$inc_deg, 1), " deg.\n",
    "phi: ", round(result$geometry$pa_deg, 1), " deg.\n",
    "phi_b: ", round(result$bar_geometry$phi_b_deg, 1), " deg.\n",
    "v_sys: ", round(result$fit$parameters$vsys[1], 1), " km/s\n",
    "mean |V2|: ", round(result$fit$parameters$mean_V2[1], 1), " km/s\n",
    "Q_kin: ", round(summary$Q_kin, 2), "\n",
    "fit: ", result$fit$fit_status
  )
  text_panel <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 1, label = param_text, hjust = 0, vjust = 1, size = 3.5, lineheight = 1.15) +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_void()

  panel <- (
    .capivara_map_plot(log10(pmax(flux, 0) + 1e-3), paste(result$plateifu), NULL, legend_position = "none", display_orientation = display_orientation) |
      .capivara_map_plot(observed, "Halpha Velocity", "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation)
  ) / (
    .capivara_map_plot(bar_model, "Bisymmetric model", "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation) |
      .capivara_map_plot(axisym_model, "Axisym. Model", "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation)
  ) / (
    .capivara_map_plot(bar_resid, "Bisymmetric residual", "km/s", diverging = TRUE, limits = res_lim, display_orientation = display_orientation) |
      .capivara_map_plot(axisym_resid, "Axisym. Residual", "km/s", diverging = TRUE, limits = res_lim, display_orientation = display_orientation)
  ) / (curve | text_panel) +
    patchwork::plot_layout(heights = c(1, 1, 1, 0.78)) +
    patchwork::plot_annotation(title = paste("Bisymmetric kinematic model:", result$plateifu)) &
    ggplot2::theme(plot.margin = ggplot2::margin(6, 6, 6, 6))

  if (!is.null(png_file)) {
    .capivara_dir_create(dirname(png_file))
    ggplot2::ggsave(png_file, panel, width = 7.4, height = 12.0, dpi = 320, bg = "white")
  }
  if (!is.null(pdf_file)) {
    .capivara_dir_create(dirname(pdf_file))
    ggplot2::ggsave(pdf_file, panel, width = 7.4, height = 12.0, bg = "white")
  }
  panel
}

#' Plot a paper-style bisymmetric component decomposition panel
#'
#' @param result Result from `run_one_galaxy()`.
#' @param png_file,pdf_file Output paths.
#' @return The patchwork plot object.
#' @export
plot_capivara_component_decomposition <- function(result, png_file = NULL, pdf_file = NULL) {
  .capivara_require(c("ggplot2", "patchwork"))

  display_orientation <- result$config$display_orientation
  if (is.null(display_orientation) || !length(display_orientation)) {
    display_orientation <- "transpose"
  }
  dims <- result$dims
  to_mat <- function(value) {
    mat <- matrix(NA_real_, dims[1], dims[2])
    idx <- cbind(result$spaxels$y, result$spaxels$x)
    mat[idx] <- value
    mat
  }

  valid_footprint <- to_mat(as.numeric(result$spaxels$valid))
  mask_to_valid <- function(mat) {
    mat[!is.finite(valid_footprint) | valid_footprint <= 0] <- NA_real_
    mat
  }

  observed <- to_mat(result$spaxels$velocity)
  model <- to_mat(result$spaxels$v_bar_model)
  vt <- to_mat(result$spaxels$v_t_component)
  v2t <- to_mat(result$spaxels$v_2t_component)
  v2r <- to_mat(result$spaxels$v_2r_component)
  resid <- to_mat(result$spaxels$v_bar_resid)
  no_v2 <- observed - v2t - v2r
  no_vt_v2r <- observed - vt - v2r
  no_vt_v2t <- observed - vt - v2t

  sigma <- result$manga$line_sigma
  if (is.null(sigma)) {
    sigma <- matrix(NA_real_, dims[1], dims[2])
  }
  sigma <- mask_to_valid(sigma)
  surface <- result$manga$line_flux
  if (is.null(surface)) {
    surface <- result$manga$flux
  }
  if (is.null(surface)) {
    surface <- matrix(NA_real_, dims[1], dims[2])
  }
  surface <- mask_to_valid(surface)
  bar_support <- to_mat(as.numeric(tolower(result$spaxels$seg_class) == "bar" & result$spaxels$valid))
  has_bar_support <- sum(bar_support > 0, na.rm = TRUE) > 0
  bar_panel <- if (has_bar_support) {
    .capivara_map_plot(bar_support, "Bar support mask", "bar", discrete = TRUE, display_orientation = display_orientation)
  } else {
    bar_text <- paste0(
      "Bar angle prior\n",
      "phi_b = ", round(result$bar_geometry$phi_b_deg, 1), " deg\n",
      "no support mask"
    )
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.55, label = bar_text, hjust = 0.5, vjust = 0.5, size = 4, lineheight = 1.2) +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, 1) +
      ggplot2::theme_void() +
      ggplot2::theme(
        panel.border = ggplot2::element_rect(colour = "grey25", fill = NA, linewidth = 0.35),
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5)
      ) +
      ggplot2::labs(title = "Bar prior")
  }

  vabs <- stats::quantile(abs(c(observed, model, vt)), 0.995, na.rm = TRUE)
  vel_lim <- c(-vabs, vabs)
  comp_abs <- stats::quantile(abs(c(v2t, v2r)), 0.995, na.rm = TRUE)
  comp_lim <- c(-comp_abs, comp_abs)
  res_abs <- stats::quantile(abs(c(resid, no_vt_v2r, no_vt_v2t)), 0.98, na.rm = TRUE)
  res_lim <- c(-res_abs, res_abs)
  sigma_lim <- c(0, stats::quantile(sigma, 0.995, na.rm = TRUE))
  surface_lim <- stats::quantile(surface, c(0.02, 0.995), na.rm = TRUE)

  panel <- (
    .capivara_map_plot(observed, "Halpha Velocity", "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation) |
      .capivara_map_plot(sigma, "Halpha Dispersion", "km/s", limits = sigma_lim, display_orientation = display_orientation) |
      .capivara_map_plot(log10(pmax(surface, 0) + 1e-4), "Halpha Surface Brightness", "log flux", limits = log10(pmax(surface_lim, 0) + 1e-4), display_orientation = display_orientation) |
      bar_panel
  ) / (
    .capivara_map_plot(model, "Model", "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation) |
      .capivara_map_plot(vt, expression(V[t]), "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation) |
      .capivara_map_plot(v2t, expression(V[2 * t]), "km/s", diverging = TRUE, limits = comp_lim, display_orientation = display_orientation) |
      .capivara_map_plot(v2r, expression(V[2 * r]), "km/s", diverging = TRUE, limits = comp_lim, display_orientation = display_orientation)
  ) / (
    .capivara_map_plot(resid, "Residual", "km/s", diverging = TRUE, limits = res_lim, display_orientation = display_orientation) |
      .capivara_map_plot(no_v2, expression(Data - (V[2 * t] + V[2 * r])), "km/s", diverging = TRUE, limits = vel_lim, display_orientation = display_orientation) |
      .capivara_map_plot(no_vt_v2r, expression(Data - (V[t] + V[2 * r])), "km/s", diverging = TRUE, limits = res_lim, display_orientation = display_orientation) |
      .capivara_map_plot(no_vt_v2t, expression(Data - (V[t] + V[2 * t])), "km/s", diverging = TRUE, limits = res_lim, display_orientation = display_orientation)
  ) +
    patchwork::plot_annotation(title = paste("Bisymmetric component decomposition:", result$plateifu)) &
    ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))

  if (!is.null(png_file)) {
    .capivara_dir_create(dirname(png_file))
    ggplot2::ggsave(png_file, panel, width = 11.6, height = 9.2, dpi = 320, bg = "white")
  }
  if (!is.null(pdf_file)) {
    .capivara_dir_create(dirname(pdf_file))
    ggplot2::ggsave(pdf_file, panel, width = 11.6, height = 9.2, bg = "white")
  }
  panel
}
