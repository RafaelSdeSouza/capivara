#' Compute residual diagnostics by Capivara segment class
#'
#' @param spaxels Spaxel data frame with `v_resid`, `seg_class`, and flags.
#' @param plateifu Optional plate-IFU identifier.
#' @param fit_parameters Optional one-row parameter data frame.
#' @param fit_status Fit status string.
#' @param geometry_status Geometry status string.
#' @return A list with `by_class` and one-row `summary`.
#' @noRd
compute_residual_diagnostics <- function(spaxels,
                                         plateifu = NA_character_,
                                         fit_parameters = NULL,
                                         fit_status = NA_character_,
                                         geometry_status = NA_character_,
                                         bar_support_available = NULL) {
  if (is.null(bar_support_available)) {
    bar_support_available <- any(tolower(spaxels$seg_class) == "bar", na.rm = TRUE)
  }
  valid <- spaxels$valid & is.finite(spaxels$v_resid)
  df <- spaxels[valid, , drop = FALSE]
  classes <- sort(unique(df$seg_class))
  by_class <- do.call(rbind, lapply(classes, function(cls) {
    z <- df$v_resid[df$seg_class == cls]
    data.frame(
      seg_class = cls,
      N = length(z),
      median_resid = stats::median(z, na.rm = TRUE),
      mad_resid = stats::mad(z, constant = 1.4826, na.rm = TRUE),
      rms_resid = sqrt(mean(z^2, na.rm = TRUE)),
      A_resid = stats::median(abs(z), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  metric <- function(cls, col) {
    hit <- match(cls, by_class$seg_class)
    if (is.na(hit)) NA_real_ else by_class[[col]][hit]
  }
  N_class <- function(cls) sum(df$seg_class == cls, na.rm = TRUE)
  bar <- df[tolower(df$seg_class) == "bar", , drop = FALSE]
  total_power <- sum(df$v_resid^2, na.rm = TRUE)
  bar_power <- sum(bar$v_resid^2, na.rm = TRUE)

  par <- if (is.null(fit_parameters)) {
    data.frame(vsys = NA_real_, vmax = NA_real_, Rt = NA_real_)
  } else {
    fit_parameters
  }
  summary <- data.frame(
    plateifu = plateifu,
    N_valid = nrow(df),
    N_disc = N_class("disc"),
    N_bar = if (bar_support_available) N_class("bar") else NA_integer_,
    N_ring = N_class("ring"),
    N_nucleus = N_class("nucleus"),
    vsys = par$vsys[1],
    vmax = par$vmax[1],
    Rt = par$Rt[1],
    A_bar = if (bar_support_available) metric("bar", "A_resid") else NA_real_,
    RMS_disc = metric("disc", "rms_resid"),
    RMS_bar = if (bar_support_available) metric("bar", "rms_resid") else NA_real_,
    Q_kin = if (bar_support_available && is.finite(metric("disc", "rms_resid")) &&
      metric("disc", "rms_resid") > 0) metric("bar", "rms_resid") / metric("disc", "rms_resid") else NA_real_,
    f_bar = if (bar_support_available && nrow(bar) > 0L && is.finite(total_power) && total_power > 0) bar_power / total_power else NA_real_,
    bar_diagnostic_qc = if (!bar_support_available) "bar_support_not_supplied" else
      if (!nrow(bar)) "no_valid_spaxels_in_supplied_bar_support" else "defined_on_supplied_support",
    fit_status = fit_status,
    geometry_status = geometry_status,
    stringsAsFactors = FALSE
  )

  list(by_class = by_class, summary = summary)
}
