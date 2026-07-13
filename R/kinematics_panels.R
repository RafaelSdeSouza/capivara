#' Extract individual kinematic plot panels
#'
#' The standard `plot()` method creates a compact multi-panel figure. This
#' helper returns the individual ggplot objects so they can be saved, restyled,
#' or assembled into a paper-specific layout without rerunning the analysis.
#'
#' For an axisymmetric model, `view = "model"` returns `footprint`, `velocity`,
#' `disc_model`, `disc_residual`, `circular_speed`, and `summary`. For a
#' bisymmetric bar model it returns the corresponding full-model, axisymmetric,
#' residual, and velocity-profile panels. `view = "components"` returns the
#' individual circular and bar-flow component maps, along with the supporting
#' line-property maps and residual comparisons.
#'
#' @param result A result from [run_kinematic_analysis()] or
#'   [run_manga_bar_model()].
#' @param view `"model"` extracts the disc or bisymmetric-model panels.
#'   `"components"` extracts the bar-component panels and is available only for
#'   `model = "bisymmetric_bar"`.
#' @return A named list of ggplot objects.
#' @export
kinematic_panels <- function(result, view = c("model", "components")) {
  view <- match.arg(view)
  if (is.null(result$model_plot) || is.null(result$model)) {
    stop("`result` must be returned by `run_kinematic_analysis()`.", call. = FALSE)
  }
  model <- .capivara_match_kinematic_model(result$model)
  if (identical(view, "components") && !identical(model, "bisymmetric_bar")) {
    stop(
      "Component panels are available only for `model = \"bisymmetric_bar\"`.",
      call. = FALSE
    )
  }

  plot_object <- if (identical(view, "model")) result$model_plot else result$component_plot
  if (is.null(plot_object)) {
    stop("The requested panels were not created for this result.", call. = FALSE)
  }
  panels <- attr(plot_object, "capivara_panels", exact = TRUE)
  if (is.null(panels)) {
    stop(
      "This result was created by an older Capivara version. Rerun the analysis to extract individual panels.",
      call. = FALSE
    )
  }
  panels
}
