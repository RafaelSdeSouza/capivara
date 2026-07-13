#' List the installed kinematic model modules
#'
#' Capivara keeps emission-line measurements and kinematic segmentation
#' independent from dynamical interpretation. A model module is selected only
#' after the velocity maps have been made. The axisymmetric disc is the neutral
#' comparison model; a bisymmetric bar is an explicit additional hypothesis.
#'
#' @return A data frame describing the available model modules.
#' @export
kinematic_models <- function() {
  data.frame(
    model = c("axisymmetric", "bisymmetric_bar"),
    label = c("Axisymmetric disc", "Bisymmetric bar"),
    requires_bar_angle = c(FALSE, TRUE),
    components_plot = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
}

.capivara_match_kinematic_model <- function(model) {
  if (is.null(model) || !length(model) || is.na(model[[1]])) {
    model <- "axisymmetric"
  }
  model <- tolower(trimws(as.character(model[[1]])))
  aliases <- c(
    disc = "axisymmetric",
    axisym = "axisymmetric",
    bisymmetric = "bisymmetric_bar",
    bar = "bisymmetric_bar"
  )
  if (model %in% names(aliases)) {
    model <- unname(aliases[[model]])
  }
  available <- kinematic_models()$model
  if (!model %in% available) {
    stop(
      "Unknown kinematic model `", model, "`. Choose one of: ",
      paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }
  model
}

.capivara_is_bar_model <- function(model) {
  identical(.capivara_match_kinematic_model(model), "bisymmetric_bar")
}
