# CAPIVARA scientific plotting defaults. These functions affect display only.

.capivara_ink <- "#322825"
.capivara_paper <- "#FCFCFA"
.capivara_invalid <- "#DEDCD5"

.capivara_sequential_palette <- function(n = 256L, direction = 1L) {
  stops <- c("#1D285A", "#2B445D", "#315E90", "#3E7DA8",
             "#709B96", "#C7AD33", "#D7D14D", "#ECE8A2")
  values <- grDevices::colorRampPalette(stops, space = "Lab")(n)
  if (direction < 0) rev(values) else values
}

.capivara_identity_palette <- function(n) {
  stops <- c("#1E295B", "#E3D75A", "#4E74A6", "#B58B2B",
             "#68A69C", "#C96442", "#80538A", "#44312A",
             "#82B7D3", "#D6A65A", "#2D6B68", "#A74E62")
  if (n <= length(stops)) return(stops[seq_len(n)])
  grDevices::colorRampPalette(stops, space = "Lab")(n)
}

.capivara_diverging_palette <- function(n = 257L, direction = 1L) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 2L) {
    stop("`n` must be an integer of at least two.", call. = FALSE)
  }
  neutral <- "#F4F0DC"
  n_left <- ceiling((n + 1L) / 2L)
  n_right <- n - n_left + 1L
  left <- grDevices::colorRampPalette(
    c("#1E295B", "#4E74A6", "#79A49E", neutral), space = "Lab"
  )(n_left)
  right <- grDevices::colorRampPalette(
    c(neutral, "#DCD655", "#BFA524", "#7A4932", "#322825"),
    space = "Lab"
  )(n_right)
  left[n_left] <- neutral
  right[1L] <- neutral
  values <- c(left, right[-1L])
  if (direction < 0) rev(values) else values
}

.capivara_continuous_scale <- function(aesthetics, colours, name, limits,
                                        oob, na.value, guide, trans, ...) {
  ggplot2::continuous_scale(
    aesthetics = aesthetics,
    palette = scales::gradient_n_pal(colours),
    name = name, limits = limits, oob = oob, na.value = na.value,
    guide = guide, transform = trans, ...
  )
}

#' CAPIVARA Van Gogh 2.0 continuous scale
#'
#' Use for continuous unsigned quantities such as continuum, flux, age, and
#' metallicity. Styling changes display only.
#'
#' @param aesthetics Aesthetic to scale, normally `"fill"` or `"colour"`.
#' @param name Legend title.
#' @param limits Optional scale limits.
#' @param direction Palette direction (`1` or `-1`).
#' @param na.value Colour for invalid or unmeasured values.
#' @param oob Out-of-bounds handler.
#' @param guide Guide specification.
#' @param trans Scale transformation.
#' @param ... Additional arguments passed to [ggplot2::continuous_scale()].
#' @return A ggplot2 continuous scale.
#' @export
scale_capivara_continuous <- function(
    aesthetics = "fill", name = ggplot2::waiver(), limits = NULL,
    direction = 1L, na.value = .capivara_invalid, oob = scales::squish,
    guide = "colourbar", trans = "identity", ...) {
  .capivara_continuous_scale(
    aesthetics, .capivara_sequential_palette(256L, direction), name,
    limits, oob, na.value, guide, trans, ...
  )
}

#' CAPIVARA Van Gogh 2.0 ordered-region scale
#'
#' @param n Number of ranks.
#' @inheritParams scale_capivara_continuous
#' @return A ggplot2 ordered continuous scale.
#' @export
scale_capivara_rank <- function(
    n, aesthetics = "fill", name = ggplot2::waiver(), direction = 1L,
    na.value = .capivara_invalid, oob = scales::squish,
    guide = "coloursteps", trans = "identity", ...) {
  n <- as.integer(n)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be one positive integer.", call. = FALSE)
  }
  .capivara_continuous_scale(
    aesthetics, .capivara_sequential_palette(n, direction), name,
    c(1, n), oob, na.value, guide, trans, ...
  )
}

#' CAPIVARA Van Gogh 2.0 region-identity scale
#'
#' @param aesthetics Aesthetic to scale.
#' @param name Legend title.
#' @param na.value Colour for invalid or unmeasured values.
#' @param guide Guide specification.
#' @param ... Additional arguments passed to [ggplot2::discrete_scale()].
#' @return A ggplot2 discrete scale.
#' @export
scale_capivara_identity <- function(
    aesthetics = "fill", name = ggplot2::waiver(),
    na.value = .capivara_invalid, guide = "legend", ...) {
  ggplot2::discrete_scale(
    aesthetics = aesthetics, palette = .capivara_identity_palette, name = name,
    na.value = na.value, guide = guide, ...
  )
}

#' CAPIVARA Van Gogh 2.0 zero-centred diverging scale
#'
#' @param midpoint Value mapped to the neutral centre.
#' @inheritParams scale_capivara_continuous
#' @return A ggplot2 diverging continuous scale.
#' @export
scale_capivara_diverging <- function(
    aesthetics = "fill", name = ggplot2::waiver(), midpoint = 0,
    limits = NULL, direction = 1L, na.value = .capivara_invalid,
    oob = scales::squish, guide = "colourbar", trans = "identity", ...) {
  scale <- .capivara_continuous_scale(
    aesthetics, .capivara_diverging_palette(257L, direction), name,
    limits, oob, na.value, guide, trans, ...
  )
  scale$rescaler <- function(x, from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to = c(0, 1), from = from, mid = midpoint)
  }
  scale
}

#' CAPIVARA publication theme
#'
#' @param style One of `"paper"`, `"talk"`, `"dark"`, or `"spatial"`.
#' @param base_size Base text size in points. The style default is used when
#'   `NULL`.
#' @param base_family Font family.
#' @return A ggplot2 theme.
#' @export
theme_capivara <- function(style = c("paper", "talk", "dark", "spatial"),
                           base_size = NULL, base_family = "sans") {
  style <- match.arg(style)
  if (is.null(base_size)) {
    base_size <- switch(style, paper = 8.5, spatial = 8.5,
                        talk = 16, dark = 13)
  }
  dark <- identical(style, "dark")
  paper <- if (dark) "#11161C" else .capivara_paper
  ink <- if (dark) "#E9EDF1" else "#20252B"
  muted <- if (dark) "#A9B1BA" else "#66717D"
  talk <- identical(style, "talk")
  linewidth <- if (talk) 0.52 else if (dark) 0.42 else 0.32
  tick_length <- if (talk) 4.2 else if (dark) 3.2 else 2.3
  ggplot2::theme(
    line = ggplot2::element_line(colour = ink, linewidth = 0.32,
                                 lineend = "square"),
    rect = ggplot2::element_rect(fill = paper, colour = ink,
                                 linewidth = linewidth),
    text = ggplot2::element_text(family = base_family, colour = ink,
                                 size = base_size, lineheight = 0.96),
    axis.line = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = ggplot2::rel(0.96),
                                       margin = ggplot2::margin(2.5, 2.5, 2.5, 2.5)),
    axis.text = ggplot2::element_text(size = ggplot2::rel(0.86), colour = ink),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 2)),
    axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
    axis.ticks = ggplot2::element_line(colour = ink, linewidth = linewidth),
    axis.ticks.length = grid::unit(tick_length, "pt"),
    panel.background = ggplot2::element_rect(fill = paper, colour = NA),
    panel.border = ggplot2::element_rect(fill = NA, colour = muted,
                                         linewidth = linewidth),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.background = ggplot2::element_rect(fill = paper, colour = NA),
    plot.margin = grid::unit(if (talk) c(9, 10, 8, 10) else c(5.5, 5.5, 4.5, 5.5), "pt"),
    strip.background = ggplot2::element_rect(fill = paper, colour = muted,
                                              linewidth = linewidth),
    strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.88),
                                       margin = ggplot2::margin(3, 3, 3, 3)),
    legend.background = ggplot2::element_rect(fill = paper, colour = NA),
    legend.key = ggplot2::element_rect(fill = paper, colour = NA),
    legend.key.height = grid::unit(if (talk) 18 else 10, "pt"),
    legend.key.width = grid::unit(if (talk) 12 else 7.5, "pt"),
    legend.title = ggplot2::element_text(size = ggplot2::rel(0.82), face = "bold"),
    legend.text = ggplot2::element_text(size = ggplot2::rel(0.76))
  )
}
