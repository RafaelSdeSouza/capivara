.starry_night_palette <- function(n) {
  grDevices::colorRampPalette(
    c("#80B7FF", "#547FFF", "#405CFF", "#263C8B",
      "#FFFAA3", "#FFDE38", "#BFA524"),
    space = "Lab"
  )(n)
}

.capivara_east_left_display_sfg <- function(geometry) {
  reflect <- function(coordinates) {
    coordinates <- as.matrix(coordinates)
    coordinates[, 1L] <- -coordinates[, 1L]
    coordinates
  }
  type <- as.character(sf::st_geometry_type(sf::st_sfc(geometry)))[1L]
  switch(type,
    POINT = sf::st_point(as.numeric(reflect(matrix(geometry, nrow = 1L)))),
    MULTIPOINT = sf::st_multipoint(reflect(geometry)),
    LINESTRING = sf::st_linestring(reflect(geometry)),
    MULTILINESTRING = sf::st_multilinestring(lapply(geometry, reflect)),
    POLYGON = sf::st_polygon(lapply(geometry, reflect)),
    MULTIPOLYGON = sf::st_multipolygon(lapply(geometry, function(polygon) {
      lapply(polygon, reflect)
    })),
    GEOMETRYCOLLECTION = sf::st_geometrycollection(
      lapply(geometry, .capivara_east_left_display_sfg)
    ),
    stop("Unsupported geometry type for East-left display: ", type,
         call. = FALSE)
  )
}

.capivara_east_left_display_sf <- function(x) {
  canonical <- sf::st_geometry(x)
  out <- x
  sf::st_geometry(out) <- sf::st_sfc(
    lapply(canonical, .capivara_east_left_display_sfg),
    crs = sf::st_crs(canonical), precision = sf::st_precision(canonical)
  )
  attr(out, "display_transform") <- list(
    name = "astronomical_east_left",
    canonical = c(x = "+East", y = "+North"),
    display = c(x = "-East", y = "+North"),
    east_left_reflections = 1L
  )
  out
}

#' Plot exact CAPIVARA region geometry
#'
#' The standard renderer draws exact cell-union `sf` regions. Scientific
#' geometry always stores positive East and positive North. For sky or physical
#' coordinates, the renderer reflects East exactly once in a temporary display
#' copy so that North is up and East is left.
#'
#' @param cluster_data A CAPIVARA segmentation result.
#' @param palette `"capivara"`, legacy `"starry_night"`, a viridis option,
#'   or a character vector. When omitted, explicit semantic results use the
#'   CAPIVARA palette and historical `representation = NULL` results retain
#'   the legacy palette.
#' @param coords One of `"auto"`, `"native"`, `"sky"`, or `"physical"`.
#'   `"auto"` selects physical coordinates when WCS and redshift are present,
#'   sky coordinates when only WCS is present, and native coordinates otherwise.
#' @param mode `"identity"` for nominal regions or `"rank"` for a declared
#'   ordered scalar.
#' @param rank_by For rank mode, a numeric vector with one value per region or
#'   the name of a numeric column in [regions_sf()].
#' @param east_left Apply the astronomical East-left display convention for sky
#'   and physical coordinates.
#' @param show_legend Display the fill legend.
#' @param boundary_colour Exact region-boundary colour.
#' @param boundary_linewidth Boundary width in millimetres.
#' @return A ggplot2 object whose first layer is [ggplot2::geom_sf()].
#' @export
plot_cluster <- function(
    cluster_data, palette = "starry_night",
    coords = c("auto", "native", "sky", "physical"),
    mode = c("identity", "rank"), rank_by = NULL,
    east_left = TRUE, show_legend = FALSE,
    boundary_colour = .capivara_ink, boundary_linewidth = 0.25) {
  if (!is.list(cluster_data) || !is.matrix(cluster_data$cluster_map)) {
    stop("`cluster_data` must contain a label matrix named `cluster_map`.",
         call. = FALSE)
  }
  coords <- match.arg(coords)
  mode <- match.arg(mode)
  if (missing(palette) && !is.null(cluster_data$representation)) {
    palette <- "capivara"
  }
  contract <- .capivara_object_contract(cluster_data)
  coordinate_mode <- .capivara_spatial_mode(contract, coords)
  regions <- regions_sf(cluster_data, coordinate_mode)
  canonical_geometry <- serialize(sf::st_geometry(regions), NULL, version = 3)

  display_reflections <- 0L
  if (isTRUE(east_left) && coordinate_mode %in% c("sky", "physical")) {
    regions <- .capivara_east_left_display_sf(regions)
    display_reflections <- 1L
  }

  if (mode == "identity") {
    regions$capivara_fill <- factor(regions$region_id)
  } else {
    if (is.character(rank_by) && length(rank_by) == 1L) {
      if (!rank_by %in% names(regions)) {
        stop("`rank_by` does not name a region geometry column.", call. = FALSE)
      }
      rank_by <- regions[[rank_by]]
    }
    if (!is.numeric(rank_by) || length(rank_by) != nrow(regions) ||
        any(!is.finite(rank_by))) {
      stop("Rank mode requires one finite declared scalar per region.",
           call. = FALSE)
    }
    regions$capivara_fill <- rank(rank_by, ties.method = "first")
  }

  plot <- ggplot2::ggplot(regions) +
    ggplot2::geom_sf(
      ggplot2::aes(fill = capivara_fill), colour = boundary_colour,
      linewidth = boundary_linewidth, linejoin = "mitre"
    )

  if (mode == "rank") {
    plot <- plot + scale_capivara_rank(nrow(regions))
  } else if (is.character(palette) && length(palette) == 1L &&
             identical(palette, "capivara")) {
    plot <- plot + scale_capivara_identity()
  } else if (is.character(palette) && length(palette) == 1L &&
             identical(palette, "starry_night")) {
    plot <- plot + ggplot2::scale_fill_manual(
      values = .starry_night_palette(nrow(regions)),
      na.value = .capivara_invalid
    )
  } else if (is.character(palette) && length(palette) == 1L) {
    plot <- plot + viridis::scale_fill_viridis(discrete = TRUE, option = palette)
  } else if (is.character(palette) && length(palette) > 1L) {
    colours <- if (length(palette) < nrow(regions)) {
      grDevices::colorRampPalette(palette)(nrow(regions))
    } else palette[seq_len(nrow(regions))]
    plot <- plot + ggplot2::scale_fill_manual(values = colours,
                                               na.value = .capivara_invalid)
  } else {
    stop("`palette` must be a CAPIVARA/viridis name or colour vector.",
         call. = FALSE)
  }

  labels <- switch(coordinate_mode,
    native = list(x = "native column", y = "native row"),
    sky = list(x = expression(Delta*E~"[arcsec]"),
               y = expression(Delta*N~"[arcsec]")),
    physical = list(x = expression(Delta*E~"[kpc]"),
                    y = expression(Delta*N~"[kpc]"))
  )
  plot <- plot +
    ggplot2::coord_sf(datum = NA, expand = FALSE) +
    ggplot2::labs(x = labels$x, y = labels$y, fill = "region") +
    theme_capivara("spatial") +
    ggplot2::theme(legend.position = if (show_legend) "right" else "none")
  if (display_reflections == 1L) {
    plot <- plot + ggplot2::scale_x_continuous(labels = function(x) -x)
  }
  attr(plot, "capivara_provenance") <- list(
    coordinate_contract = contract,
    canonical_coordinates = if (coordinate_mode == "native")
      c(x = "column", y = "row") else c(x = "+East", y = "+North"),
    coordinate_mode = coordinate_mode,
    east_left_reflections = display_reflections,
    scientific_geometry_unchanged = identical(
      canonical_geometry,
      serialize(sf::st_geometry(regions_sf(cluster_data, coordinate_mode)),
                NULL, version = 3)
    ),
    renderer = "exact sf cell-union geometry"
  )
  plot
}
