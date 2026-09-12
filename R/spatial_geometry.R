# Exact spatial geometry for CAPIVARA partitions.
#
# Scientific arrays remain matrices indexed [row, column]. Cell corners are
# transformed explicitly; no transpose, rotation, smoothing, interpolation, or
# display reflection enters the stored geometry.

.capivara_default_cosmology <- function() {
  list(name = "flat_LCDM_H0_70_Om0_0.3", H0_km_s_Mpc = 70,
       Omega_m = 0.3, Omega_lambda = 0.7, Omega_k = 0)
}

.capivara_validate_cosmology <- function(cosmology) {
  required <- c("name", "H0_km_s_Mpc", "Omega_m", "Omega_lambda")
  if (!is.list(cosmology) || !all(required %in% names(cosmology)) ||
      !is.character(cosmology$name) || length(cosmology$name) != 1L ||
      !nzchar(cosmology$name)) {
    stop("`cosmology` must declare name, H0_km_s_Mpc, Omega_m, and Omega_lambda.",
         call. = FALSE)
  }
  values <- unlist(cosmology[c("H0_km_s_Mpc", "Omega_m", "Omega_lambda")])
  if (any(!is.finite(values)) || cosmology$H0_km_s_Mpc <= 0 ||
      cosmology$Omega_m < 0 || cosmology$Omega_lambda < 0) {
    stop("Cosmological parameters must be finite and physically admissible.",
         call. = FALSE)
  }
  cosmology$Omega_k <- 1 - cosmology$Omega_m - cosmology$Omega_lambda
  cosmology
}

.capivara_angular_diameter_distance_mpc <- function(redshift, cosmology) {
  if (!is.finite(redshift) || redshift <= 0) return(NA_real_)
  ez <- function(z) sqrt(cosmology$Omega_m * (1 + z)^3 +
                           cosmology$Omega_k * (1 + z)^2 +
                           cosmology$Omega_lambda)
  c_km_s <- 299792.458
  chi <- (c_km_s / cosmology$H0_km_s_Mpc) *
    stats::integrate(function(z) 1 / ez(z), 0, redshift,
                     rel.tol = 1e-10)$value
  ok <- cosmology$Omega_k
  dm <- if (abs(ok) < 1e-12) {
    chi
  } else if (ok > 0) {
    (c_km_s / cosmology$H0_km_s_Mpc) / sqrt(ok) *
      sinh(sqrt(ok) * chi * cosmology$H0_km_s_Mpc / c_km_s)
  } else {
    (c_km_s / cosmology$H0_km_s_Mpc) / sqrt(-ok) *
      sin(sqrt(-ok) * chi * cosmology$H0_km_s_Mpc / c_km_s)
  }
  dm / (1 + redshift)
}

#' Declare the spatial coordinate contract for a cube
#'
#' The native matrix convention is always `[row, column]`. A celestial WCS is
#' represented by its local two-dimensional CD matrix and an explicit mapping
#' from native row/column offsets to FITS axes. Canonical celestial coordinates
#' increase toward East and North; East-left is applied only by the renderer.
#'
#' @param centre Named numeric vector or list containing native `row` and
#'   `column` coordinates of the adopted centre.
#' @param wcs Optional list containing `cd1_1_deg`, `cd1_2_deg`, `cd2_1_deg`,
#'   and `cd2_2_deg`.
#' @param native_wcs_axes Named character vector mapping FITS `axis1` and
#'   `axis2` to native `row` and `column`.
#' @param redshift Object redshift. A finite positive value permits physical
#'   coordinates when WCS is also supplied.
#' @param redshift_source Source of the adopted redshift.
#' @param cosmology Named cosmology parameters.
#' @param provenance Additional named provenance retained verbatim.
#' @return A `capivara_spatial_contract`.
#' @export
capivara_spatial_contract <- function(
    centre, wcs = NULL,
    native_wcs_axes = c(axis1 = "column", axis2 = "row"),
    redshift = NA_real_, redshift_source = "not supplied",
    cosmology = .capivara_default_cosmology(), provenance = list()) {
  if (is.list(centre)) centre <- unlist(centre[c("row", "column")])
  if (!is.numeric(centre) || length(centre) != 2L ||
      any(!is.finite(centre))) {
    stop("`centre` must contain finite native row and column coordinates.",
         call. = FALSE)
  }
  if (is.null(names(centre))) names(centre) <- c("row", "column")
  if (!all(c("row", "column") %in% names(centre))) {
    stop("`centre` must be named with row and column.", call. = FALSE)
  }
  if (!all(c("axis1", "axis2") %in% names(native_wcs_axes)) ||
      !setequal(unname(native_wcs_axes), c("row", "column"))) {
    stop("`native_wcs_axes` must map axis1/axis2 to row/column.", call. = FALSE)
  }
  if (!is.null(wcs)) {
    required <- c("cd1_1_deg", "cd1_2_deg", "cd2_1_deg", "cd2_2_deg")
    if (!is.list(wcs) || !all(required %in% names(wcs)) ||
        any(!is.finite(unlist(wcs[required])))) {
      stop("`wcs` must contain a finite two-dimensional celestial CD matrix.",
           call. = FALSE)
    }
  }
  if (!is.numeric(redshift) || length(redshift) != 1L ||
      (!is.na(redshift) && (!is.finite(redshift) || redshift < 0))) {
    stop("`redshift` must be NA or one finite non-negative value.", call. = FALSE)
  }
  if (!is.character(redshift_source) || length(redshift_source) != 1L ||
      is.na(redshift_source) || !nzchar(redshift_source)) {
    stop("`redshift_source` must be one non-empty string.", call. = FALSE)
  }
  cosmology <- .capivara_validate_cosmology(cosmology)
  distance <- .capivara_angular_diameter_distance_mpc(redshift, cosmology)
  kpc_per_arcsec <- distance * 1000 / 206264.80624709636
  pixel_scale_arcsec <- if (is.null(wcs)) NA_real_ else {
    cd <- matrix(c(wcs$cd1_1_deg, wcs$cd1_2_deg,
                   wcs$cd2_1_deg, wcs$cd2_2_deg), nrow = 2L, byrow = TRUE)
    sqrt(abs(det(cd))) * 3600
  }
  out <- list(
    version = "capivara_spatial_contract_v1",
    native_indexing = "matrix[row, column]",
    centre = list(row = unname(centre[["row"]]),
                  column = unname(centre[["column"]])),
    wcs = wcs, native_wcs_axes = as.list(native_wcs_axes),
    redshift = redshift, redshift_source = redshift_source,
    cosmology = cosmology,
    angular_diameter_distance_mpc = distance,
    kpc_per_arcsec = kpc_per_arcsec,
    pixel_scale_arcsec = pixel_scale_arcsec,
    kpc_per_pixel = pixel_scale_arcsec * kpc_per_arcsec,
    coordinate_centre = list(row = unname(centre[["row"]]),
                             column = unname(centre[["column"]])),
    canonical_sky_coordinates = c(x = "Delta East", y = "Delta North"),
    display_reflections = 0L,
    provenance = provenance
  )
  class(out) <- c("capivara_spatial_contract", "list")
  out$contract_id <- .capivara_stable_id("spatial", out)
  out
}

.capivara_native_spatial_contract <- function(dimensions) {
  capivara_spatial_contract(
    centre = c(row = (dimensions[1L] + 1) / 2,
               column = (dimensions[2L] + 1) / 2),
    provenance = list(source = "native matrix dimensions")
  )
}

.capivara_spatial_mode <- function(contract, coords = "auto") {
  coords <- match.arg(coords, c("auto", "native", "sky", "physical"))
  has_wcs <- !is.null(contract$wcs)
  has_physical <- has_wcs && is.finite(contract$kpc_per_arcsec) &&
    contract$kpc_per_arcsec > 0
  if (coords == "auto") {
    if (has_physical) return("physical")
    if (has_wcs) return("sky")
    return("native")
  }
  if (coords == "sky" && !has_wcs) {
    stop("Sky coordinates require a declared celestial WCS.", call. = FALSE)
  }
  if (coords == "physical" && !has_physical) {
    stop("Physical coordinates require WCS and a finite positive redshift.",
         call. = FALSE)
  }
  coords
}

.capivara_wcs_axis_matrix <- function(contract) {
  axes <- unlist(contract$native_wcs_axes)
  rbind(
    axis1 = c(row = as.numeric(axes[["axis1"]] == "row"),
              column = as.numeric(axes[["axis1"]] == "column")),
    axis2 = c(row = as.numeric(axes[["axis2"]] == "row"),
              column = as.numeric(axes[["axis2"]] == "column"))
  )
}

.capivara_transform_spatial <- function(native_row, native_column, contract,
                                         coords) {
  coords <- .capivara_spatial_mode(contract, coords)
  if (coords == "native") {
    return(data.frame(x = native_column, y = native_row))
  }
  delta <- rbind(row = native_row - contract$centre$row,
                 column = native_column - contract$centre$column)
  fits_delta <- .capivara_wcs_axis_matrix(contract) %*% delta
  wcs <- contract$wcs
  cd <- matrix(c(wcs$cd1_1_deg, wcs$cd1_2_deg,
                 wcs$cd2_1_deg, wcs$cd2_2_deg), nrow = 2L, byrow = TRUE)
  tangent <- cd %*% fits_delta
  out <- data.frame(x = tangent[1L, ] * 3600,
                    y = tangent[2L, ] * 3600)
  if (coords == "physical") {
    out$x <- out$x * contract$kpc_per_arcsec
    out$y <- out$y * contract$kpc_per_arcsec
  }
  out
}

.capivara_native_table <- function(values) {
  data.frame(
    native_row = rep(seq_len(nrow(values)), times = ncol(values)),
    native_column = rep(seq_len(ncol(values)), each = nrow(values)),
    value = as.vector(values)
  )
}

.capivara_cell_sf <- function(dimensions, contract, coords) {
  template <- matrix(0, dimensions[1L], dimensions[2L])
  native <- .capivara_native_table(template)
  offsets <- rbind(c(-0.5, -0.5), c(-0.5, 0.5), c(0.5, 0.5),
                   c(0.5, -0.5), c(-0.5, -0.5))
  geometry <- lapply(seq_len(nrow(native)), function(index) {
    corners <- .capivara_transform_spatial(
      native$native_row[index] + offsets[, 1L],
      native$native_column[index] + offsets[, 2L], contract, coords
    )
    ring <- as.matrix(corners[c("x", "y")])
    ring[nrow(ring), ] <- ring[1L, ]
    sf::st_polygon(list(ring))
  })
  out <- sf::st_sf(
    cell_id = seq_len(nrow(native)),
    native_row = native$native_row,
    native_column = native$native_column,
    geometry = sf::st_sfc(geometry, crs = sf::NA_crs_)
  )
  if (any(!sf::st_is_valid(out))) stop("Spatial cell footprints are invalid.")
  attr(out, "coordinate_contract") <- contract
  attr(out, "coordinate_mode") <- coords
  out
}

.capivara_polygon_parts <- function(geometry) {
  type <- as.character(sf::st_geometry_type(geometry))
  raw <- sf::st_geometry(geometry)[[1L]]
  if (type == "POLYGON") {
    c(components = 1L, holes = length(raw) - 1L)
  } else if (type == "MULTIPOLYGON") {
    c(components = length(raw),
      holes = sum(vapply(raw, length, integer(1L)) - 1L))
  } else {
    stop("Exact cell unions must be POLYGON or MULTIPOLYGON.", call. = FALSE)
  }
}

.capivara_geometry_metadata <- function(x) {
  representation <- x$representation
  contract <- tryCatch(.capivara_object_contract(x), error = function(e) NULL)
  spatial_provenance <- if (is.null(contract)) list() else contract$provenance
  version <- tryCatch(as.character(getNamespaceVersion("capivara")),
                      error = function(e) NA_character_)
  list(
    galaxy_id = if (is.null(spatial_provenance$galaxy_id)) NA_character_ else
      as.character(spatial_provenance$galaxy_id),
    representation = if (is.list(representation) &&
      !is.null(representation$name)) representation$name else NA_character_,
    support_source = if (is.null(x$support_source)) "historical_default" else
      x$support_source,
    validity_contract = if (is.null(x$validity_contract)) "historical" else
      "explicit representation validity",
    eligibility_contract = if (is.null(x$eligibility_contract)) "historical" else
      "explicit representation eligibility",
    representation_id = if (is.null(x$representation_id)) NA_character_ else
      x$representation_id,
    support_id = if (is.null(x$support_id)) NA_character_ else x$support_id,
    capivara_version = version,
    capivara_commit = if (is.null(spatial_provenance$capivara_commit))
      NA_character_ else as.character(spatial_provenance$capivara_commit)
  )
}

.capivara_labels_to_sf <- function(labels, cells, metadata, coords) {
  native <- .capivara_native_table(labels)
  ids <- sort(unique(stats::na.omit(native$value)))
  if (!length(ids)) stop("The label matrix contains no regions.", call. = FALSE)
  rows <- lapply(ids, function(id) {
    members <- which(is.finite(native$value) & native$value == id)
    geometry <- sf::st_union(sf::st_geometry(cells)[members])
    if (!as.character(sf::st_geometry_type(geometry)) %in%
        c("POLYGON", "MULTIPOLYGON")) {
      stop("A cell union did not produce polygonal geometry.", call. = FALSE)
    }
    centroid <- suppressWarnings(sf::st_centroid(geometry))
    centre <- sf::st_coordinates(centroid)[1L, ]
    parts <- .capivara_polygon_parts(sf::st_sf(geometry = geometry))
    list(region_id = as.integer(id), n_spaxels = length(members),
         area = as.numeric(sf::st_area(geometry)),
         centroid_x = unname(centre[["X"]]),
         centroid_y = unname(centre[["Y"]]),
         n_components = unname(parts[["components"]]),
         n_holes = unname(parts[["holes"]]), geometry = geometry[[1L]])
  })
  attributes <- data.frame(
    region_id = vapply(rows, `[[`, integer(1L), "region_id"),
    n_spaxels = vapply(rows, `[[`, integer(1L), "n_spaxels"),
    area = vapply(rows, `[[`, numeric(1L), "area"),
    centroid_x = vapply(rows, `[[`, numeric(1L), "centroid_x"),
    centroid_y = vapply(rows, `[[`, numeric(1L), "centroid_y"),
    n_components = vapply(rows, `[[`, integer(1L), "n_components"),
    n_holes = vapply(rows, `[[`, integer(1L), "n_holes"),
    stringsAsFactors = FALSE
  )
  for (name in names(metadata)) attributes[[name]] <- metadata[[name]]
  if (coords == "physical") {
    attributes$area_kpc2 <- attributes$area
    attributes$centroid_E_kpc <- attributes$centroid_x
    attributes$centroid_N_kpc <- attributes$centroid_y
  } else {
    attributes$area_kpc2 <- NA_real_
    attributes$centroid_E_kpc <- NA_real_
    attributes$centroid_N_kpc <- NA_real_
  }
  out <- sf::st_sf(
    attributes,
    geometry = sf::st_sfc(lapply(rows, `[[`, "geometry"), crs = sf::NA_crs_)
  )
  if (any(sf::st_is_empty(out)) || any(!sf::st_is_valid(out))) {
    stop("Exact region geometry is empty or invalid.", call. = FALSE)
  }
  attr(out, "coordinate_mode") <- coords
  attr(out, "scientific_coordinates") <- switch(
    coords, native = c(x = "column", y = "row"),
    sky = c(x = "Delta East [arcsec]", y = "Delta North [arcsec]"),
    physical = c(x = "Delta East [kpc]", y = "Delta North [kpc]")
  )
  attr(out, "display_reflections") <- 0L
  out
}

.capivara_mask_to_sf <- function(mask, cells, metadata, coords, domain) {
  native <- .capivara_native_table(mask * 1L)
  members <- which(native$value == 1L)
  if (!length(members)) {
    return(sf::st_sf(domain = domain, n_spaxels = 0L, area = 0,
                     geometry = sf::st_sfc(sf::st_geometrycollection(),
                                           crs = sf::NA_crs_)))
  }
  geometry <- sf::st_union(sf::st_geometry(cells)[members])
  centroid <- suppressWarnings(sf::st_centroid(geometry))
  centre <- sf::st_coordinates(centroid)[1L, ]
  parts <- .capivara_polygon_parts(sf::st_sf(geometry = geometry))
  attributes <- data.frame(
    domain = domain, n_spaxels = length(members),
    area = as.numeric(sf::st_area(geometry)),
    centroid_x = unname(centre[["X"]]), centroid_y = unname(centre[["Y"]]),
    n_components = unname(parts[["components"]]),
    n_holes = unname(parts[["holes"]]), stringsAsFactors = FALSE
  )
  for (name in names(metadata)) attributes[[name]] <- metadata[[name]]
  if (coords == "physical") {
    attributes$area_kpc2 <- attributes$area
    attributes$centroid_E_kpc <- attributes$centroid_x
    attributes$centroid_N_kpc <- attributes$centroid_y
  }
  out <- sf::st_sf(attributes, geometry = geometry)
  attr(out, "coordinate_mode") <- coords
  attr(out, "display_reflections") <- 0L
  out
}

.capivara_native_adjacency <- function(labels) {
  append_pairs <- function(first, second) {
    keep <- is.finite(first) & is.finite(second) & first != second
    if (!any(keep)) return(matrix(integer(), ncol = 2L))
    unique(cbind(pmin(first[keep], second[keep]),
                 pmax(first[keep], second[keep])))
  }
  pairs <- matrix(integer(), ncol = 2L)
  if (nrow(labels) > 1L) pairs <- rbind(
    pairs, append_pairs(labels[-nrow(labels), , drop = FALSE],
                        labels[-1L, , drop = FALSE])
  )
  if (ncol(labels) > 1L) pairs <- rbind(
    pairs, append_pairs(labels[, -ncol(labels), drop = FALSE],
                        labels[, -1L, drop = FALSE])
  )
  if (!nrow(pairs)) return(data.frame(region_id_1 = integer(),
                                      region_id_2 = integer()))
  pairs <- unique(pairs)
  pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  data.frame(region_id_1 = as.integer(pairs[, 1L]),
             region_id_2 = as.integer(pairs[, 2L]))
}

.capivara_region_adjacency <- function(regions, labels) {
  tolerance <- sqrt(min(regions$area / regions$n_spaxels)) * 1e-8
  touches <- sf::st_touches(regions)
  pairs <- do.call(rbind, lapply(seq_along(touches), function(first) {
    second <- touches[[first]]
    second <- second[second > first]
    if (!length(second)) return(NULL)
    cbind(first = first, second = second)
  }))
  if (is.null(pairs)) pairs <- matrix(integer(), ncol = 2L)
  out <- if (nrow(pairs)) do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    first <- pairs[i, "first"]; second <- pairs[i, "second"]
    common <- suppressWarnings(sf::st_intersection(
      sf::st_boundary(sf::st_geometry(regions)[first]),
      sf::st_boundary(sf::st_geometry(regions)[second])
    ))
    data.frame(
      region_id_1 = regions$region_id[first],
      region_id_2 = regions$region_id[second], touches = TRUE,
      shared_boundary_length = sum(as.numeric(sf::st_length(common))),
      stringsAsFactors = FALSE
    )
  })) else data.frame(region_id_1 = integer(), region_id_2 = integer(),
                      touches = logical(), shared_boundary_length = numeric())
  out$adjacent <- out$shared_boundary_length > tolerance
  native <- .capivara_native_adjacency(labels)
  key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = ":")
  native_keys <- key(native$region_id_1, native$region_id_2)
  out$native_four_neighbour <- key(out$region_id_1, out$region_id_2) %in%
    native_keys
  out$agreement <- out$adjacent == out$native_four_neighbour
  polygon_keys <- key(out$region_id_1[out$adjacent],
                      out$region_id_2[out$adjacent])
  native_match <- setequal(polygon_keys, native_keys)
  graph <- list(
    directed = FALSE,
    vertices = regions$region_id,
    edges = out[out$adjacent,
                c("region_id_1", "region_id_2", "shared_boundary_length"),
                drop = FALSE]
  )
  class(graph) <- c("capivara_adjacency_graph", "list")
  attr(out, "native_adjacency") <- native
  attr(out, "native_adjacency_match") <- native_match
  attr(out, "graph") <- graph
  attr(out, "length_tolerance") <- tolerance
  out
}

.capivara_rasterize_regions <- function(regions, cells, dimensions) {
  centres <- sf::st_centroid(sf::st_geometry(cells))
  membership <- sf::st_intersects(centres, regions)
  if (any(lengths(membership) > 1L)) {
    stop("A cell centre belongs to more than one region.", call. = FALSE)
  }
  reconstructed <- rep(NA_integer_, nrow(cells))
  assigned <- lengths(membership) == 1L
  reconstructed[assigned] <- regions$region_id[
    vapply(membership[assigned], `[[`, integer(1L), 1L)
  ]
  matrix(reconstructed, dimensions[1L], dimensions[2L])
}

.capivara_object_contract <- function(x) {
  contract <- if (!is.null(x$spatial$contract)) x$spatial$contract else
    if (!is.null(x$original_cube$spatial_contract))
      x$original_cube$spatial_contract else NULL
  if (is.null(contract)) contract <- .capivara_native_spatial_contract(
    dim(x$cluster_map)
  )
  if (!inherits(contract, "capivara_spatial_contract")) {
    stop("The spatial contract is not a `capivara_spatial_contract`.",
         call. = FALSE)
  }
  contract
}

.capivara_domain_mask <- function(x, domain) {
  switch(domain,
    analysis = is.finite(x$cluster_map),
    validity = if (!is.null(x$representation_validity))
      x$representation_validity else is.finite(x$cluster_map),
    eligibility = if (!is.null(x$representation_eligibility))
      x$representation_eligibility else is.finite(x$cluster_map)
  )
}

#' Exact polygon geometry for CAPIVARA regions
#'
#' @param x A CAPIVARA segmentation result.
#' @param coords Coordinate mode: `"native"`, tangent-plane `"sky"` in
#'   arcsec, `"physical"` in kpc, or `"auto"`.
#' @return One `sf` POLYGON/MULTIPOLYGON feature per region.
#' @export
regions_sf <- function(x, coords = c("auto", "native", "sky", "physical")) {
  if (!is.list(x) || !is.matrix(x$cluster_map)) {
    stop("`x` must be a CAPIVARA segmentation result.", call. = FALSE)
  }
  coords <- match.arg(coords)
  contract <- .capivara_object_contract(x)
  mode <- .capivara_spatial_mode(contract, coords)
  if (!is.null(x$spatial$geometry_exact) &&
      identical(x$spatial$canonical_coords, mode)) {
    return(x$spatial$geometry_exact)
  }
  cells <- .capivara_cell_sf(dim(x$cluster_map), contract, mode)
  metadata <- .capivara_geometry_metadata(x)
  metadata$K <- length(unique(stats::na.omit(as.vector(x$cluster_map))))
  metadata$coordinate_contract <- contract$contract_id
  .capivara_labels_to_sf(x$cluster_map, cells, metadata, mode)
}

#' Exact polygon geometry for a CAPIVARA analysis domain
#'
#' @param x A CAPIVARA segmentation result.
#' @param domain One of `"analysis"`, `"validity"`, or `"eligibility"`.
#' @param coords Coordinate mode; see [regions_sf()].
#' @return A one-feature `sf` domain object.
#' @export
domain_sf <- function(x, domain = c("analysis", "validity", "eligibility"),
                      coords = c("auto", "native", "sky", "physical")) {
  domain <- match.arg(domain); coords <- match.arg(coords)
  contract <- .capivara_object_contract(x)
  mode <- .capivara_spatial_mode(contract, coords)
  if (!is.null(x$spatial$domains_exact[[domain]]) &&
      identical(x$spatial$canonical_coords, mode)) {
    return(x$spatial$domains_exact[[domain]])
  }
  cells <- .capivara_cell_sf(dim(x$cluster_map), contract, mode)
  metadata <- .capivara_geometry_metadata(x)
  metadata$coordinate_contract <- contract$contract_id
  .capivara_mask_to_sf(.capivara_domain_mask(x, domain), cells,
                       metadata, mode, domain)
}

#' Region adjacency from exact polygon topology
#'
#' @param x A CAPIVARA segmentation result.
#' @param coords Coordinate mode; see [regions_sf()].
#' @return A data frame of topological touches and positive shared-boundary
#'   adjacency, checked against native four-neighbour labels. Its `graph`
#'   attribute contains vertices and the positive shared-boundary edge table.
#' @export
region_adjacency <- function(x, coords = c("auto", "native", "sky", "physical")) {
  coords <- match.arg(coords)
  contract <- .capivara_object_contract(x)
  mode <- .capivara_spatial_mode(contract, coords)
  if (!is.null(x$spatial$adjacency) &&
      identical(x$spatial$canonical_coords, mode)) return(x$spatial$adjacency)
  .capivara_region_adjacency(regions_sf(x, mode), x$cluster_map)
}

.capivara_attach_spatial_products <- function(out, raw) {
  if (is.null(out$cluster_map)) return(out)
  contract <- raw$spatial_contract
  if (is.null(contract)) contract <- .capivara_native_spatial_contract(
    dim(out$cluster_map)
  )
  if (!inherits(contract, "capivara_spatial_contract")) {
    stop("`input$spatial_contract` must be made by capivara_spatial_contract().",
         call. = FALSE)
  }
  mode <- .capivara_spatial_mode(contract, "auto")
  cells <- .capivara_cell_sf(dim(out$cluster_map), contract, mode)
  metadata <- .capivara_geometry_metadata(out)
  metadata$K <- length(unique(stats::na.omit(as.vector(out$cluster_map))))
  metadata$coordinate_contract <- contract$contract_id
  geometry <- .capivara_labels_to_sf(out$cluster_map, cells, metadata, mode)
  domains <- lapply(c("analysis", "validity", "eligibility"), function(domain) {
    .capivara_mask_to_sf(.capivara_domain_mask(out, domain), cells, metadata,
                         mode, domain)
  })
  names(domains) <- c("analysis", "validity", "eligibility")
  adjacency <- .capivara_region_adjacency(geometry, out$cluster_map)
  if (any(!adjacency$agreement) ||
      !isTRUE(attr(adjacency, "native_adjacency_match"))) {
    stop("Polygon adjacency disagrees with native four-neighbour adjacency.",
         call. = FALSE)
  }
  if (!identical(.capivara_rasterize_regions(
    geometry, cells, dim(out$cluster_map)), out$cluster_map)) {
    stop("Exact region geometry failed the labelled-cell round trip.",
         call. = FALSE)
  }
  out$spatial <- list(
    contract = contract, canonical_coords = mode,
    geometry_exact = geometry, domains_exact = domains,
    adjacency = adjacency, adjacency_graph = attr(adjacency, "graph"),
    coordinate_chain = paste(
      "native [row,column] -> FITS celestial WCS ->",
      "tangent-plane East/North -> physical coordinates"
    ),
    scientific_geometry_smoothing = FALSE,
    display_east_left_applied = FALSE
  )
  out$regions_exact <- geometry
  out$domain_exact <- domains$analysis
  out$region_adjacency_exact <- adjacency
  class(out) <- unique(c("capivara_segmentation", class(out)))
  out
}
