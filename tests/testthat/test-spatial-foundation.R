spatial_fixture_result <- function(labels, wcs = NULL, redshift = NA_real_) {
  contract <- capivara_spatial_contract(
    centre = c(row = (nrow(labels) + 1) / 2,
               column = (ncol(labels) + 1) / 2),
    wcs = wcs, redshift = redshift, redshift_source = "test fixture",
    provenance = list(galaxy_id = "asymmetric-fixture",
                      capivara_commit = "fixture")
  )
  validity <- matrix(TRUE, nrow(labels), ncol(labels))
  eligibility <- validity
  eligibility[1L, 1L] <- FALSE
  structure(list(
    cluster_map = labels,
    original_cube = list(spatial_contract = contract),
    representation = list(name = "fixture"),
    representation_id = "representation-fixture",
    support_id = "support-fixture",
    support_source = "representation_domain_default",
    validity_contract = list(rule = "fixture validity"),
    eligibility_contract = list(rule = "fixture eligibility"),
    representation_validity = validity,
    representation_eligibility = eligibility
  ), class = "capivara_segmentation")
}

test_that("exact cell unions preserve components, holes, domains, and labels", {
  labels <- matrix(c(
    1, 1, 1, 1, 2, 2,
    1, 3, 3, 1, 2, 2,
    1, 3, 3, 1, 4, 4,
    1, 1, 1, 1, 4, 2
  ), nrow = 4L, byrow = TRUE)
  storage.mode(labels) <- "integer"
  result <- spatial_fixture_result(labels)
  regions <- regions_sf(result, "native")
  cells <- .capivara_cell_sf(dim(labels),
                             result$original_cube$spatial_contract, "native")

  expect_s3_class(regions, "sf")
  expect_identical(regions$region_id, 1:4)
  expect_identical(regions$n_spaxels, as.integer(tabulate(labels, nbins = 4L)))
  expect_identical(regions$n_holes[regions$region_id == 1L], 1L)
  expect_identical(regions$n_components[regions$region_id == 2L], 2L)
  expect_identical(
    as.character(sf::st_geometry_type(regions[regions$region_id == 2L, ])),
    "MULTIPOLYGON"
  )
  expect_true(all(sf::st_is_valid(regions)))
  expect_true(all(lengths(sf::st_overlaps(regions)) == 0L))
  expect_identical(.capivara_rasterize_regions(regions, cells, dim(labels)),
                   labels)

  valid <- domain_sf(result, "validity", "native")
  eligible <- domain_sf(result, "eligibility", "native")
  analysis <- domain_sf(result, "analysis", "native")
  expect_identical(valid$n_spaxels, length(labels))
  expect_identical(eligible$n_spaxels, length(labels) - 1L)
  expect_gt(as.numeric(sf::st_area(valid)), as.numeric(sf::st_area(eligible)))
  difference <- sf::st_sym_difference(sf::st_geometry(analysis),
                                      sf::st_union(regions))
  expect_equal(sum(as.numeric(sf::st_area(difference))), 0, tolerance = 1e-14)

  expected_area <- vapply(regions$region_id, function(region_id) {
    sum(as.numeric(sf::st_area(cells))[as.vector(labels) == region_id])
  }, numeric(1L))
  expect_equal(regions$area, expected_area, tolerance = 1e-12)
})

test_that("polygon topology reproduces native four-neighbour adjacency", {
  labels <- matrix(c(1L, 1L, 2L, 3L, 3L, 2L, 4L, 4L, 2L), 3L, 3L)
  result <- spatial_fixture_result(labels)
  adjacency <- region_adjacency(result, "native")
  topology <- adjacency[adjacency$adjacent,
                        c("region_id_1", "region_id_2"), drop = FALSE]
  native <- .capivara_native_adjacency(labels)
  rownames(topology) <- NULL
  rownames(native) <- NULL
  expect_true(all(adjacency$agreement))
  expect_true(attr(adjacency, "native_adjacency_match"))
  expect_s3_class(attr(adjacency, "graph"), "capivara_adjacency_graph")
  expect_identical(topology, native)
  expect_true(all(adjacency$shared_boundary_length[adjacency$adjacent] > 0))
})

test_that("physical geometry records scale, centre, and provenance", {
  labels <- matrix(c(1L, 1L, 2L, 2L), 2L, 2L)
  wcs <- list(cd1_1_deg = 1 / 3600, cd1_2_deg = 0,
              cd2_1_deg = 0, cd2_2_deg = 1 / 3600)
  result <- spatial_fixture_result(labels, wcs, redshift = 0.02)
  contract <- result$original_cube$spatial_contract
  regions <- regions_sf(result, "physical")
  expect_true(is.finite(contract$kpc_per_arcsec))
  expect_equal(contract$pixel_scale_arcsec, 1, tolerance = 1e-14)
  expect_equal(contract$kpc_per_pixel, contract$kpc_per_arcsec,
               tolerance = 1e-14)
  expect_true(all(is.finite(regions$area_kpc2)))
  expect_true(all(is.finite(regions$centroid_E_kpc)))
  expect_true(all(is.finite(regions$centroid_N_kpc)))
  expect_identical(unique(regions$galaxy_id), "asymmetric-fixture")
  expect_identical(unique(regions$capivara_commit), "fixture")
})

test_that("15-group asymmetric fixture applies East-left exactly once", {
  labels <- matrix(1L, nrow = 5L, ncol = 7L)
  labels[5L, 3L] <- 2L              # north
  labels[3L, 1L] <- 3L              # west
  labels[2L, 7L] <- 4L              # east
  labels[1L, 5L] <- 5L              # south
  extra <- rbind(
    c(1L, 1L), c(1L, 3L), c(2L, 2L), c(2L, 5L), c(3L, 3L),
    c(3L, 6L), c(4L, 2L), c(4L, 4L), c(4L, 7L), c(5L, 6L)
  )
  labels[extra] <- 6:15
  wcs <- list(cd1_1_deg = 1 / 3600, cd1_2_deg = 0,
              cd2_1_deg = 0, cd2_2_deg = 1 / 3600)
  result <- spatial_fixture_result(labels, wcs, redshift = 0.02)
  canonical <- regions_sf(result, "sky")
  before <- serialize(sf::st_geometry(canonical), NULL, version = 3)
  displayed <- .capivara_east_left_display_sf(canonical)

  east <- canonical[canonical$region_id == 4L, ]
  north <- canonical[canonical$region_id == 2L, ]
  east_display <- displayed[displayed$region_id == 4L, ]
  north_display <- displayed[displayed$region_id == 2L, ]
  expect_gt(sf::st_coordinates(sf::st_centroid(sf::st_geometry(east)))[1L, "X"], 0)
  expect_gt(sf::st_coordinates(sf::st_centroid(sf::st_geometry(north)))[1L, "Y"], 0)
  expect_lt(sf::st_coordinates(sf::st_centroid(sf::st_geometry(east_display)))[1L, "X"], 0)
  expect_gt(sf::st_coordinates(sf::st_centroid(sf::st_geometry(north_display)))[1L, "Y"], 0)

  plot <- plot_cluster(result, coords = "sky")
  provenance <- attr(plot, "capivara_provenance")
  expect_identical(provenance$east_left_reflections, 1L)
  expect_true(provenance$scientific_geometry_unchanged)
  expect_s3_class(ggplot2::ggplot_build(plot), "ggplot_built")
  expect_identical(serialize(sf::st_geometry(regions_sf(result, "sky")),
                             NULL, version = 3), before)

  cells <- .capivara_cell_sf(dim(labels),
                             result$original_cube$spatial_contract, "sky")
  reconstructed <- .capivara_rasterize_regions(canonical, cells, dim(labels))
  expect_identical(reconstructed, labels)
  rebuilt <- .capivara_labels_to_sf(
    reconstructed, cells, .capivara_geometry_metadata(result), "sky"
  )
  hausdorff <- as.numeric(sf::st_distance(
    sf::st_boundary(sf::st_union(canonical)),
    sf::st_boundary(sf::st_union(rebuilt)), which = "Hausdorff"
  ))
  expect_equal(hausdorff, 0, tolerance = 1e-14)
  expect_identical(nrow(canonical), 15L)
})

test_that("Van Gogh scales and theme are display-only ggplot components", {
  continuous <- scale_capivara_continuous()
  diverging <- scale_capivara_diverging()
  identity <- scale_capivara_identity()
  ordered <- scale_capivara_rank(8L)
  expect_s3_class(continuous, "ScaleContinuous")
  expect_s3_class(diverging, "ScaleContinuous")
  expect_s3_class(identity, "ScaleDiscrete")
  expect_s3_class(ordered, "ScaleContinuous")
  expect_s3_class(theme_capivara(), "theme")
  expect_identical(.capivara_diverging_palette(257L)[129L], "#F4F0DC")
})

test_that("arbitrary polygons intersect exact regions in the same coordinates", {
  labels <- matrix(rep(1:2, each = 8L), 4L, 4L)
  result <- spatial_fixture_result(labels)
  regions <- regions_sf(result, "native")
  supplied <- sf::st_sfc(sf::st_polygon(list(matrix(c(
    1.5, 1.5, 3.5, 1.5, 3.5, 3.5, 1.5, 3.5, 1.5, 1.5
  ), ncol = 2L, byrow = TRUE))), crs = sf::NA_crs_)
  intersection <- suppressWarnings(sf::st_intersection(regions, supplied))
  expect_s3_class(intersection, "sf")
  expect_true(nrow(intersection) > 0L)
  expect_gt(sum(as.numeric(sf::st_area(intersection))), 0)
})
