test_that("frozen MaNGA atlas keeps native counts and incomplete apertures explicit", {
  skip_if_not_installed("jsonlite")
  path <- system.file("atlas", "objects", "8602-12705.json",
                      package = "capivara", mustWork = TRUE)
  d <- jsonlite::read_json(path, simplifyVector = FALSE)
  expect_identical(d$schema, "capivara-semantic-atlas-v1")
  expect_identical(d$object, "8602-12705")
  expect_match(d$source, "Observed")
  expect_identical(d$wavelength_frame, "observed")
  expect_identical(d$wavelength_medium, "vacuum")
  expect_length(d$wavelength, 2201L)
  expect_equal(c(d$wavelength[[1L]], tail(d$wavelength, 1L)[[1L]]),
               c(4800, 7000))
  expect_length(d$labels, 74L)
  expect_true(all(lengths(d$labels) == 74L))
  expect_length(d$regions, 8L)
  ids <- vapply(d$regions, function(r) r$region, integer(1))
  sizes <- vapply(d$regions, function(r) r$n_spaxels, integer(1))
  native <- unlist(d$labels, use.names = FALSE)
  expect_equal(unname(tabulate(native[!is.na(native)], nbins = 8L)), sizes)
  expect_equal(sum(sizes), 2761L)
  for (r in d$regions) {
    expect_length(r$spectrum, 2201L)
    expect_length(r$variance_spatial, 2201L)
    expect_length(r$variance_diagonal, 2201L)
    counts <- unlist(r$valid_counts, use.names = FALSE)
    complete <- counts == r$n_spaxels
    expect_true(all(vapply(r$spectrum[!complete], is.null, logical(1))))
    expect_true(all(vapply(r$variance_spatial[!complete], is.null, logical(1))))
    expect_true(all(vapply(r$variance_diagonal[!complete], is.null, logical(1))))
    expect_true(all(is.finite(unlist(r$spectrum[complete]))))
  }
  expect_match(d$provenance$uncertainty, "spectral covariance not reconstructed")
  expect_match(d$provenance$interpretation, "no new physical fit")
})

test_that("published MaNGA collection resolves to measured regional products", {
  skip_if_not_installed("jsonlite")
  root <- system.file("atlas", "objects", package = "capivara", mustWork = TRUE)
  index <- jsonlite::read_json(file.path(root, "collection.json"),
                               simplifyVector = FALSE)
  expect_identical(index$schema, "capivara-observed-collection-v1")
  expect_length(index$entries, 15L)
  ids <- vapply(index$entries, function(e) e$object, character(1))
  expect_false(anyDuplicated(ids) > 0L)
  expect_true("8602-12705" %in% ids)
  for (entry in index$entries) {
    object <- entry$object
    expect_true(file.exists(file.path(root, paste0(object, "-regions.svg"))))
    expect_true(file.exists(file.path(root, paste0(object, "-region_features.csv"))))
    expect_true(file.exists(file.path(root, paste0(object, "-region_adjacency.csv"))))
    d <- jsonlite::read_json(file.path(root, paste0(object, ".json")),
                             simplifyVector = FALSE)
    expect_identical(d$schema, "capivara-semantic-atlas-v1")
    expect_identical(d$object, object)
    expect_identical(d$representation, "full_flux")
    expect_identical(d$wavelength_medium, "vacuum")
    expect_length(d$wavelength, entry$n_channels)
    expect_length(d$regions, entry$n_regions)
    expect_equal(sum(vapply(d$regions, function(r) r$n_spaxels, integer(1))),
                 entry$n_spaxels)
    expect_match(d$provenance$interpretation, "no new physical fit")
    expect_match(d$provenance$independent_check, "Every selected native-channel")
  }
})
