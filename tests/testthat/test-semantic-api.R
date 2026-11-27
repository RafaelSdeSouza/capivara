# Regression checks for additive interfaces. No simulated performance claims.
semantic_api_fixture <- function() {
  wave <- seq(4800, 7400, by = 20)
  base <- 2 + 0.0002 * (wave - 6000)
  flux <- array(outer(c(1, 1.1, 3, 3.1, 5, 5.1), base),
                c(2L, 3L, length(wave)))
  variance <- array(0.01, dim(flux))
  valid <- array(TRUE, dim(flux))
  valid[1, 2, 10] <- FALSE
  list(cube = list(imDat = flux, wavelength = wave,
                   wavelength_frame = "rest", wavelength_medium = "vacuum"),
       variance = variance, valid = valid)
}

test_that("canonical segmentation uses the same declared hierarchy", {
  d <- semantic_api_fixture()
  rep <- full_flux_representation()
  new <- segment_spectra(d$cube, 3L, rep, variance = d$variance,
                         sample_validity = d$valid)
  old <- segment_large(d$cube, Ncomp = 3L, representation = rep,
                        var_cube = d$variance, sample_validity = d$valid)
  expect_identical(new$cluster_map, old$cluster_map)
  expect_equal(new$representation, old$representation)
  expect_equal(new$hierarchy, old$hierarchy)
  expect_identical(new$original_cube$imDat, old$original_cube$imDat)
  expect_identical(new$coverage_qc, old$coverage_qc)
})

test_that("canonical segmentation requires scientific choices", {
  d <- semantic_api_fixture()
  expect_error(segment_spectra(d$cube, representation = full_flux_representation()),
               "positive integer n_regions")
  expect_error(segment_spectra(d$cube, 2), "explicit representation")
  expect_error(segment_spectra(d$cube, 2, NULL), "explicit representation")
  for (k in list(0, -1, 1.5, Inf, NA_real_, "2", c(2, 3))) {
    expect_error(segment_spectra(d$cube, k, full_flux_representation()),
                 "positive integer n_regions")
  }
})

test_that("regional extraction and exact spatial geometry preserve old results", {
  d <- semantic_api_fixture()
  s <- segment_spectra(d$cube, 3L, full_flux_representation(),
                       variance = d$variance, sample_validity = d$valid)
  s$original_cube$imDat[!d$valid] <- NA_real_
  expect_identical(extract_region_spectra(s, variance = d$variance),
                   summarize_cluster_spectra(s, var_cube = d$variance))
  expect_identical(extract_region_spectra(s, variance = d$variance,
                                         variance_inflation = 2),
                   summarize_cluster_spectra(s, var_cube = d$variance,
                                             variance_inflation = 2))
  expect_identical(region_geometry(s, "native"), regions_sf(s, "native"))
  expect_identical(region_geometry(s), regions_sf(s))
  a <- ggplot2::ggplot_build(plot_regions(s, coords = "native", show_legend = TRUE))
  b <- ggplot2::ggplot_build(plot_cluster(s, coords = "native", show_legend = TRUE))
  expect_equal(a$data, b$data)
})

test_that("optical line registry is identical in either wavelength convention", {
  expect_identical(list_emission_lines(), emission_lines())
  expect_identical(list_emission_lines("air"), emission_lines("air"))
})

test_that("legacy function names and representative defaults remain available", {
  legacy <- c("build_adaptive_support", "build_capivara_support",
              "build_quality_support", "build_starlet_mask",
              "capivara_amplitude_functional", "capivara_representation",
              "capivara_spatial_contract", "choose_ncomp_by_snr",
              "convert_wavelength_medium", "domain_sf", "emission_lines",
              "estimate_segment_memory", "full_flux_representation",
              "kinematic_models", "kinematic_panels",
              "manga_sandra_spectral_shape", "plot_cluster", "plot_cluster_spectra",
              "prepare_capivara_representation", "prepare_segment_spectra",
              "read_manga_lsf", "reconstruct_cluster_cube",
              "reconstruct_flux_preserving_cube", "region_adjacency", "regions_sf",
              "representation_coverage_qc", "run_kinematic_analysis",
              "run_manga_bar_model", "scale_capivara_continuous",
              "scale_capivara_diverging", "scale_capivara_identity",
              "scale_capivara_rank", "segment", "segment_emission_lines",
              "segment_kinematics", "segment_large", "select_manga_lsf",
              "spectral_shape_diagnostic", "spectral_shape_representation",
              "summarize_cluster_spectra", "theme_capivara")
  expect_true(all(legacy %in% getNamespaceExports("capivara")))
  expect_identical(formals(segment)$Ncomp, 15)
  expect_identical(formals(segment_large)$Ncomp, 15)
  expect_identical(formals(segment_emission_lines)$Ncomp, 50)
  expect_identical(formals(segment_kinematics)$n_segments, 25)
  expect_null(formals(segment)$representation)
  expect_null(formals(segment_large)$representation)
  expect_identical(formals(segment_large)$spatial_weight, 0)
  expect_identical(formals(summarize_cluster_spectra)$variance_inflation, 1)
})
