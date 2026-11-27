test_that("native bandpass recovers a known excess and its uncertainty", {
  wave <- 4980:5020
  flux <- 10 + 0.01 * (wave - 5000)
  flux[wave == 5000] <- flux[wave == 5000] + 5
  sidebands <- rbind(c(4980, 4984), c(5016, 5020))
  args <- list(wavelength = wave, flux = flux,
               variance = rep(4, length(wave)),
               line_window = c(4998.5, 5001.5),
               continuum_windows = sidebands,
               flux_unit = "test flux density",
               lsf_sigma_angstrom = 1.2)
  measured <- do.call(measure_spectral_bandpass, args)
  expect_identical(measured$status, "ok")
  expect_equal(measured$line_flux, 5, tolerance = 1e-12)
  expect_equal(measured$continuum_density, 10, tolerance = 1e-12)
  expect_equal(measured$equivalent_width, 0.5, tolerance = 1e-12)
  expect_equal(measured$measurement_covariance[1, 1], 15.6,
               tolerance = 1e-10)
  expect_equal(measured$measurement_covariance[2, 2], 0.4,
               tolerance = 1e-10)
  expect_equal(measured$measurement_covariance[1, 2], -1.2,
               tolerance = 1e-10)
  expect_equal(measured$equivalent_width_se^2, 0.169,
               tolerance = 1e-10)
  expect_equal(measured$instrumental_sigma_angstrom, 1.2)
  expect_equal(measured$window_width_over_lsf_fwhm,
               3 / (2 * sqrt(2 * log(2)) * 1.2))
  expect_match(measured$resolution_policy, "not deconvolved")

  rho <- 0.2
  C <- 4 * ((1 - rho) * diag(length(wave)) +
              rho * matrix(1, length(wave), length(wave)))
  measured_cov <- do.call(measure_spectral_bandpass,
                          c(args, list(spectral_covariance = C)))
  expect_equal(measured_cov$line_flux, measured$line_flux)
  expect_equal(measured_cov$measurement_covariance[1, 1], 12.48,
               tolerance = 1e-10)
  expect_equal(measured_cov$measurement_covariance[2, 2], 1.12,
               tolerance = 1e-10)
  expect_equal(measured_cov$measurement_covariance[1, 2], -0.96,
               tolerance = 1e-10)
  expect_equal(measured_cov$equivalent_width_se^2, 0.1372,
               tolerance = 1e-10)
})

test_that("missing support and invalid covariance cannot create a measurement", {
  wave <- 4980:5020
  flux <- rep(10, length(wave))
  args <- list(wavelength = wave, flux = flux,
               variance = rep(1, length(wave)),
               line_window = c(4998.5, 5001.5),
               continuum_windows = rbind(c(4980, 4984), c(5016, 5020)),
               flux_unit = "test flux density")
  args$flux[wave == 5000] <- NA_real_
  incomplete <- do.call(measure_spectral_bandpass, args)
  expect_identical(incomplete$status, "incomplete_support")
  expect_true(is.na(incomplete$line_flux))
  expect_true(is.na(incomplete$equivalent_width))
  args$flux[wave == 5000] <- 10
  C <- diag(length(wave))
  C[1, 1] <- 2
  expect_error(do.call(measure_spectral_bandpass,
                       c(args, list(spectral_covariance = C))),
               "match `variance`")
  args$lsf_sigma_angstrom <- 0
  expect_error(do.call(measure_spectral_bandpass, args),
               "positive")
})
