relaxed_full_flux <- function(min_shared = 0.7, min_contributor = 0.5) {
  full_flux_representation(required_support = list(
    min_valid_fraction = 0.7,
    min_reliable_shared_fraction = min_shared,
    min_cluster_contributor_fraction = min_contributor,
    feature_windows = list(),
    min_feature_fraction = 0
  ))
}

shape_fixture <- function(nr = 2L, nc = 2L, gains = rep(1, nr * nc),
                          amplitude_snr = rep(50, nr * nc)) {
  wave <- seq(4800, 7400, by = 10)
  base <- 2 + 0.0002 * (wave - 6000) -
    0.25 * exp(-0.5 * ((wave - 5175) / 10)^2) +
    0.15 * exp(-0.5 * ((wave - 6563) / 7)^2)
  flux_matrix <- matrix(NA_real_, nr * nc, length(wave))
  variance_matrix <- flux_matrix
  profile <- manga_sandra_spectral_shape()
  anchor <- Reduce(`|`, lapply(profile$required_support$amplitude_windows,
                               function(w) wave >= w[1] & wave <= w[2]))
  for (i in seq_len(nr * nc)) {
    spectrum <- gains[i] * base
    amplitude <- mean(spectrum[anchor])
    channel_sigma <- amplitude * sqrt(sum(anchor)) / amplitude_snr[i]
    flux_matrix[i, ] <- spectrum
    variance_matrix[i, ] <- channel_sigma^2
  }
  flux <- array(flux_matrix, c(nr, nc, length(wave)))
  variance <- array(variance_matrix, c(nr, nc, length(wave)))
  input <- list(imDat = flux, wavelength = wave,
                wavelength_frame = "rest", wavelength_medium = "vacuum")
  list(input = input, variance = list(imDat = variance),
       validity = array(TRUE, dim(flux)), profile = profile, anchor = anchor)
}

complete_support <- function(nr, nc, representation_eligibility = NULL) {
  q <- matrix(TRUE, nr, nc)
  build_capivara_support(
    q, q, q, representation_eligibility = representation_eligibility,
    construction_method = "deterministic semantic test fixture",
    source = "package regression", provenance = list(fixture = "semantic-v4")
  )
}

direct_observed_sse <- function(ids, x, q) {
  sum(vapply(seq_len(ncol(x)), function(j) {
    observed <- ids[q[ids, j]]
    if (!length(observed)) return(0)
    sum((x[observed, j] - mean(x[observed, j]))^2)
  }, numeric(1)))
}

test_that("semantic declarations expose scientific meaning and validation scope", {
  shape <- manga_sandra_spectral_shape()
  flux <- full_flux_representation()
  required <- c("name", "type", "coordinate_domain", "invariances", "retains",
                "removes", "required_support", "eligibility", "units",
                "validation_domain", "provenance", "representation_id")
  expect_true(all(required %in% names(shape)))
  expect_identical(shape$type, "spectral_shape")
  expect_identical(flux$type, "full_flux")
  expect_identical(shape$required_support$amplitude_windows,
                   list(c(5050, 5100), c(5400, 5500),
                        c(6000, 6100), c(6800, 6900)))
  expect_identical(shape$amplitude_functional$required_support$weighting,
                   "equal native samples over complete union")
  expect_equal(shape$eligibility$measured_amplitude_snr_min, 30)
  expect_match(shape$validation_domain$status, "VALIDATED_FOR_SN_GE_30")
  expect_match(flux$validation_domain$status, "NOT_VALIDATED_BY_EXPERIMENTS")
  expect_false("evaluator" %in% names(.capivara_declaration(shape)))
})

test_that("new semantic arguments preserve every historical positional slot", {
  segment_prefix <- c(
    "input", "Ncomp", "redshift", "scale_fn", "target_snr", "var_cube",
    "k_values", "wavelength_range", "feature_wavelength_range", "snr_stat",
    "variance_inflation", "use_starlet_mask", "support_method", "support_args",
    "collapse_fn", "starlet_J", "starlet_scales", "include_coarse", "denoise_k",
    "starlet_mode", "positive_only", "mask_mode", "wavelength_frame",
    "feature_wavelength_frame"
  )
  large_prefix <- c(
    segment_prefix[1:22], "knn_k", "auto_k", "max_k", "feature_scale",
    "spatial_weight", "mask", "valid_mode", "return_details", "verbose",
    "wavelength_frame", "feature_wavelength_frame"
  )
  expect_identical(head(names(formals(segment)), length(segment_prefix)),
                   segment_prefix)
  expect_identical(head(names(formals(segment_large)), length(large_prefix)),
                   large_prefix)
  expect_identical(tail(names(formals(segment)), 3),
                   c("support", "representation", "sample_validity"))
  expect_identical(tail(names(formals(segment_large)), 3),
                   c("support", "representation", "sample_validity"))
})

test_that("generic amplitude functionals support non-wavelength coordinates without validation claims", {
  evaluator <- function(flux, variance, validity, coordinate) {
    complete <- all(validity)
    list(value = if (complete) mean(flux) else NA_real_,
         variance = if (complete) sum(variance) / length(flux)^2 else NA_real_,
         complete = complete, samples = length(flux))
  }
  amplitude <- capivara_amplitude_functional(
    "frequency_mean", evaluator,
    coordinate_domain = list(kind = "frequency", frame = "native"),
    validation_domain = list(status = "UNVALIDATED", instrument = "example")
  )
  representation <- spectral_shape_representation(
    amplitude,
    coordinate_domain = list(kind = "frequency", frame = "native"),
    required_support = list(
      min_valid_fraction = 1, min_reliable_shared_fraction = 1,
      min_cluster_contributor_fraction = 1, feature_windows = list(),
      min_feature_fraction = 0
    ),
    eligibility = list(measured_amplitude_snr_min = 1),
    validation_domain = list(status = "UNVALIDATED")
  )
  input <- list(imDat = array(rbind(1:5, 3 * (1:5)), c(1, 2, 5)),
                coordinate = seq(100, 500, 100))
  p <- prepare_capivara_representation(
    input, representation, var_cube = array(0.01, c(1, 2, 5))
  )
  expect_equal(p$features[1, ], p$features[2, ], tolerance = 1e-15)
  expect_equal(p$amplitude[2], 3 * p$amplitude[1])
  expect_identical(representation$validation_domain$status, "UNVALIDATED")

  full <- prepare_capivara_representation(input, full_flux_representation())
  expect_equal(full$features[2, ], 3 * full$features[1, ])
})

test_that("spectral shape is invariant to positive amplitude and returns brightness", {
  a <- shape_fixture(gains = c(1, 2, 3, 4))
  b <- shape_fixture(gains = c(7, 14, 21, 28))
  pa <- prepare_capivara_representation(a$input, a$profile, a$variance,
                                        a$validity)
  pb <- prepare_capivara_representation(b$input, b$profile, b$variance,
                                        b$validity)
  expect_equal(pa$features, pb$features, tolerance = 2e-15)
  expect_equal(pb$amplitude, 7 * pa$amplitude, tolerance = 2e-15)
  expect_equal(pb$amplitude_snr, pa$amplitude_snr, tolerance = 2e-13)
  expect_true(all(pa$eligible) && all(pb$eligible))
})

test_that("amplitude eligibility is measured at S/N 30 and needs a complete union", {
  x <- shape_fixture(amplitude_snr = c(29.9, 30.1, 31, 50))
  p <- prepare_capivara_representation(x$input, x$profile, x$variance,
                                       x$validity)
  expect_identical(p$eligible, c(FALSE, TRUE, TRUE, TRUE))
  incomplete <- x$validity
  anchor_channel <- which(x$anchor)[1]
  incomplete[2, 1, anchor_channel] <- FALSE
  q <- prepare_capivara_representation(x$input, x$profile, x$variance,
                                       incomplete)
  expect_false(q$eligible[2])
  expect_true(q$eligibility_causes$incomplete_amplitude[2, 1])
})

test_that("explicit representations default to their own validity and eligibility domains", {
  x <- shape_fixture(amplitude_snr = c(29.9, 30.1, 31, 50))
  x$validity[1, 2, 20] <- FALSE
  prepared <- prepare_capivara_representation(
    x$input, x$profile, x$variance, x$validity
  )
  support <- .capivara_default_representation_support(prepared)
  analysis <- .capivara_analysis_support(support, prepared)
  expect_identical(support$source, "representation_domain_default")
  expect_false(support$configuration$starlet_intersection)
  expect_identical(support$analysis_mask, prepared$validity_map)
  expect_identical(analysis$final_analysis_mask, prepared$eligible_map)
  expect_identical(analysis$representation$name, "spectral_shape")
  expect_identical(analysis$validity_contract$missing_samples,
                   "retained as missing; never zero-filled")

  fit <- segment(
    x$input, Ncomp = 2L, var_cube = x$variance,
    representation = x$profile, sample_validity = x$validity
  )
  expect_identical(fit$support_source, "representation_domain_default")
  expect_false(fit$backend_info$starlet_intersection)
  expect_identical(fit$analysis_support$final_analysis_mask,
                   prepared$eligible_map)
  expect_true(all(c("support_source", "validity_contract",
                    "eligibility_contract", "representation") %in% names(fit)))
})

test_that("one bad non-anchor channel remains missing without rejecting its spaxel", {
  x <- shape_fixture()
  channel <- which(!x$anchor)[20]
  x$validity[1, 1, channel] <- FALSE
  p <- prepare_capivara_representation(x$input, x$profile, x$variance,
                                       x$validity)
  expect_true(p$eligible[1])
  selected_channel <- match(channel, p$selected_channel_indices)
  expect_true(is.na(p$features[1, selected_channel]))
  expect_false(p$sample_validity[1, selected_channel])
})

test_that("admissible shared wavelengths give exact missingness neutrality", {
  x <- rbind(seq_len(10), seq_len(10))
  q <- matrix(TRUE, 2, 10)
  q[1, 1] <- FALSE
  q[2, 10] <- FALSE
  stored <- x
  stored[!q] <- 0
  attr(stored, "coordinate") <- seq_len(10)
  attr(stored, "input_error_gamma") <- .capivara_gamma(1)
  h <- .observed_entry_hierarchy(
    stored, q, matrix(TRUE, 1, 2), 1:2,
    relaxed_full_flux(min_shared = 0.7), 1
  )
  expect_equal(h$cost, 0)
  expect_equal(h$actual_k, 1)
  expect_identical(h$qc$imputed_samples, 0L)
})

test_that("complete-data Ward agrees with Ward.D2 and direct SSE", {
  x <- rbind(c(0, 0), c(1, 1), c(10, 10), c(11, 11))
  q <- matrix(TRUE, nrow(x), ncol(x))
  attr(x, "coordinate") <- 1:2
  attr(x, "input_error_gamma") <- 0
  h <- .observed_entry_hierarchy(
    x, q, matrix(TRUE, 1, 4), 1:4,
    relaxed_full_flux(min_shared = 1, min_contributor = 1), 1
  )
  conventional <- stats::hclust(stats::dist(x), method = "ward.D2")
  expect_equal(h$cost, conventional$height^2 / 2, tolerance = 1e-13)
  expect_equal(h$cost, c(1, 1, 200), tolerance = 1e-13)
  expect_equal(h$labels, rep(1L, 4))
})

test_that("every observed-entry merge cost is the direct SSE increase", {
  x <- matrix(c(
    0, 1, 2, 3,
    1, 2, 3, 4,
    5, 4, 3, 2,
    6, 5, 4, 1
  ), nrow = 4, byrow = TRUE)
  q <- matrix(c(
    TRUE, TRUE, TRUE, FALSE,
    TRUE, TRUE, FALSE, TRUE,
    TRUE, FALSE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE
  ), nrow = 4, byrow = TRUE)
  stored <- x
  stored[!q] <- 0
  attr(stored, "coordinate") <- 1:4
  attr(stored, "input_error_gamma") <- 0
  h <- .observed_entry_hierarchy(
    stored, q, matrix(TRUE, 1, 4), 1:4,
    relaxed_full_flux(min_shared = 0.25, min_contributor = 0.25), 1
  )
  members <- lapply(seq_len(nrow(x)), function(i) i)
  for (i in seq_len(nrow(h$children))) {
    left <- members[[h$children[i, 1]]]
    right <- members[[h$children[i, 2]]]
    union <- c(left, right)
    direct <- direct_observed_sse(union, x, q) -
      direct_observed_sse(left, x, q) - direct_observed_sse(right, x, q)
    expect_equal(h$cost[i], direct, tolerance = 1e-13)
    members[[nrow(x) + i]] <- union
  }
  expect_true(all(h$cost >= 0))
})

test_that("insufficient overlap leaves spatial neighbours unresolved", {
  x <- rbind(1:10, 11:20)
  q <- matrix(FALSE, 2, 10)
  q[1, 1:7] <- TRUE
  q[2, 4:10] <- TRUE
  stored <- x
  stored[!q] <- 0
  attr(stored, "coordinate") <- 1:10
  attr(stored, "input_error_gamma") <- 0
  h <- .observed_entry_hierarchy(
    stored, q, matrix(TRUE, 1, 2), 1:2,
    relaxed_full_flux(min_shared = 0.7), 1
  )
  expect_equal(nrow(h$children), 0)
  expect_equal(h$actual_k, 2)
  expect_gt(h$qc$unsupported_comparisons, 0)
  expect_false(h$admission$spatial_adjacency_supplies_evidence)
  expect_false(h$admission$cost_penalty)
})

test_that("the IX numerical tie rule is deterministic under edge order", {
  x <- matrix(1, 4, 8)
  q <- matrix(TRUE, 4, 8)
  attr(x, "coordinate") <- 1:8
  attr(x, "input_error_gamma") <- .capivara_gamma(1)
  mask <- matrix(TRUE, 2, 2)
  rep <- relaxed_full_flux(min_shared = 1, min_contributor = 1)
  edges <- .four_neighbour_edges(mask)
  fits <- lapply(list(seq_len(nrow(edges)), rev(seq_len(nrow(edges))), c(3, 1, 4, 2)),
                 function(order) .observed_entry_hierarchy(x, q, mask, 1:4,
                                                            rep, 1, order))
  expect_identical(fits[[1]]$children, fits[[2]]$children)
  expect_identical(fits[[1]]$children, fits[[3]]$children)
  expect_identical(fits[[1]]$cost, fits[[3]]$cost)
  expect_gt(fits[[1]]$qc$exact_demonstrated_tie_groups, 0)
  expect_match(fits[[1]]$qc$tie_rule, "binary64 forward-error")
})

test_that("non-exact numerical ties use the IX key while separated costs keep order", {
  q <- matrix(TRUE, 3, 8)
  rep <- relaxed_full_flux(min_shared = 1, min_contributor = 1)
  evaluate <- function(last_value, order = NULL) {
    x <- rbind(rep(0, 8), rep(1, 8), rep(last_value, 8))
    attr(x, "coordinate") <- 1:8
    attr(x, "input_error_gamma") <- .capivara_gamma(185)
    .observed_entry_hierarchy(x, q, matrix(TRUE, 1, 3), 1:3,
                              rep, 1, order)
  }
  tied <- evaluate(2 - 2e-15)
  tied_reverse <- evaluate(2 - 2e-15, 2:1)
  expect_identical(tied$children, tied_reverse$children)
  expect_identical(tied$children[1, ], c(1L, 2L))
  expect_identical(tied$qc$numerical_tie_groups, 1L)
  expect_identical(tied$qc$equivalent_selected_not_strict_minimum, 1L)
  expect_false(as.logical(tied$tie_groups$exact_demonstrated[1]))

  separated <- evaluate(1.9)
  expect_identical(separated$children[1, ], c(2L, 3L))
  expect_identical(separated$qc$decisions_with_ties, 0)
})

test_that("raw hierarchy inversions remain visible and unmodified", {
  x <- structure(c(0.23, 0.11, 1.15, 0.9, NA, -0.15, -0.56, NA, 0, 0,
                   1.64, 0.86, 0.17, NA, NA, 0.48, -1.88, 0.19, 0.09, NA,
                   0.56, 0.62, NA, 0.59, 1.63), dim = c(5L, 5L))
  q <- is.finite(x)
  x[!q] <- 0
  attr(x, "coordinate") <- 1:5
  attr(x, "input_error_gamma") <- 0
  h <- .observed_entry_hierarchy(
    x, q, matrix(TRUE, 1, 5), 1:5,
    relaxed_full_flux(min_shared = 0.25, min_contributor = 0.25), 1
  )
  expect_equal(h$cost, c(0.03625, 0.5408, 3.18205, 2.61065), tolerance = 1e-12)
  expect_identical(h$qc$raw_parent_child_inversions, 1L)
  expect_false(h$qc$height_monotonicization)
  expect_false(h$qc$costs_modified)
})

test_that("semantic segmentation preserves disconnected support and measured flux", {
  x <- shape_fixture(nr = 1, nc = 3, gains = c(1, 2, 4))
  support_map <- matrix(c(TRUE, FALSE, TRUE), 1, 3)
  support <- build_capivara_support(
    support_map, support_map, support_map,
    construction_method = "disconnected flux fixture", source = "package regression"
  )
  out <- segment(x$input, Ncomp = 2, var_cube = x$variance,
                 representation = x$profile, sample_validity = x$validity,
                 support = support)
  expect_true(is.na(out$cluster_map[1, 2]))
  expect_equal(out$Ncomp, 2)
  expect_true(out$support$qc$disconnected_support_legal)
  regional <- summarize_cluster_spectra(out)
  flux_matrix <- matrix(x$input$imDat, nrow = 3)
  for (cluster in sort(unique(out$cluster_map[is.finite(out$cluster_map)]))) {
    expected <- colSums(flux_matrix[which(as.vector(out$cluster_map) == cluster), , drop = FALSE])
    got <- unname(regional$sum_spectra[as.character(cluster), ])
    expect_equal(got, expected)
  }
})

test_that("support and analysis provenance hashes are stable and sensitive", {
  x <- shape_fixture()
  support <- complete_support(2, 2)
  p <- prepare_capivara_representation(x$input, x$profile, x$variance,
                                       x$validity, support)
  a <- .capivara_analysis_support(support, p)
  copy <- unserialize(serialize(a, NULL, version = 3))
  expect_identical(copy$analysis_support_id, a$analysis_support_id)
  p$eligible_map[1] <- FALSE
  b <- .capivara_analysis_support(support, p)
  expect_false(identical(a$analysis_support_id, b$analysis_support_id))
  expect_identical(a$sample_validity, p$sample_validity)
  expect_match(a$sample_validity_provenance$selected_validity_id,
               "^validity-sha256-")
  expect_identical(a$sample_validity_provenance$source,
                   "explicit sample_validity argument")
})

test_that("10224-like coverage is restricted while 8602-like coverage is whole", {
  blank <- matrix(0, 21, 21)
  radius <- sqrt((row(blank) - 11)^2 + (col(blank) - 11)^2)
  requested <- radius <= 9
  continuum <- exp(-radius / 2) + 0.01
  support <- build_capivara_support(
    requested, requested, requested, centre = c(11, 11),
    construction_method = "radial coverage fixture", source = "package regression"
  )
  profile <- manga_sandra_spectral_shape()
  make_prepared <- function(eligible) structure(list(
    eligible_map = eligible, measured_continuum = continuum,
    representation = profile
  ), class = c("capivara_prepared_representation", "list"))
  restricted <- requested & (radius <= 7 | (radius > 7 & row(blank) > 16))
  whole <- requested & !(radius > 8.5 & col(blank) < 9)
  q10224 <- representation_coverage_qc(make_prepared(restricted), support,
                                       radius = radius)
  q8602 <- representation_coverage_qc(make_prepared(whole), support,
                                      radius = radius)
  expect_identical(q10224$applicability, "RESTRICTED_DOMAIN")
  expect_identical(q10224$radius_units, "supplied coordinate units")
  expect_lt(q10224$fraction_requested_support_eligible, 0.75)
  expect_gte(q10224$fraction_positive_measured_continuum_retained, 0.90)
  expect_identical(q8602$applicability, "WHOLE_GALAXY")
  expect_gte(q8602$outer_quartile_fraction, 0.5)
})

test_that("the first-order diagnostic is QC and never a Ward cost", {
  x <- shape_fixture()
  support <- complete_support(2, 2)
  p <- prepare_capivara_representation(x$input, x$profile, x$variance,
                                       x$validity, support)
  diagnostic <- spectral_shape_diagnostic(p, 1, 2)
  expect_true(all(is.finite(unlist(diagnostic[c("p_value", "observed", "mean",
                                                "variance", "df", "scale")]))))
  expect_match(diagnostic$role, "excluded from Ward cost")
  anchor <- p$amplitude_channel_mask
  pooled <- 0.5 * (p$features[1, ] + p$features[2, ])
  endpoint_covariance <- function(amplitude, variance) {
    w <- as.numeric(anchor) / sum(anchor)
    derivative <- (diag(length(pooled)) - outer(pooled, w)) / amplitude
    derivative %*% diag(variance) %*% t(derivative)
  }
  covariance <- endpoint_covariance(p$amplitude[1], p$selected_variance[1, ]) +
    endpoint_covariance(p$amplitude[2], p$selected_variance[2, ])
  dense_mean <- sum(diag(covariance)) / nrow(covariance)
  dense_variance <- 2 * sum(covariance^2) / nrow(covariance)^2
  expect_equal(diagnostic$mean, dense_mean, tolerance = 2e-13)
  expect_equal(diagnostic$variance, dense_variance, tolerance = 2e-13)
  out <- segment(x$input, Ncomp = 2, var_cube = x$variance,
                 representation = x$profile, sample_validity = x$validity,
                 support = support)
  expect_false(out$backend_info$diagnostic_in_merge_cost)

  without_frame <- x$input
  without_frame$wavelength_frame <- NULL
  framed <- segment(without_frame, Ncomp = 2, var_cube = x$variance,
                    wavelength_frame = "rest", representation = x$profile,
                    sample_validity = x$validity, support = support)
  expect_identical(framed$cluster_map, out$cluster_map)
})
