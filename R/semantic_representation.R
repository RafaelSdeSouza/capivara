.capivara_stable_id <- function(prefix, payload) {
  paste0(prefix, "-sha256-", digest::digest(
    payload, algo = "sha256", serialize = TRUE, serializeVersion = 3
  ))
}

.capivara_declaration <- function(x) {
  x$evaluator <- NULL
  if (!is.null(x$amplitude_functional)) x$amplitude_functional$evaluator <- NULL
  x[setdiff(names(x), c("amplitude_id", "representation_id"))]
}

#' Declare a positive-homogeneous amplitude functional
#'
#' The functional supplies the amplitude in `X = F / A[F]`. Its declaration is
#' independent of any one instrument. A user evaluator receives one flux
#' vector, variance vector, validity vector and coordinate vector, and returns
#' `value`, `variance` and `complete`.
#'
#' @param name Stable name for the functional.
#' @param evaluator Function implementing the functional.
#' @param coordinate_domain Named description of the coordinate domain.
#' @param required_support Named list describing required samples or intervals.
#' @param units Units of the returned amplitude.
#' @param validation_domain Named list describing what has been validated.
#' @param provenance Named construction provenance.
#' @return A `capivara_amplitude_functional` declaration.
#' @export
capivara_amplitude_functional <- function(name, evaluator, coordinate_domain,
                                          required_support = list(),
                                          units = "same as input flux",
                                          validation_domain = list(status = "UNVALIDATED"),
                                          provenance = list()) {
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be one non-empty character value.", call. = FALSE)
  }
  if (!is.function(evaluator)) stop("`evaluator` must be a function.", call. = FALSE)
  if (!is.list(coordinate_domain) || is.null(coordinate_domain$kind)) {
    stop("`coordinate_domain` must declare at least `kind`.", call. = FALSE)
  }
  out <- list(
    name = name, method = "user_declared", coordinate_domain = coordinate_domain,
    homogeneity = "A[a F] = a A[F] for a > 0",
    required_support = required_support, units = units,
    validation_domain = validation_domain, provenance = provenance,
    evaluator = evaluator
  )
  class(out) <- c("capivara_amplitude_functional", "list")
  out$amplitude_id <- .capivara_stable_id("amplitude", .capivara_declaration(out))
  out
}

#' Declare a CAPIVARA semantic representation
#'
#' @param type Representation type. `spectral_shape` must be requested
#'   explicitly and requires an amplitude functional.
#' @param coordinate_domain Coordinate kind, frame, medium and optional range.
#' @param invariances Mathematical invariances of the representation.
#' @param retains Quantities intentionally retained.
#' @param removes Quantities intentionally removed.
#' @param required_support Channel and feature-support requirements.
#' @param amplitude_functional Optional `capivara_amplitude_functional`.
#' @param eligibility Named eligibility contract.
#' @param units Named input, output and amplitude units.
#' @param validation_domain Scope in which the declaration is validated.
#' @param provenance Named scientific provenance.
#' @param qc Optional representation-specific QC declaration.
#' @return A serializable `capivara_representation` declaration.
#' @export
capivara_representation <- function(type = c("full_flux", "spectral_shape"),
                                    coordinate_domain, invariances, retains,
                                    removes, required_support,
                                    amplitude_functional = NULL,
                                    eligibility, units, validation_domain,
                                    provenance = list(), qc = list()) {
  type <- match.arg(type)
  if (!is.list(coordinate_domain) || is.null(coordinate_domain$kind)) {
    stop("`coordinate_domain` must declare at least `kind`.", call. = FALSE)
  }
  if (identical(coordinate_domain$kind, "wavelength")) {
    frame <- coordinate_domain$frame
    medium <- coordinate_domain$medium
    if (!is.character(frame) || length(frame) != 1L ||
        !(frame %in% c("observed", "rest")) ||
        !is.character(medium) || length(medium) != 1L ||
        !(medium %in% c("air", "vacuum"))) {
      stop("A wavelength domain must declare an observed/rest frame and air/vacuum medium.", call. = FALSE)
    }
  }
  coordinate_range <- coordinate_domain$range
  if (!is.null(coordinate_range) &&
      (!is.numeric(coordinate_range) || length(coordinate_range) != 2L ||
       any(!is.finite(coordinate_range)) || coordinate_range[2L] < coordinate_range[1L])) {
    stop("A coordinate-domain range must be an increasing finite pair.", call. = FALSE)
  }
  if (type == "spectral_shape" &&
      !inherits(amplitude_functional, "capivara_amplitude_functional")) {
    stop("`spectral_shape` requires a declared amplitude functional.", call. = FALSE)
  }
  fields <- list(invariances = invariances, retains = retains, removes = removes,
                 required_support = required_support, eligibility = eligibility,
                 units = units, validation_domain = validation_domain)
  if (any(!vapply(fields, function(x) length(x) > 0L, logical(1)))) {
    stop("Representation declarations cannot omit semantic contract fields.", call. = FALSE)
  }
  support_fields <- c("min_valid_fraction", "min_reliable_shared_fraction",
                      "min_cluster_contributor_fraction", "feature_windows",
                      "min_feature_fraction")
  if (!all(support_fields %in% names(required_support))) {
    stop("`required_support` must declare validity, overlap, contributor and feature guards.", call. = FALSE)
  }
  fractions <- unlist(required_support[support_fields[support_fields != "feature_windows"]],
                      use.names = FALSE)
  if (!is.numeric(fractions) || any(lengths(required_support[
      support_fields[support_fields != "feature_windows"]]) != 1L) ||
      any(!is.finite(fractions) | fractions < 0 | fractions > 1)) {
    stop("Representation support fractions must be finite scalar values in [0, 1].", call. = FALSE)
  }
  if (!is.list(required_support$feature_windows) ||
      any(!vapply(required_support$feature_windows, function(w) {
        is.numeric(w) && length(w) == 2L && all(is.finite(w)) && w[2L] >= w[1L]
      }, logical(1)))) {
    stop("`feature_windows` must be a list of increasing finite coordinate pairs.", call. = FALSE)
  }
  if (type == "spectral_shape") {
    threshold <- eligibility$measured_amplitude_snr_min
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold < 0) {
      stop("A spectral-shape eligibility contract must declare a non-negative measured amplitude S/N threshold.", call. = FALSE)
    }
  }
  out <- list(
    name = type, type = type, coordinate_domain = coordinate_domain,
    invariances = invariances, retains = retains, removes = removes,
    required_support = required_support,
    amplitude_functional = amplitude_functional,
    eligibility = eligibility, units = units,
    validation_domain = validation_domain, provenance = provenance, qc = qc
  )
  class(out) <- c("capivara_representation", "list")
  out$representation_id <- .capivara_stable_id("representation", .capivara_declaration(out))
  out
}

#' Declare the full-flux semantic representation
#'
#' @param coordinate_domain Coordinate declaration; defaults to native channels.
#' @param required_support Support and comparability requirements.
#' @param units Input/output flux units.
#' @param provenance Named provenance.
#' @return A `full_flux` representation declaration.
#' @export
full_flux_representation <- function(
    coordinate_domain = list(kind = "channel", frame = "native", medium = "native", range = NULL),
    required_support = list(min_valid_fraction = 0.7,
                            min_reliable_shared_fraction = 0.7,
                            min_cluster_contributor_fraction = 0.5,
                            feature_windows = list(), min_feature_fraction = 0.5),
    units = list(input = "declared native flux", output = "same as input",
                 amplitude = "retained in representation"),
    provenance = list()) {
  capivara_representation(
    "full_flux", coordinate_domain,
    invariances = "none declared",
    retains = c("absolute amplitude", "continuum shape",
                "absorption structure", "emission structure"),
    removes = "nothing by construction", required_support = required_support,
    eligibility = list(finite_measured_samples = TRUE,
                       missing_samples_remain_missing = TRUE),
    units = units,
    validation_domain = list(status = "DECLARED_NOT_VALIDATED_BY_EXPERIMENTS_VIII_X"),
    provenance = provenance
  )
}

#' Declare a generic spectral-shape representation
#'
#' @param amplitude_functional Positive-homogeneous amplitude declaration.
#' @param coordinate_domain Coordinate declaration.
#' @param required_support Channel and feature-support requirements.
#' @param eligibility Eligibility declaration.
#' @param units Input, output and amplitude units.
#' @param validation_domain Explicit validation scope.
#' @param provenance Named provenance.
#' @param qc Optional QC declaration.
#' @return A `spectral_shape` representation declaration.
#' @export
spectral_shape_representation <- function(amplitude_functional,
                                          coordinate_domain,
                                          required_support,
                                          eligibility,
                                          units = list(input = "declared flux density",
                                                       output = "dimensionless",
                                                       amplitude = "same as input"),
                                          validation_domain = list(status = "UNVALIDATED"),
                                          provenance = list(), qc = list()) {
  capivara_representation(
    "spectral_shape", coordinate_domain,
    invariances = "F(xi) ~ a F(xi), a > 0",
    retains = c("continuum shape", "relative absorption structure",
                "relative emission structure"),
    removes = "one positive multiplicative amplitude, returned separately",
    required_support = required_support,
    amplitude_functional = amplitude_functional,
    eligibility = eligibility, units = units,
    validation_domain = validation_domain, provenance = provenance, qc = qc
  )
}

.equal_native_union_evaluator <- function(windows) {
  force(windows)
  function(flux, variance, validity, coordinate) {
    anchor <- Reduce(`|`, lapply(windows, function(w) coordinate >= w[1] & coordinate <= w[2]))
    if (!any(anchor)) return(list(value = NA_real_, variance = NA_real_, complete = FALSE, samples = 0L))
    complete <- all(validity[anchor]) && all(is.finite(flux[anchor])) &&
      !is.null(variance) && all(is.finite(variance[anchor]) & variance[anchor] > 0)
    if (!complete) return(list(value = NA_real_, variance = NA_real_, complete = FALSE,
                               samples = sum(anchor)))
    n <- sum(anchor)
    list(value = mean(flux[anchor]), variance = sum(variance[anchor]) / n^2,
         complete = TRUE, samples = n)
  }
}

#' Validated MaNGA/Sandra optical spectral-shape profile
#'
#' This named profile preserves the exact Experiment-X equal-native-sample
#' definition. The four windows are configuration, not CAPIVARA ontology.
#'
#' @return A `spectral_shape` declaration validated for measured amplitude
#'   S/N at least 30 under the stated error model.
#' @export
manga_sandra_spectral_shape <- function() {
  windows <- list(c(5050, 5100), c(5400, 5500), c(6000, 6100), c(6800, 6900))
  domain <- list(kind = "wavelength", frame = "rest", medium = "vacuum",
                 range = c(4800, 7400))
  amplitude <- capivara_amplitude_functional(
    name = "manga_sandra_four_window_equal_native_mean",
    evaluator = .equal_native_union_evaluator(windows),
    coordinate_domain = domain,
    required_support = list(complete_union = windows,
                            weighting = "equal native samples over complete union"),
    units = "same flux-density units as input",
    validation_domain = list(
      status = "VALIDATED", instruments = c("MaNGA", "Sandra optical configuration"),
      measured_amplitude_snr_min = 30,
      spectral_error_model = "independent Gaussian channel errors; supplied diagonal variance"
    ),
    provenance = list(experiment8 = "9c23b7d4", experiment9 = "19d0c60a",
                      experiment10 = "5335c878")
  )
  spectral_shape_representation(
    amplitude, domain,
    required_support = list(
      amplitude_windows = windows,
      feature_windows = list(hbeta = c(4840, 4885), mgb = c(5150, 5200),
                             nad = c(5875, 5910), halpha_nii = c(6540, 6600)),
      min_valid_fraction = 0.7, min_reliable_shared_fraction = 0.7,
      min_cluster_contributor_fraction = 0.5, min_feature_fraction = 0.5,
      missing_samples = "remain missing; no imputation"
    ),
    eligibility = list(amplitude_complete = TRUE, amplitude_positive = TRUE,
                       measured_amplitude_snr_min = 30,
                       fixed_threshold_per_object = TRUE),
    validation_domain = list(
      status = "SHAPE_REPRESENTATION_VALIDATED_FOR_SN_GE_30",
      configuration = "MaNGA/Sandra optical", coordinate_sampling = "native",
      measured_amplitude_snr_min = 30,
      error_model = "independent Gaussian spectral errors; no MaNGA spatial covariance",
      other_modalities = "require their own amplitude functional and validation"
    ),
    provenance = list(experiment10_commit = "5335c878",
                      tie_contract_commit = "19d0c60a"),
    qc = list(diagnostic = "first_order_satterthwaite",
              diagnostic_role = "validation/QC only; never a Ward cost",
              continuum_interval = c(5050, 5500),
              whole_galaxy_gates = c(support = 0.75, positive_continuum = 0.90,
                                     outer_quartile = 0.50, r90 = 0.85),
              catastrophe_gates = c(support = 0.50, positive_continuum = 0.75,
                                     outer_quartile = 0.25, r90 = 0.70))
  )
}

.representation_coordinate <- function(cubedat, representation, redshift) {
  n_wave <- dim(cubedat$imDat)[3]
  domain <- representation$coordinate_domain
  if (identical(domain$kind, "channel")) return(seq_len(n_wave))
  if (!identical(domain$kind, "wavelength")) {
    coordinate <- cubedat$coordinate
    if (is.null(coordinate) || length(coordinate) != n_wave || any(!is.finite(coordinate))) {
      stop("Supply one finite input$coordinate value per channel for this domain.", call. = FALSE)
    }
    return(as.numeric(coordinate))
  }
  wave <- .wavelength_axis(cubedat, n_wave)
  input_frame <- .input_wavelength_frame(cubedat, required = TRUE)
  input_medium <- cubedat$wavelength_medium
  if (is.null(input_medium)) stop("A wavelength representation requires input$wavelength_medium.", call. = FALSE)
  input_medium <- match.arg(input_medium, c("air", "vacuum"))
  target_medium <- match.arg(domain$medium, c("air", "vacuum"))
  if (input_medium != target_medium) wave <- convert_wavelength_medium(wave, input_medium, target_medium)
  if (input_frame != domain$frame) {
    if (!.valid_systemic_redshift(redshift)) {
      stop("Changing representation wavelength frame requires an explicit finite redshift > -1.", call. = FALSE)
    }
    wave <- if (input_frame == "observed") wave / (1 + redshift) else wave * (1 + redshift)
  }
  as.numeric(wave)
}

.representation_variance <- function(var_cube, dims) {
  if (is.null(var_cube)) return(NULL)
  variance <- .as_cubedat(var_cube)$imDat
  if (!identical(dim(variance), dims)) stop("`var_cube` must match the input cube.", call. = FALSE)
  variance
}

.representation_validity <- function(flux, variance, sample_validity, support) {
  if (is.null(sample_validity) && !is.null(support) &&
      inherits(support$quality, "capivara_quality_support")) {
    sample_validity <- support$quality$sample_validity
  }
  if (is.null(sample_validity)) sample_validity <- array(TRUE, dim(flux))
  if (!is.logical(sample_validity) || !identical(dim(sample_validity), dim(flux))) {
    stop("`sample_validity` must be a logical array matching the input cube.", call. = FALSE)
  }
  sample_validity[is.na(sample_validity)] <- FALSE
  sample_validity <- sample_validity & is.finite(flux)
  if (!is.null(variance)) sample_validity <- sample_validity & is.finite(variance) & variance > 0
  sample_validity
}

.window_fraction_matrix <- function(validity, coordinate, windows) {
  if (!length(windows)) return(matrix(TRUE, nrow(validity), 0L))
  out <- vapply(windows, function(w) {
    q <- coordinate >= w[1] & coordinate <= w[2]
    if (!any(q)) rep(FALSE, nrow(validity)) else rowMeans(validity[, q, drop = FALSE])
  }, numeric(nrow(validity)))
  if (is.null(dim(out))) out <- matrix(out, ncol = 1L)
  colnames(out) <- names(windows)
  out
}

#' Prepare flux under an explicit semantic representation
#'
#' @param input FITS-like cube or row-column-channel array.
#' @param representation A `capivara_representation` declaration.
#' @param var_cube Optional variance cube. It is required when eligibility uses
#'   measured amplitude S/N.
#' @param sample_validity Optional exact logical voxel-validity array.
#' @param support Optional `capivara_support` supplying validity provenance.
#' @param redshift Systemic redshift when a frame conversion is required.
#' @return A `capivara_prepared_representation` with missing entries retained as
#'   `NA`, separate amplitudes and explicit eligibility causes.
#' @export
prepare_capivara_representation <- function(input, representation, var_cube = NULL,
                                            sample_validity = NULL, support = NULL,
                                            redshift = NA_real_) {
  if (!inherits(representation, "capivara_representation")) {
    stop("`representation` must be a `capivara_representation` object.", call. = FALSE)
  }
  cubedat <- .as_cubedat(input); flux <- cubedat$imDat; dims <- dim(flux)
  if (!is.array(flux) || length(dims) != 3L) stop("Input must be a row x column x channel cube.", call. = FALSE)
  variance <- .representation_variance(var_cube, dims)
  validity_source <- if (!is.null(sample_validity)) "explicit sample_validity argument" else
    if (!is.null(support) && inherits(support$quality, "capivara_quality_support"))
      "capivara_support quality object" else "derived from finite flux and supplied variance"
  validity <- .representation_validity(flux, variance, sample_validity, support)
  coordinate <- .representation_coordinate(cubedat, representation, redshift)
  request <- representation$coordinate_domain$range
  selected <- if (is.null(request)) seq_along(coordinate) else which(coordinate >= request[1] & coordinate <= request[2])
  if (!length(selected)) stop("No native sample lies in the representation coordinate domain.", call. = FALSE)
  f <- matrix(flux[, , selected, drop = FALSE], ncol = length(selected))
  q <- matrix(validity[, , selected, drop = FALSE], ncol = length(selected))
  v <- if (is.null(variance)) NULL else matrix(variance[, , selected, drop = FALSE], ncol = length(selected))
  coord <- coordinate[selected]
  required <- representation$required_support
  valid_fraction <- rowMeans(q)
  feature_fraction <- .window_fraction_matrix(q, coord, required$feature_windows)
  feature_ok <- if (!ncol(feature_fraction)) rep(TRUE, nrow(q)) else
    apply(feature_fraction >= required$min_feature_fraction, 1L, all)
  support_ok <- valid_fraction >= required$min_valid_fraction & feature_ok
  amplitude <- amplitude_variance <- rep(NA_real_, nrow(f)); complete <- positive <- precision <- rep(TRUE, nrow(f))
  amplitude_samples <- 0L
  if (representation$type == "spectral_shape") {
    if (is.null(v)) stop("`spectral_shape` amplitude S/N eligibility requires `var_cube`.", call. = FALSE)
    evaluator <- representation$amplitude_functional$evaluator
    evaluated <- lapply(seq_len(nrow(f)), function(i) evaluator(f[i, ], v[i, ], q[i, ], coord))
    valid_evaluation <- vapply(evaluated, function(z) {
      is.list(z) && all(c("value", "variance", "complete", "samples") %in% names(z)) &&
        is.numeric(z$value) && length(z$value) == 1L &&
        is.numeric(z$variance) && length(z$variance) == 1L &&
        is.logical(z$complete) && length(z$complete) == 1L && !is.na(z$complete) &&
        is.numeric(z$samples) && length(z$samples) == 1L &&
        is.finite(z$samples) && z$samples >= 0 &&
        (!z$complete || (is.finite(z$value) && is.finite(z$variance) && z$variance >= 0))
    }, logical(1))
    if (!all(valid_evaluation)) {
      stop("The amplitude evaluator must return scalar `value`, `variance`, `complete` and `samples` fields.", call. = FALSE)
    }
    amplitude <- vapply(evaluated, `[[`, numeric(1), "value")
    amplitude_variance <- vapply(evaluated, `[[`, numeric(1), "variance")
    complete <- vapply(evaluated, `[[`, logical(1), "complete")
    amplitude_samples <- max(vapply(evaluated, `[[`, numeric(1), "samples"))
    positive <- complete & is.finite(amplitude) & amplitude > 0
    amplitude_snr <- amplitude / sqrt(amplitude_variance)
    threshold <- representation$eligibility$measured_amplitude_snr_min
    precision <- positive & is.finite(amplitude_snr) & amplitude_snr >= threshold
    features <- sweep(f, 1L, amplitude, "/")
  } else {
    amplitude_snr <- rep(NA_real_, nrow(f)); features <- f
  }
  features[!q] <- NA_real_
  eligible <- support_ok & complete & positive & precision
  anchor <- rep(FALSE, length(coord))
  if (representation$type == "spectral_shape" &&
      length(required$amplitude_windows)) {
    anchor <- Reduce(`|`, lapply(required$amplitude_windows,
                                 function(w) coord >= w[1] & coord <= w[2]))
  }
  continuum <- rep(NA_real_, nrow(f))
  if (!is.null(representation$qc$continuum_interval)) {
    w <- representation$qc$continuum_interval
    use <- coord >= w[1] & coord <= w[2]
    if (any(use)) {
      count <- rowSums(q[, use, drop = FALSE])
      continuum[count > 0] <- rowSums(ifelse(q[, use, drop = FALSE], f[, use, drop = FALSE], 0))[count > 0] / count[count > 0]
    }
  }
  map <- function(x) matrix(x, dims[1], dims[2])
  input_gamma <- if (representation$type == "spectral_shape")
    .capivara_gamma(amplitude_samples + 6L) else .capivara_gamma(1L)
  out <- list(
    representation = representation, representation_id = representation$representation_id,
    features = features, sample_validity = q, selected_flux = f,
    selected_variance = v, selected_coordinate = coord,
    selected_channel_indices = selected, amplitude_channel_mask = anchor,
    amplitude = amplitude, amplitude_variance = amplitude_variance,
    amplitude_snr = amplitude_snr, amplitude_samples = amplitude_samples,
    valid_fraction = valid_fraction, feature_fraction = feature_fraction,
    eligible = eligible, eligible_map = map(eligible),
    measured_continuum = map(continuum), input_error_gamma = input_gamma,
    sample_validity_provenance = list(
      source = validity_source,
      parent_quality_id = if (!is.null(support) &&
        inherits(support$quality, "capivara_quality_support")) support$quality$quality_id else NULL,
      finite_flux_intersection = TRUE,
      positive_finite_variance_intersection = !is.null(variance),
      selected_validity_id = .capivara_stable_id("validity", list(
        selected_channel_indices = selected, validity = q
      ))
    ),
    eligibility_causes = list(
      incomplete_amplitude = map(!complete), nonpositive_amplitude = map(complete & !positive),
      insufficient_amplitude_snr = map(positive & !precision),
      insufficient_valid_fraction = map(valid_fraction < required$min_valid_fraction),
      insufficient_feature_support = map(!feature_ok)
    ),
    wavelength_provenance = list(
      coordinate_kind = representation$coordinate_domain$kind,
      representation_frame = representation$coordinate_domain$frame,
      representation_medium = representation$coordinate_domain$medium,
      requested_range = request, selected_channel_indices = selected,
      selected_coordinate_range = range(coord), resampled = FALSE,
      missing_samples_imputed = FALSE
    ),
    source_dimensions = dims
  )
  class(out) <- c("capivara_prepared_representation", "list")
  out
}

.capivara_gamma <- function(k) {
  u <- 2^-53; ku <- k * u
  ifelse(ku < 1, ku / (1 - ku), Inf)
}

.capivara_pair_admissible <- function(qa, qb, coordinate, required) {
  reliable <- qa & qb
  if (mean(reliable) < required$min_reliable_shared_fraction) return(FALSE)
  windows <- required$feature_windows
  if (length(windows)) for (w in windows) {
    in_window <- coordinate >= w[1] & coordinate <= w[2]
    if (!any(in_window) || mean(reliable[in_window]) < required$min_feature_fraction) return(FALSE)
  }
  TRUE
}

.endpoint_satterthwaite_terms <- function(flux, variance, anchor) {
  A <- mean(flux[anchor]); x <- flux / A; m <- sum(anchor)
  w <- as.numeric(anchor) / m; r <- w * variance
  a2 <- A^2; varA <- sum(variance[anchor]) / m^2
  list(diagonal = variance / a2,
       terms = list(list(c = -1/a2, u = r, v = x),
                    list(c = -1/a2, u = x, v = r),
                    list(c = varA/a2, u = x, v = x)))
}

.satterthwaite_moments <- function(fa, fb, va, vb, anchor) {
  a <- .endpoint_satterthwaite_terms(fa, va, anchor)
  b <- .endpoint_satterthwaite_terms(fb, vb, anchor)
  diagonal <- a$diagonal + b$diagonal; terms <- c(a$terms, b$terms)
  trace <- sum(diagonal) + sum(vapply(terms, function(z) z$c * sum(z$u * z$v), numeric(1)))
  frob <- sum(diagonal^2) + 2 * sum(vapply(
    terms, function(z) z$c * sum(diagonal * z$u * z$v), numeric(1)
  ))
  for (z1 in terms) for (z2 in terms) {
    frob <- frob + z1$c * z2$c * sum(z1$u * z2$u) * sum(z1$v * z2$v)
  }
  p <- length(diagonal); mean <- trace / p; variance <- 2 * frob / p^2
  list(mean = mean, variance = variance, df = 2 * mean^2 / variance,
       scale = variance / (2 * mean))
}

#' First-order Satterthwaite spectral-shape diagnostic
#'
#' This QC statistic is available only for the validated MaNGA/Sandra profile.
#' It never enters the Ward objective and assumes independent Gaussian channel
#' errors with the supplied diagonal variances.
#'
#' @param prepared Result of `prepare_capivara_representation()`.
#' @param index_a,index_b One-based flattened spatial indices.
#' @return Observed mean-square difference, moment-matched parameters and tail
#'   probability, with the diagnostic scope recorded.
#' @export
spectral_shape_diagnostic <- function(prepared, index_a, index_b) {
  if (!inherits(prepared, "capivara_prepared_representation") ||
      prepared$representation$type != "spectral_shape" ||
      prepared$representation$qc$diagnostic != "first_order_satterthwaite" ||
      prepared$representation$validation_domain$status !=
        "SHAPE_REPRESENTATION_VALIDATED_FOR_SN_GE_30") {
    stop("`first_order_satterthwaite` is restricted to the validated optical profile.", call. = FALSE)
  }
  ids <- c(index_a, index_b)
  if (any(ids != as.integer(ids)) || any(ids < 1L) || any(ids > nrow(prepared$features))) {
    stop("Diagnostic indices must identify two prepared spatial samples.", call. = FALSE)
  }
  if (!all(prepared$eligible[ids])) stop("Both samples must satisfy representation eligibility.", call. = FALSE)
  q <- prepared$sample_validity[index_a, ] & prepared$sample_validity[index_b, ]
  if (!.capivara_pair_admissible(prepared$sample_validity[index_a, ],
                                 prepared$sample_validity[index_b, ],
                                 prepared$selected_coordinate,
                                 prepared$representation$required_support)) {
    stop("The two samples are not comparable under the overlap contract.", call. = FALSE)
  }
  anchor <- prepared$amplitude_channel_mask[q]
  fa <- prepared$selected_flux[index_a, q]; fb <- prepared$selected_flux[index_b, q]
  va <- prepared$selected_variance[index_a, q]; vb <- prepared$selected_variance[index_b, q]
  Aa <- prepared$amplitude[index_a]; Ab <- prepared$amplitude[index_b]
  xa <- fa / Aa; xb <- fb / Ab; observed <- mean((xa - xb)^2)
  pooled <- 0.5 * (xa + xb)
  moments <- .satterthwaite_moments(Aa * pooled, Ab * pooled, va, vb, anchor)
  c(list(p_value = stats::pchisq(observed / moments$scale, moments$df,
                                 lower.tail = FALSE), observed = observed), moments,
    list(diagnostic = "first_order_satterthwaite",
         role = "validation/QC only; excluded from Ward cost",
         error_model = "independent Gaussian spectral errors; no spatial covariance"))
}

#' Audit representation coverage on a requested support
#'
#' @param prepared Prepared representation.
#' @param support `capivara_support` defining the requested analysis support.
#' @param continuum Optional measured-continuum map; defaults to the prepared
#'   profile diagnostic.
#' @param radius Optional radius map. If absent, pixel radius from the support
#'   centre is used.
#' @param radius_units Units of a supplied radius map. Defaults to `"pixels"`
#'   for an internally constructed radius and `"supplied coordinate units"`
#'   otherwise.
#' @return Coverage fractions, frozen gates and an applicability decision.
#' @export
representation_coverage_qc <- function(prepared, support, continuum = NULL,
                                       radius = NULL, radius_units = NULL) {
  if (!inherits(prepared, "capivara_prepared_representation") ||
      !inherits(support, "capivara_support")) stop("Prepared representation and support are required.", call. = FALSE)
  requested <- support$analysis_mask
  if (!identical(dim(requested), dim(prepared$eligible_map))) stop("Support dimensions differ from the prepared cube.", call. = FALSE)
  if (!any(requested)) stop("Coverage QC requires non-empty requested support.", call. = FALSE)
  eligible <- requested & prepared$eligible_map
  if (is.null(continuum)) continuum <- prepared$measured_continuum
  if (!is.matrix(continuum) || !identical(dim(continuum), dim(requested))) stop("`continuum` must match support.", call. = FALSE)
  radius_supplied <- !is.null(radius)
  if (!radius_supplied) {
    centre <- support$centre
    if (is.null(centre) || length(centre) != 2L || any(!is.finite(centre))) {
      idx <- which(requested, arr.ind = TRUE); centre <- colMeans(idx)
    }
    radius <- sqrt((row(requested) - centre[1])^2 + (col(requested) - centre[2])^2)
  }
  if (!is.matrix(radius) || !identical(dim(radius), dim(requested))) stop("`radius` must match support.", call. = FALSE)
  if (is.null(radius_units)) radius_units <- if (radius_supplied)
    "supplied coordinate units" else "pixels"
  if (!is.character(radius_units) || length(radius_units) != 1L ||
      is.na(radius_units) || !nzchar(radius_units)) {
    stop("`radius_units` must be one non-empty character value.", call. = FALSE)
  }
  positive <- ifelse(requested & is.finite(continuum), pmax(continuum, 0), 0)
  spatial <- sum(eligible) / sum(requested)
  continuum_fraction <- if (sum(positive) > 0) sum(positive[eligible]) / sum(positive) else NA_real_
  threshold <- stats::quantile(radius[requested], 0.75, names = FALSE)
  outer <- requested & radius >= threshold
  outer_fraction <- sum(eligible & outer) / sum(outer)
  r90 <- stats::quantile(radius[requested], 0.9, names = FALSE)
  eligible_r90 <- if (any(eligible)) stats::quantile(radius[eligible], 0.9, names = FALSE) else 0
  r90_ratio <- if (r90 > 0) eligible_r90 / r90 else as.numeric(any(eligible))
  metrics <- c(support = spatial, positive_continuum = continuum_fraction,
               outer_quartile = outer_fraction, r90 = r90_ratio)
  gates <- prepared$representation$qc$whole_galaxy_gates
  catastrophe <- prepared$representation$qc$catastrophe_gates
  whole <- length(gates) == 4L && all(is.finite(metrics)) && all(metrics >= gates)
  catastrophic <- length(catastrophe) == 4L &&
    (any(!is.finite(metrics)) || any(metrics < catastrophe))
  applicability <- if (whole) "WHOLE_GALAXY" else if (catastrophic) "NOT_CERTIFIED" else "RESTRICTED_DOMAIN"
  list(
    requested_support = sum(requested), eligible = sum(eligible),
    fraction_requested_support_eligible = unname(spatial),
    fraction_positive_measured_continuum_retained = unname(continuum_fraction),
    outer_quartile_spaxels = sum(outer), outer_quartile_retained = sum(eligible & outer),
    outer_quartile_fraction = unname(outer_fraction), r90_support = unname(r90),
    r90_eligible = unname(eligible_r90), r90_retention_ratio = unname(r90_ratio),
    radius_units = radius_units,
    whole_galaxy_gates = gates, catastrophe_gates = catastrophe,
    applicability = applicability
  )
}
