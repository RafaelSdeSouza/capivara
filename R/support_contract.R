#' Build wavelength-dependent voxel validity
#'
#' This constructor preserves the native row-column-channel validity array and
#' each exclusion cause. Its spatial quality mask is only a requested-support
#' diagnostic; invalid samples remain missing and do not automatically reject a
#' spaxel from every representation.
#'
#' @param flux Numeric row-column-channel array.
#' @param wavelength Native wavelength coordinate.
#' @param min_valid_fraction Required usable fraction in the selected interval.
#' @param variance Optional variance array matching `flux`.
#' @param ivar Optional inverse-variance array matching `flux`.
#' @param dq Optional integer data-quality array matching `flux`.
#' @param donotuse_bit Integer bit value marking unusable samples.
#' @param wavelength_interval Closed interval evaluated in
#'   `interval_wavelength_frame`.
#' @param wavelength_frame Frame of `wavelength`, `"observed"` or `"rest"`.
#' @param interval_wavelength_frame Frame of `wavelength_interval`.
#' @param redshift Systemic redshift when the two frames differ.
#' @param max_donotuse_fraction Optional additional spatial diagnostic limit.
#' @param provenance Named native-input provenance.
#' @return A `capivara_quality_support` with exact sample validity and cause
#'   masks in canonical row-column-channel order.
#' @export
build_quality_support <- function(flux, wavelength, min_valid_fraction,
                                  variance = NULL, ivar = NULL, dq = NULL,
                                  donotuse_bit = 1024L,
                                  wavelength_interval = range(wavelength),
                                  wavelength_frame = c("observed", "rest"),
                                  interval_wavelength_frame = wavelength_frame,
                                  redshift = NA_real_,
                                  max_donotuse_fraction = NULL,
                                  provenance = list()) {
  wavelength_frame <- match.arg(wavelength_frame)
  interval_wavelength_frame <- match.arg(interval_wavelength_frame, c("observed", "rest"))
  if (!is.array(flux) || length(dim(flux)) != 3L) stop("`flux` must be a row x column x channel array.", call. = FALSE)
  if (length(wavelength) != dim(flux)[3L] || any(!is.finite(wavelength)) || any(diff(wavelength) <= 0)) {
    stop("`wavelength` must be one increasing finite value per channel.", call. = FALSE)
  }
  if (!is.numeric(min_valid_fraction) || length(min_valid_fraction) != 1L ||
      !is.finite(min_valid_fraction) || min_valid_fraction < 0 || min_valid_fraction > 1) {
    stop("`min_valid_fraction` must be one finite number in [0, 1].", call. = FALSE)
  }
  if (!is.null(variance) && !is.null(ivar)) stop("Supply at most one of `variance` and `ivar`.", call. = FALSE)
  for (x in list(variance, ivar, dq)) if (!is.null(x) && !identical(dim(x), dim(flux))) {
    stop("Quality arrays must have the same dimensions as `flux`.", call. = FALSE)
  }
  if (length(wavelength_interval) != 2L || any(!is.finite(wavelength_interval)) ||
      wavelength_interval[2L] < wavelength_interval[1L]) {
    stop("`wavelength_interval` must be an increasing finite pair.", call. = FALSE)
  }
  coordinate <- as.numeric(wavelength)
  if (wavelength_frame != interval_wavelength_frame) {
    if (!.valid_systemic_redshift(redshift)) stop("A finite redshift > -1 is required to change wavelength frame.", call. = FALSE)
    coordinate <- if (wavelength_frame == "observed") coordinate / (1 + redshift) else coordinate * (1 + redshift)
  }
  selected <- is.finite(coordinate) & coordinate >= wavelength_interval[1L] & coordinate <= wavelength_interval[2L]
  if (!any(selected)) stop("No native channel lies in `wavelength_interval`.", call. = FALSE)

  finite_flux <- is.finite(flux)
  valid_variance <- array(TRUE, dim(flux))
  if (!is.null(variance)) valid_variance <- is.finite(variance) & variance > 0
  if (!is.null(ivar)) valid_variance <- is.finite(ivar) & ivar > 0
  donotuse <- array(FALSE, dim(flux))
  if (!is.null(dq)) {
    donotuse <- !is.finite(dq) | bitwAnd(as.integer(dq), as.integer(donotuse_bit)) != 0L
    dim(donotuse) <- dim(flux)
  }
  requested <- array(rep(selected, each = dim(flux)[1L] * dim(flux)[2L]), dim(flux))
  sample_validity <- requested & finite_flux & valid_variance & !donotuse
  selected_index <- which(selected)
  fraction <- function(x) apply(x[, , selected_index, drop = FALSE], c(1L, 2L), mean)
  finite_fraction <- fraction(finite_flux)
  variance_fraction <- fraction(valid_variance)
  donotuse_fraction <- fraction(donotuse)
  valid_fraction <- fraction(sample_validity)
  quality_mask <- valid_fraction >= min_valid_fraction
  if (!is.null(max_donotuse_fraction)) {
    if (length(max_donotuse_fraction) != 1L || !is.finite(max_donotuse_fraction) ||
        max_donotuse_fraction < 0 || max_donotuse_fraction > 1) {
      stop("`max_donotuse_fraction` must be NULL or one number in [0, 1].", call. = FALSE)
    }
    quality_mask <- quality_mask & donotuse_fraction <= max_donotuse_fraction
  }
  out <- list(
    sample_validity = sample_validity,
    exclusions_by_cause = list(
      outside_request = !requested,
      nonfinite_flux = requested & !finite_flux,
      invalid_variance = requested & finite_flux & !valid_variance,
      donotuse = requested & finite_flux & valid_variance & donotuse
    ),
    requested_channel_mask = selected, quality_mask = quality_mask,
    valid_fraction = valid_fraction, finite_fraction = finite_fraction,
    variance_fraction = variance_fraction, donotuse_fraction = donotuse_fraction,
    channel_index = selected_index,
    native_wavelength_limits = range(wavelength[selected]),
    evaluated_wavelength_limits = range(coordinate[selected]),
    wavelength_interval = as.numeric(wavelength_interval),
    wavelength_frame = wavelength_frame,
    interval_wavelength_frame = interval_wavelength_frame,
    redshift = redshift, min_valid_fraction = min_valid_fraction,
    max_donotuse_fraction = max_donotuse_fraction,
    donotuse_bit = as.integer(donotuse_bit),
    axis_order = "row, column, channel",
    conditions = list(finite_flux = TRUE,
                      positive_finite_variance = !is.null(variance) || !is.null(ivar),
                      donotuse_excluded = !is.null(dq)),
    provenance = provenance
  )
  class(out) <- c("capivara_quality_support", "list")
  out$quality_id <- .capivara_stable_id("quality", out)
  out
}

.capivara_components <- function(mask) {
  mask <- as.matrix(mask); mask[is.na(mask)] <- FALSE
  nr <- nrow(mask); nc <- ncol(mask); labels <- matrix(0L, nr, nc); component <- 0L
  for (start in which(mask & labels == 0L)) {
    if (labels[start] != 0L) next
    component <- component + 1L; queue <- start; labels[start] <- component; head <- 1L
    while (head <= length(queue)) {
      k <- queue[head]; head <- head + 1L; ij <- arrayInd(k, .dim = c(nr, nc))
      candidate <- rbind(c(ij[1L] - 1L, ij[2L]), c(ij[1L] + 1L, ij[2L]),
                         c(ij[1L], ij[2L] - 1L), c(ij[1L], ij[2L] + 1L))
      inside <- candidate[, 1L] >= 1L & candidate[, 1L] <= nr &
        candidate[, 2L] >= 1L & candidate[, 2L] <= nc
      candidate <- candidate[inside, , drop = FALSE]
      neighbours <- candidate[, 1L] + (candidate[, 2L] - 1L) * nr
      neighbours <- neighbours[mask[neighbours] & labels[neighbours] == 0L]
      if (length(neighbours)) { labels[neighbours] <- component; queue <- c(queue, neighbours) }
    }
  }
  labels[!mask] <- NA_integer_; labels
}

.capivara_support_id <- function(x) {
  .capivara_stable_id("support", x[setdiff(names(x), c("support_id", "analysis_support_id"))])
}

#' Construct an explicit CAPIVARA support contract
#'
#' Detection, host association, ambiguity, representation eligibility and the
#' final analysis domain remain separate. No connectivity filter, closing,
#' dilation or hole filling is applied; disconnected support is legal.
#'
#' @param quality A `capivara_quality_support` or logical spatial quality map.
#'   For the wavelength-dependent object, its spatial fraction screen remains a
#'   diagnostic and the representation applies its own eligibility rule. A
#'   logical map is treated as an explicit spatial constraint.
#' @param detection_mask Logical detected-signal map.
#' @param host_mask Optional logical host-association map.
#' @param analysis_mask Optional requested final spatial domain. It must lie in
#'   quality and detection support.
#' @param ambiguous_mask Optional unresolved-association map.
#' @param representation_eligibility Optional representation-specific map.
#' @param component_states Optional named states for detected components.
#' @param evidence Optional detection-statistic or surface-brightness map.
#' @param centre Optional adopted row/column centre.
#' @param background_estimate,noise_estimate Recorded detector estimates.
#' @param construction_method Non-empty construction description.
#' @param source Non-empty source identifier.
#' @param configuration,provenance Named lists retained verbatim.
#' @return A serializable `capivara_support` object with stable identities.
#' @export
build_capivara_support <- function(quality, detection_mask, host_mask = NULL,
                                   analysis_mask = NULL, ambiguous_mask = NULL,
                                   representation_eligibility = NULL,
                                   component_states = NULL, evidence = NULL,
                                   centre = NULL,
                                   background_estimate = NA_real_,
                                   noise_estimate = NA_real_,
                                   construction_method, source = "external",
                                   configuration = list(), provenance = list()) {
  if (!is.character(source) || length(source) != 1L || is.na(source) || !nzchar(source)) stop("`source` must be non-empty.", call. = FALSE)
  if (missing(construction_method) || !is.character(construction_method) ||
      length(construction_method) != 1L || is.na(construction_method) || !nzchar(construction_method)) {
    stop("`construction_method` must be non-empty.", call. = FALSE)
  }
  quality_object <- if (inherits(quality, "capivara_quality_support")) quality else NULL
  quality_diagnostic_mask <- if (is.null(quality_object)) quality else quality_object$quality_mask
  # A wavelength-dependent quality object supplies voxel validity. Its spatial
  # fraction screen is diagnostic; the selected representation decides whether
  # the remaining measurements suffice. A directly supplied logical map is an
  # explicit spatial constraint and retains its historical meaning.
  quality_mask <- if (is.null(quality_object)) quality else
    matrix(TRUE, nrow(quality_diagnostic_mask), ncol(quality_diagnostic_mask))
  dims <- dim(quality_mask)
  if (length(dims) != 2L || !is.logical(quality_mask)) stop("`quality` must contain a logical spatial quality mask.", call. = FALSE)
  masks <- list(detection_mask = detection_mask, host_mask = host_mask,
                analysis_mask = analysis_mask, ambiguous_mask = ambiguous_mask,
                representation_eligibility = representation_eligibility)
  for (name in names(masks)) if (!is.null(masks[[name]]) &&
      (!is.logical(masks[[name]]) || !identical(dim(masks[[name]]), dims))) {
    stop("Every supplied support mask must be logical and match `quality`.", call. = FALSE)
  }
  clean <- function(x, default = FALSE) {
    if (is.null(x)) x <- matrix(default, dims[1L], dims[2L])
    x[is.na(x)] <- FALSE; x
  }
  quality_mask <- clean(quality_mask)
  quality_diagnostic_mask <- clean(quality_diagnostic_mask)
  detection_mask <- clean(detection_mask)
  host_supplied <- !is.null(host_mask); host_mask <- clean(host_mask)
  if (any(host_mask & !detection_mask)) stop("`host_mask` must be a subset of `detection_mask`.", call. = FALSE)
  detected_quality <- detection_mask & quality_mask
  if (is.null(ambiguous_mask)) ambiguous_mask <- detected_quality & !host_mask
  ambiguous_mask <- clean(ambiguous_mask)
  if (any(ambiguous_mask & (!detection_mask | host_mask))) stop("Ambiguous pixels must be detected and cannot also be host pixels.", call. = FALSE)
  if (is.null(analysis_mask)) analysis_mask <- if (host_supplied) detected_quality & host_mask else matrix(FALSE, dims[1L], dims[2L])
  analysis_mask <- clean(analysis_mask)
  if (any(analysis_mask & (!quality_mask | !detection_mask))) stop("`analysis_mask` must be a subset of quality and detection support.", call. = FALSE)
  if (!is.null(representation_eligibility)) representation_eligibility <- clean(representation_eligibility)
  component_map <- .capivara_components(detection_mask)
  ids <- sort(unique(stats::na.omit(as.vector(component_map))))
  if (is.null(centre)) {
    idx <- which(analysis_mask, arr.ind = TRUE)
    if (!nrow(idx)) idx <- which(detected_quality, arr.ind = TRUE)
    centre <- if (nrow(idx)) colMeans(idx) else c(NA_real_, NA_real_)
  }
  if (!is.null(evidence) && !identical(dim(evidence), dims)) stop("`evidence` must match support.", call. = FALSE)
  component_table <- do.call(rbind, lapply(ids, function(id) {
    m <- !is.na(component_map) & component_map == id; idx <- which(m, arr.ind = TRUE)
    inferred <- if (any(m & host_mask)) "HOST_ASSOCIATED" else if (any(m & ambiguous_mask)) "AMBIGUOUS" else "NONHOST"
    state <- if (!is.null(component_states) && as.character(id) %in% names(component_states)) unname(component_states[as.character(id)]) else inferred
    data.frame(component_id = id, area = nrow(idx), row_centroid = mean(idx[, 1L]),
      column_centroid = mean(idx[, 2L]), distance_to_centre = sqrt(sum((colMeans(idx) - centre)^2)),
      touches_image_edge = any(idx[, 1L] %in% c(1L, dims[1L]) | idx[, 2L] %in% c(1L, dims[2L])),
      host_pixels = sum(m & host_mask), analysis_pixels = sum(m & analysis_mask),
      ambiguous_pixels = sum(m & ambiguous_mask), state = state,
      median_evidence = if (is.null(evidence)) NA_real_ else stats::median(evidence[m], na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  if (is.null(component_table)) component_table <- data.frame(component_id = integer(), area = integer())
  out <- list(
    quality_mask = quality_mask, quality_diagnostic_mask = quality_diagnostic_mask,
    detection_mask = detection_mask,
    host_mask = host_mask, ambiguous_mask = ambiguous_mask,
    representation_eligibility = representation_eligibility,
    analysis_mask = analysis_mask, connected_component_map = component_map,
    component_table = component_table, centre = as.numeric(centre),
    exclusions_by_cause = list(
      outside_detection = !detection_mask,
      failed_quality = detection_mask & !quality_mask,
      nonhost = detected_quality & !host_mask & !ambiguous_mask,
      ambiguous_association = ambiguous_mask,
      supplied_representation_ineligible = if (is.null(representation_eligibility)) matrix(FALSE, dims[1L], dims[2L]) else analysis_mask & !representation_eligibility
    ),
    background_estimate = background_estimate, noise_estimate = noise_estimate,
    quality = quality_object, construction_method = construction_method,
    source = source, configuration = configuration,
    qc = list(nonempty_analysis = any(analysis_mask),
              analysis_subset_quality = !any(analysis_mask & !quality_mask),
              analysis_subset_detection = !any(analysis_mask & !detection_mask),
              disconnected_support_legal = TRUE,
              component_connectivity = 4L, morphology_operations = "none"),
    provenance = provenance
  )
  class(out) <- c("capivara_support", "list")
  out$support_id <- .capivara_support_id(out)
  out$analysis_support_id <- .capivara_stable_id("analysis", list(
    support_id = out$support_id, analysis_mask = out$analysis_mask,
    representation_eligibility = out$representation_eligibility
  ))
  out
}

.capivara_analysis_support <- function(support, prepared) {
  declared <- if (is.null(support$representation_eligibility))
    matrix(TRUE, nrow(support$analysis_mask), ncol(support$analysis_mask)) else
    support$representation_eligibility
  final <- support$analysis_mask & declared & prepared$eligible_map
  out <- list(
    parent_support_id = support$support_id,
    representation_id = prepared$representation_id,
    requested_analysis_mask = support$analysis_mask,
    supplied_representation_eligibility = support$representation_eligibility,
    representation_eligibility = prepared$eligible_map,
    final_analysis_mask = final,
    exclusions_by_cause = c(support$exclusions_by_cause,
      list(computed_representation_ineligible = support$analysis_mask & !prepared$eligible_map)),
    selected_channel_indices = prepared$selected_channel_indices,
    sample_validity = prepared$sample_validity,
    sample_validity_provenance = prepared$sample_validity_provenance,
    construction = "requested support intersected with frozen representation eligibility"
  )
  out$analysis_support_id <- .capivara_stable_id("analysis", out)
  out
}

.apply_capivara_support <- function(input, support) {
  if (is.null(support)) return(.as_cubedat(input))
  if (!inherits(support, "capivara_support")) stop("`support` must be a `capivara_support` object.", call. = FALSE)
  if (!isTRUE(support$qc$nonempty_analysis)) stop("`support$analysis_mask` is empty; resolve host/analysis support before segmentation.", call. = FALSE)
  cubedat <- .as_cubedat(input)
  if (!identical(dim(support$analysis_mask), dim(cubedat$imDat)[1:2])) stop("Support dimensions do not match the input cube.", call. = FALSE)
  cubedat$imDat[!support$analysis_mask] <- NA_real_
  cubedat
}
