.four_neighbour_edges <- function(mask) {
  mask <- as.matrix(mask); mask[is.na(mask)] <- FALSE
  ids <- matrix(0L, nrow(mask), ncol(mask)); ids[mask] <- seq_len(sum(mask))
  edges <- matrix(integer(), 0L, 2L)
  if (nrow(mask) > 1L) {
    q <- mask[-nrow(mask), , drop = FALSE] & mask[-1L, , drop = FALSE]
    if (any(q)) edges <- rbind(edges, cbind(ids[-nrow(mask), , drop = FALSE][q], ids[-1L, , drop = FALSE][q]))
  }
  if (ncol(mask) > 1L) {
    q <- mask[, -ncol(mask), drop = FALSE] & mask[, -1L, drop = FALSE]
    if (any(q)) edges <- rbind(edges, cbind(ids[, -ncol(mask), drop = FALSE][q], ids[, -1L, drop = FALSE][q]))
  }
  matrix(as.integer(edges), ncol = 2L)
}

.representation_feature_groups <- function(coordinate, windows) {
  group <- integer(length(coordinate))
  if (length(windows)) for (i in seq_along(windows)) {
    w <- windows[[i]]; group[coordinate >= w[1L] & coordinate <= w[2L]] <- i
  }
  group
}

.leaf_signatures <- function(features, validity) {
  vapply(seq_len(nrow(features)), function(i) digest::digest(
    list(validity = validity[i, ], values = ifelse(validity[i, ], features[i, ], 0)),
    algo = "sha256", serialize = TRUE, serializeVersion = 3
  ), character(1))
}

.observed_entry_hierarchy <- function(features, validity, final_mask,
                                      original_ids, representation, Ncomp,
                                      edge_order = NULL) {
  if (!is.matrix(features) || !is.matrix(validity) || !identical(dim(features), dim(validity))) {
    stop("Feature and validity matrices must agree.", call. = FALSE)
  }
  edges <- .four_neighbour_edges(final_mask)
  if (!is.null(edge_order)) edges <- edges[edge_order, , drop = FALSE]
  required <- representation$required_support
  fit <- capivara_observed_ward_cpp(
    features, validity, edges, as.integer(original_ids),
    .representation_feature_groups(
      attr(features, "coordinate"), required$feature_windows
    ),
    required$min_reliable_shared_fraction,
    required$min_cluster_contributor_fraction,
    required$min_feature_fraction,
    attr(features, "input_error_gamma"),
    .leaf_signatures(features, validity), as.integer(Ncomp)
  )
  fit$admission <- list(
    min_reliable_shared_fraction = required$min_reliable_shared_fraction,
    min_cluster_contributor_fraction = required$min_cluster_contributor_fraction,
    feature_windows = required$feature_windows,
    min_feature_fraction = required$min_feature_fraction,
    cost_penalty = FALSE, spatial_adjacency_supplies_evidence = FALSE
  )
  fit$objective <- list(
    name = "observed_entry_cluster_sse",
    merge_increment = "sum_j n_Aj*n_Bj/(n_Aj+n_Bj)*(mu_Aj-mu_Bj)^2",
    sufficient_statistics = c("channel count", "channel sum", "absolute-sum error accumulator"),
    imputation = "none", raw_heights = TRUE, monotonicized = FALSE
  )
  class(fit) <- c("capivara_observed_hierarchy", "list")
  fit
}

.segment_semantic <- function(input, Ncomp, redshift, var_cube,
                              representation, sample_validity, support,
                              wavelength_frame = NULL,
                              return_details = FALSE) {
  if (!inherits(support, "capivara_support")) {
    stop("Explicit semantic segmentation requires a `capivara_support` object.", call. = FALSE)
  }
  if (!isTRUE(support$qc$nonempty_analysis)) stop("`support$analysis_mask` is empty.", call. = FALSE)
  if (!is.numeric(Ncomp) || length(Ncomp) != 1L || !is.finite(Ncomp) ||
      Ncomp != as.integer(Ncomp) || Ncomp < 1L) stop("`Ncomp` must be a positive integer.", call. = FALSE)
  raw <- .as_cubedat(input)
  if (!is.null(wavelength_frame)) {
    wavelength_frame <- match.arg(wavelength_frame, c("observed", "rest"))
    if (!is.null(raw$wavelength_frame) &&
        !identical(raw$wavelength_frame, wavelength_frame)) {
      stop("`wavelength_frame` conflicts with input$wavelength_frame; correct the metadata explicitly.", call. = FALSE)
    }
    raw$wavelength_frame <- wavelength_frame
  }
  prepared <- prepare_capivara_representation(
    raw, representation, var_cube = var_cube,
    sample_validity = sample_validity, support = support, redshift = redshift
  )
  analysis_support <- .capivara_analysis_support(support, prepared)
  final <- analysis_support$final_analysis_mask
  original_ids <- which(final)
  if (!length(original_ids)) stop("No representation-eligible samples remain on the requested support.", call. = FALSE)
  if (Ncomp > length(original_ids)) stop("`Ncomp` exceeds the representation-eligible support.", call. = FALSE)
  feature_rows <- prepared$features[original_ids, , drop = FALSE]
  validity_rows <- prepared$sample_validity[original_ids, , drop = FALSE]
  attr(feature_rows, "coordinate") <- prepared$selected_coordinate
  attr(feature_rows, "input_error_gamma") <- prepared$input_error_gamma
  hierarchy <- .observed_entry_hierarchy(
    feature_rows, validity_rows, final, original_ids,
    representation, Ncomp
  )
  cluster_map <- matrix(NA_integer_, dim(raw$imDat)[1L], dim(raw$imDat)[2L])
  cluster_map[original_ids] <- hierarchy$labels
  if (hierarchy$actual_k > Ncomp) warning(
    "The overlap contract leaves ", hierarchy$actual_k,
    " unresolved roots; spatial adjacency cannot supply missing spectral evidence.",
    call. = FALSE
  )
  # Preserve historical cluster-S/N output without using it in the hierarchy.
  sn <- .compute_signal_noise(cube_to_matrix(raw))
  cluster_snr <- .compute_cluster_snr(hierarchy$labels,
                                      sn$signal[original_ids], sn$noise[original_ids])
  declaration <- .capivara_declaration(representation)
  coverage <- if (length(representation$qc$whole_galaxy_gates))
    representation_coverage_qc(prepared, support) else NULL
  out <- list(
    cluster_map = cluster_map, header = raw$hdr, axDat = raw$axDat,
    cluster_snr = cluster_snr, Ncomp = hierarchy$actual_k,
    requested_Ncomp = as.integer(Ncomp), original_cube = raw,
    representation = declaration,
    representation_id = representation$representation_id,
    amplitude = matrix(prepared$amplitude, dim(raw$imDat)[1L], dim(raw$imDat)[2L]),
    amplitude_variance = matrix(prepared$amplitude_variance, dim(raw$imDat)[1L], dim(raw$imDat)[2L]),
    amplitude_snr = matrix(prepared$amplitude_snr, dim(raw$imDat)[1L], dim(raw$imDat)[2L]),
    representation_eligibility = prepared$eligible_map,
    support = support, support_id = support$support_id,
    analysis_support = analysis_support,
    analysis_support_id = analysis_support$analysis_support_id,
    support_provenance = support$provenance,
    coverage_qc = coverage, hierarchy = hierarchy,
    wavelength_provenance = prepared$wavelength_provenance,
    backend = "observed_entry_spatial_ward",
    backend_info = c(hierarchy$qc, list(
      representation_id = representation$representation_id,
      support_id = support$support_id,
      analysis_support_id = analysis_support$analysis_support_id,
      valid_pixels = length(original_ids),
      diagnostic_in_merge_cost = FALSE,
      fixed_k_role = "backward-compatible chronological cut; no physical privilege"
    ))
  )
  if (return_details) out$prepared_representation <- prepared
  out
}
