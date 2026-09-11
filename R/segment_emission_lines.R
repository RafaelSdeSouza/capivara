.capivara_emission_line_table <- function() {
  data.frame(
    name = c(
      "oii3727", "neiii3869", "hdelta", "hgamma", "hbeta",
      "oiii4959", "oiii5007", "oi6300", "halpha",
      "nii6548", "nii6583", "sii6716", "sii6731"
    ),
    label = c(
      "[O II] 3727", "[Ne III] 3869", "Hdelta", "Hgamma", "Hbeta",
      "[O III] 4959", "[O III] 5007", "[O I] 6300", "Halpha",
      "[N II] 6548", "[N II] 6583", "[S II] 6716", "[S II] 6731"
    ),
    rest_wavelength = c(
      3727.09, 3869.86, 4101.74, 4340.47, 4861.33,
      4958.91, 5006.84, 6300.30, 6562.80,
      6548.05, 6583.45, 6716.44, 6730.82
    ),
    family = c(
      "blue", "blue", "balmer", "balmer", "balmer",
      "agn", "agn", "agn", "balmer",
      "agn", "agn", "agn", "agn"
    ),
    stringsAsFactors = FALSE
  )
}

.capivara_match_emission_lines <- function(lines) {
  tab <- .capivara_emission_line_table()
  aliases <- c(
    oii = "oii3727", oii3726 = "oii3727", oii3729 = "oii3727",
    neiii = "neiii3869",
    hd = "hdelta", hdelta4102 = "hdelta",
    hg = "hgamma", hgamma4340 = "hgamma",
    hb = "hbeta", hbeta4861 = "hbeta",
    oiii = "oiii5007", oiii5007 = "oiii5007", oiii4959 = "oiii4959",
    oi = "oi6300", oi6300 = "oi6300",
    ha = "halpha", halpha6563 = "halpha",
    nii = "nii6583", nii6583 = "nii6583", nii6548 = "nii6548",
    sii = "sii6716", sii6716 = "sii6716", sii6731 = "sii6731"
  )

  if (length(lines) == 1L && tolower(lines) %in% c("agn", "default")) {
    lines <- c("hbeta", "oiii4959", "oiii5007", "oi6300", "halpha", "nii6548", "nii6583", "sii6716", "sii6731")
  } else if (length(lines) == 1L && tolower(lines) %in% c("strong", "optical")) {
    lines <- c("oii3727", "hbeta", "oiii5007", "halpha", "nii6583", "sii6716", "sii6731")
  }

  keys <- tolower(gsub("[^a-z0-9]+", "", lines))
  keys <- unname(ifelse(keys %in% names(aliases), aliases[keys], keys))

  unknown <- setdiff(keys, tab$name)
  if (length(unknown)) {
    stop(
      "Unknown emission line(s): ", paste(unknown, collapse = ", "),
      ". Use `emission_lines()` to see available names.",
      call. = FALSE
    )
  }

  tab[match(keys, tab$name), , drop = FALSE]
}

.capivara_line_window_index <- function(wavelengths, center, half_width_kms) {
  c_kms <- 299792.458
  dv <- c_kms * (wavelengths / center - 1)
  which(is.finite(dv) & abs(dv) <= half_width_kms)
}

.capivara_continuum_window_index <- function(wavelengths, center, line_half_width_kms,
                                             continuum_inner_kms, continuum_outer_kms) {
  c_kms <- 299792.458
  dv <- c_kms * (wavelengths / center - 1)
  which(
    is.finite(dv) &
      abs(dv) >= max(line_half_width_kms, continuum_inner_kms) &
      abs(dv) <= continuum_outer_kms
  )
}

.capivara_safe_row_median <- function(x) {
  out <- matrixStats::rowMedians(x, na.rm = TRUE)
  out[!is.finite(out)] <- 0
  out
}

.capivara_scale_columns <- function(x) {
  center <- apply(x, 2, stats::median, na.rm = TRUE)
  scale <- vapply(seq_len(ncol(x)), function(j) {
    stats::mad(x[, j], center = center[j], constant = 1.4826, na.rm = TRUE)
  }, numeric(1))
  fallback <- stats::median(scale[is.finite(scale) & scale > 0], na.rm = TRUE)
  if (!is.finite(fallback) || fallback <= 0) fallback <- 1
  scale[!is.finite(scale) | scale <= 0] <- fallback
  x <- sweep(x, 2, center, "-")
  x <- sweep(x, 2, scale, "/")
  x[!is.finite(x)] <- 0
  x
}

.capivara_emission_feature_cube <- function(input,
                                            redshift,
                                            lines = "agn",
                                            line_window_kms = 500,
                                            continuum_inner_kms = 800,
                                            continuum_outer_kms = 2200,
                                            profile_bins = 9,
                                            feature_mode = c("windows", "window_derivative", "hybrid", "moments", "profile", "profile_derivative"),
                                            continuum_subtract = FALSE,
                                            positive_only = FALSE,
                                            line_weights = NULL) {
  feature_mode <- match.arg(feature_mode)

  if (!is.finite(redshift)) {
    stop("`redshift` must be a finite numeric value for emission-line segmentation.", call. = FALSE)
  }

  cubedat <- .as_cubedat(input)
  cube <- cubedat$imDat
  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("`input$imDat` must be a 3D array with dimensions (n_row, n_col, n_wave).", call. = FALSE)
  }

  dims <- dim(cube)
  wavelengths <- .wavelength_axis(cubedat, dims[3])
  tab <- .capivara_match_emission_lines(lines)
  tab$observed_wavelength <- tab$rest_wavelength * (1 + redshift)
  covered <- tab$observed_wavelength >= min(wavelengths, na.rm = TRUE) &
    tab$observed_wavelength <= max(wavelengths, na.rm = TRUE)
  tab <- tab[covered, , drop = FALSE]
  if (!nrow(tab)) {
    stop("None of the requested redshifted emission lines are covered by the cube wavelength axis.", call. = FALSE)
  }

  if (is.null(line_weights)) {
    line_weights <- rep(1, nrow(tab))
  } else {
    if (is.null(names(line_weights))) {
      if (length(line_weights) != nrow(tab)) {
        stop("Unnamed `line_weights` must have one value per covered line.", call. = FALSE)
      }
    } else {
      lw <- rep(1, nrow(tab))
      names(lw) <- tab$name
      matched <- intersect(names(line_weights), names(lw))
      lw[matched] <- line_weights[matched]
      line_weights <- lw
    }
  }

  mat <- matrix(cube, nrow = dims[1] * dims[2], ncol = dims[3])
  features <- list()
  measured_valid <- rep(FALSE, nrow(mat))
  feature_info <- data.frame()

  for (i in seq_len(nrow(tab))) {
    line_idx <- .capivara_line_window_index(wavelengths, tab$observed_wavelength[i], line_window_kms)
    cont_idx <- .capivara_continuum_window_index(
      wavelengths,
      tab$observed_wavelength[i],
      line_window_kms,
      continuum_inner_kms,
      continuum_outer_kms
    )
    if (length(line_idx) < 2L) next

    line_flux <- mat[, line_idx, drop = FALSE]
    measured_valid <- measured_valid | rowSums(is.finite(line_flux)) >= 2L
    if (isTRUE(continuum_subtract)) {
      continuum <- if (length(cont_idx)) {
        .capivara_safe_row_median(mat[, cont_idx, drop = FALSE])
      } else {
        .capivara_safe_row_median(line_flux)
      }
      resid <- sweep(line_flux, 1, continuum, "-")
    } else {
      resid <- line_flux
    }

    if (isTRUE(positive_only)) {
      resid[!is.finite(resid) | resid < 0] <- 0
    } else {
      resid[!is.finite(resid)] <- 0
    }

    dv <- 299792.458 * (wavelengths[line_idx] / tab$observed_wavelength[i] - 1)
    moments <- NULL
    if (feature_mode %in% c("hybrid", "moments")) {
      total <- rowSums(resid, na.rm = TRUE)
      centroid <- rowSums(sweep(resid, 2, dv, "*"), na.rm = TRUE) / pmax(total, .Machine$double.eps)
      centroid[!is.finite(centroid) | total <= 0] <- 0
      width <- sqrt(abs(rowSums(sweep(resid, 2, dv^2, "*"), na.rm = TRUE) / pmax(total, .Machine$double.eps)))
      width[!is.finite(width) | total <= 0] <- 0
      blue <- rowSums(resid[, dv < 0, drop = FALSE], na.rm = TRUE)
      red <- rowSums(resid[, dv > 0, drop = FALSE], na.rm = TRUE)
      asym <- (red - blue) / pmax(red + blue, .Machine$double.eps)
      asym[!is.finite(asym)] <- 0
      peak <- matrixStats::rowMaxs(resid, na.rm = TRUE)
      peak[!is.finite(peak)] <- 0

      moments <- cbind(
        log_flux = log1p(pmax(total, 0)),
        centroid = centroid,
        width = width,
        asymmetry = asym,
        peak = log1p(pmax(peak, 0))
      )
    }

    profile <- NULL
    if (feature_mode %in% c("windows", "window_derivative", "hybrid", "profile", "profile_derivative")) {
      grid <- seq(min(dv), max(dv), length.out = profile_bins)
      profile <- t(vapply(seq_len(nrow(resid)), function(j) {
        stats::approx(dv, resid[j, ], xout = grid, rule = 2, ties = "ordered")$y
      }, numeric(profile_bins)))
      profile[!is.finite(profile)] <- 0
      if (feature_mode %in% c("profile", "profile_derivative")) {
        profile <- profile / pmax(matrixStats::rowSums2(abs(profile), na.rm = TRUE), .Machine$double.eps)
        profile[!is.finite(profile)] <- 0
      }
      if (feature_mode %in% c("window_derivative", "profile_derivative")) {
        profile <- t(apply(profile, 1, diff))
      }
    }

    block <- switch(
      feature_mode,
      windows = profile,
      window_derivative = profile,
      moments = moments,
      profile = profile,
      profile_derivative = profile,
      hybrid = cbind(moments, profile)
    )
    block <- .capivara_scale_columns(block) * line_weights[i]
    names <- paste(tab$name[i], seq_len(ncol(block)), sep = "_")
    colnames(block) <- names
    features[[length(features) + 1L]] <- block
    feature_info <- rbind(
      feature_info,
      data.frame(
        line = tab$name[i],
        label = tab$label[i],
        rest_wavelength = tab$rest_wavelength[i],
        observed_wavelength = tab$observed_wavelength[i],
        n_line_channels = length(line_idx),
        n_continuum_channels = length(cont_idx),
        weight = line_weights[i],
        stringsAsFactors = FALSE
      )
    )
  }

  if (!length(features)) {
    stop("No requested emission line had enough wavelength channels for feature construction.", call. = FALSE)
  }

  feature_mat <- do.call(cbind, features)
  feature_mat[!is.finite(feature_mat)] <- 0
  feature_cube <- array(feature_mat, dim = c(dims[1], dims[2], ncol(feature_mat)))

  list(
    cube = list(imDat = feature_cube, hdr = cubedat$hdr, axDat = NULL),
    feature_matrix = feature_mat,
    feature_info = feature_info,
    line_table = tab,
    wavelengths = wavelengths,
    measured_valid = measured_valid,
    original_cube = cubedat
  )
}

#' List built-in emission lines for Capivara line-sensitive segmentation
#'
#' @return A data frame with line names, labels, rest wavelengths, and families.
#' @export
emission_lines <- function() {
  .capivara_emission_line_table()
}

#' Segment an IFU cube using emission-line features
#'
#' `segment_emission_lines()` is an off-the-shelf segmentation mode for cubes
#' where small line-emitting structures matter more than the full continuum
#' shape. By default it does not line-fit, continuum-subtract, or build a
#' support mask: it simply clusters the cube using spectral channels around
#' redshifted optical/AGN emission lines, then runs [segment_large()].
#'
#' @param input FITS-like cube object with `imDat` and `axDat`, or a raw cube.
#' @param redshift Numeric redshift used to place the rest-frame lines.
#' @param lines Character vector of line names, or `"agn"` for H beta, [O III],
#'   [O I], H alpha, [N II], and [S II] where covered.
#' @param Ncomp Number of output segments.
#' @param line_window_kms Half-width of each line window in km/s.
#' @param continuum_inner_kms Inner half-width excluded from continuum sidebands.
#' @param continuum_outer_kms Outer half-width for continuum sidebands.
#' @param feature_mode `"windows"` clusters flux samples in the redshifted line
#'   windows and is the default. `"window_derivative"` clusters local spectral
#'   slopes. `"hybrid"` combines line moments and normalized line profile
#'   shape. `"moments"` is fastest. `"profile"` emphasizes line shape.
#'   `"profile_derivative"` emphasizes profile changes/asymmetries.
#' @param profile_bins Number of bins used for normalized profile features.
#' @param continuum_subtract Optional local side-band continuum subtraction.
#'   Defaults to `FALSE` so the off-the-shelf behaviour is only line-window
#'   selection, not preprocessing.
#' @param positive_only If `TRUE`, negative continuum residuals are clipped to
#'   zero before features are computed. Defaults to `FALSE`.
#' @param line_weights Optional named or positional weights for covered lines.
#' @param knn_k Number of neighbours used by sparse Ward.
#' @param spatial_weight Spatial regularization passed to [segment_large()].
#' @param valid_mode Valid-pixel rule passed to [segment_large()]. `"finite"`
#'   is the default because line-window features are column-scaled and may be
#'   negative even for useful spectra. No support mask is supplied.
#' @param max_pixels Maximum number of valid spaxels used to learn the sparse
#'   Ward graph. If the cube has more valid spaxels, Capivara learns regions on
#'   a reproducible sample and assigns every valid spaxel to the nearest learned
#'   region. This keeps no-mask emission-line segmentation usable on large IFUs.
#'   Use `Inf` to force the full graph.
#' @param seed Random seed used only when `max_pixels` triggers sampling.
#' @param ... Additional arguments passed to [segment_large()].
#' @return A `segment_large` result with emission-line feature metadata.
#' @export
segment_emission_lines <- function(input,
                                   redshift,
                                   lines = "agn",
                                   Ncomp = 50,
                                   line_window_kms = 500,
                                   continuum_inner_kms = 800,
                                   continuum_outer_kms = 2200,
                                   feature_mode = c("windows", "window_derivative", "hybrid", "moments", "profile", "profile_derivative"),
                                   profile_bins = 9,
                                   continuum_subtract = FALSE,
                                   positive_only = FALSE,
                                   line_weights = NULL,
                                   knn_k = 30,
                                   spatial_weight = 0.05,
                                   valid_mode = c("finite", "signal", "sagui"),
                                   max_pixels = 30000,
                                   seed = 1L,
                                   ...) {
  feature_mode <- match.arg(feature_mode)
  valid_mode <- match.arg(valid_mode)

  features <- .capivara_emission_feature_cube(
    input = input,
    redshift = redshift,
    lines = lines,
    line_window_kms = line_window_kms,
    continuum_inner_kms = continuum_inner_kms,
    continuum_outer_kms = continuum_outer_kms,
    profile_bins = profile_bins,
    feature_mode = feature_mode,
    continuum_subtract = continuum_subtract,
    positive_only = positive_only,
    line_weights = line_weights
  )

  n_row <- dim(features$cube$imDat)[1]
  n_col <- dim(features$cube$imDat)[2]
  feature_mat <- features$feature_matrix
  finite_counts <- rowSums(is.finite(feature_mat))
  finite_frac <- finite_counts / ncol(feature_mat)
  row_energy <- rowSums(feature_mat^2, na.rm = TRUE)
  if (identical(valid_mode, "finite")) {
    valid <- finite_counts == ncol(feature_mat)
  } else if (identical(valid_mode, "signal")) {
    valid <- rowSums(abs(feature_mat), na.rm = TRUE) > 0
  } else {
    row_mad <- apply(feature_mat, 1, function(v) {
      vv <- v[is.finite(v)]
      if (!length(vv)) return(NA_real_)
      stats::mad(vv, na.rm = TRUE)
    })
    valid <- finite_counts >= pmin(10L, ncol(feature_mat)) &
      finite_frac >= 0.8 &
      is.finite(row_energy) & row_energy > 0 &
      is.finite(row_mad) & row_mad > 0
  }
  valid <- valid & features$measured_valid
  valid_indices <- which(valid)
  if (!length(valid_indices)) {
    stop("No valid pixels after emission-line feature construction.", call. = FALSE)
  }

  # Match the pre-existing unsampled row scaling in both computational paths.
  feature_mat <- .sparse_ward_scale_features(feature_mat, scale_fn = NULL)
  features$cube$imDat <- array(feature_mat, dim(features$cube$imDat))
  sampled <- is.finite(max_pixels) && length(valid_indices) > max_pixels
  if (!sampled) {
    out <- segment_large(
      features$cube,
      Ncomp = Ncomp,
      use_starlet_mask = FALSE,
      mask = matrix(valid, n_row, n_col),
      valid_mode = valid_mode,
      knn_k = knn_k,
      spatial_weight = spatial_weight,
      feature_scale = "none",
      scale_fn = identity,
      ...
    )
  } else {
    set.seed(seed)
    sample_n <- as.integer(max_pixels)
    sample_pos <- sort(sample(seq_along(valid_indices), sample_n))
    sample_indices <- valid_indices[sample_pos]

    x <- feature_mat[valid_indices, , drop = FALSE]
    x[!is.finite(x)] <- 0
    if (spatial_weight > 0) {
      ij <- arrayInd(valid_indices, .dim = c(n_row, n_col))
      x_sd <- stats::sd(ij[, 2])
      y_sd <- stats::sd(ij[, 1])
      if (!is.finite(x_sd) || x_sd <= 0) x_sd <- 1
      if (!is.finite(y_sd) || y_sd <= 0) y_sd <- 1
      xy <- cbind(
        x = (ij[, 2] - mean(ij[, 2])) / x_sd,
        y = (ij[, 1] - mean(ij[, 1])) / y_sd
      )
      x <- cbind(x, spatial_weight * xy)
    }

    sample_features <- x[sample_pos, , drop = FALSE]
    fit <- .sparse_ward_cluster_matrix(
      features = sample_features,
      Ncomp = Ncomp,
      knn_k = knn_k,
      auto_k = FALSE,
      verbose = isTRUE(list(...)$verbose)
    )
    sample_labels <- fit$labels
    label_levels <- sort(unique(sample_labels))
    centers <- rowsum(sample_features, sample_labels) /
      as.vector(table(factor(sample_labels, levels = label_levels)))
    centers <- centers[match(label_levels, as.integer(rownames(centers))), , drop = FALSE]

    assigned <- integer(nrow(x))
    center_norm <- rowSums(centers^2)
    block <- 5000L
    for (start in seq(1L, nrow(x), by = block)) {
      idx <- start:min(nrow(x), start + block - 1L)
      d <- sweep(-2 * tcrossprod(x[idx, , drop = FALSE], centers), 2, center_norm, "+")
      assigned[idx] <- label_levels[max.col(-d, ties.method = "first")]
    }

    cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
    cluster_map[valid_indices] <- assigned
    cluster_snr <- .compute_cluster_snr(
      clusters = assigned,
      signal_valid = rowSums(abs(feature_mat[valid_indices, , drop = FALSE]), na.rm = TRUE),
      noise_valid = sqrt(pmax(rowSums(abs(feature_mat[valid_indices, , drop = FALSE]), na.rm = TRUE), .Machine$double.eps))
    )
    out <- list(
      cluster_map = cluster_map,
      header = features$original_cube$hdr,
      axDat = features$original_cube$axDat,
      cluster_snr = cluster_snr,
      Ncomp = length(unique(assigned)),
      snr_grid = NULL,
      requested_Ncomp = as.integer(Ncomp),
      original_cube = features$original_cube,
      backend = "sparse_ward_sampled_projection",
      backend_info = list(
        algorithm = "sparse_ward_sampled_projection",
        knn_k = fit$knn_k,
        requested_Ncomp = fit$requested_Ncomp,
        actual_Ncomp = length(unique(assigned)),
        sampled_pixels = sample_n,
        valid_pixels = length(valid_indices),
        disconnected = fit$disconnected,
        feature_scale = "emission_line_columns",
        spatial_weight = spatial_weight,
        valid_mode = valid_mode
      )
    )
  }

  out$original_cube <- features$original_cube
  out$axDat <- features$original_cube$axDat
  out$header <- features$original_cube$hdr
  out$emission_line_features <- list(
    mode = feature_mode,
    redshift = redshift,
    lines = features$feature_info,
    line_window_kms = line_window_kms,
    continuum_inner_kms = continuum_inner_kms,
    continuum_outer_kms = continuum_outer_kms,
    profile_bins = profile_bins,
    continuum_subtract = continuum_subtract,
    positive_only = positive_only
  )
  out$backend_info$feature_family <- "emission_lines"
  out$backend_info$emission_feature_mode <- feature_mode
  out$backend_info$max_pixels <- max_pixels
  out$backend_info$sampled_projection <- sampled
  class(out) <- c("capivara_emission_line_segmentation", class(out))
  out
}
