.as_cubedat <- function(input) {
  if (is.list(input) && !is.null(input$imDat)) {
    return(input)
  }

  list(imDat = input, hdr = NULL, axDat = NULL)
}

.compute_signal_noise <- function(IFU2D) {
  signal <- rowSums(IFU2D, na.rm = TRUE)
  signal[is.na(signal) | signal <= 0] <- 0
  noise <- sqrt(signal)
  noise[is.na(noise) | noise == 0] <- Inf

  list(signal = signal, noise = noise)
}

.wavelength_axis <- function(cubedat, n_wave) {
  if (!is.null(cubedat$axDat)) {
    wavelengths <- tryCatch(
      FITSio::axVec(3, cubedat$axDat),
      error = function(e) NULL
    )
    if (!is.null(wavelengths) && length(wavelengths) == n_wave) {
      return(wavelengths)
    }
  }

  seq_len(n_wave)
}

.wavelength_range_index <- function(cubedat,
                                    n_wave,
                                    wavelength_range = NULL,
                                    arg_name = "wavelength_range") {
  if (is.null(wavelength_range)) {
    return(seq_len(n_wave))
  }

  if (length(wavelength_range) != 2) {
    stop("`", arg_name, "` must have length 2.")
  }

  wavelengths <- .wavelength_axis(cubedat, n_wave)
  wave_idx <- which(wavelengths >= min(wavelength_range) & wavelengths <= max(wavelength_range))

  if (!length(wave_idx)) {
    stop("No wavelengths fall inside `", arg_name, "`.")
  }

  wave_idx
}

.subset_cubedat_wavelength_range <- function(cubedat, feature_wavelength_range = NULL) {
  cubedat <- .as_cubedat(cubedat)
  cube <- cubedat$imDat

  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("`input$imDat` must be a 3D array with dimensions (n_row, n_col, n_wave).")
  }

  n_wave <- dim(cube)[3]
  wave_idx <- .wavelength_range_index(
    cubedat = cubedat,
    n_wave = n_wave,
    wavelength_range = feature_wavelength_range,
    arg_name = "feature_wavelength_range"
  )

  wavelengths <- .wavelength_axis(cubedat, n_wave)

  if (is.null(feature_wavelength_range)) {
    return(list(
      cubedat = cubedat,
      wave_idx = wave_idx,
      selected_wavelengths = wavelengths
    ))
  }

  out <- cubedat
  out$imDat <- cube[, , wave_idx, drop = FALSE]

  if (!is.null(out$axDat) && is.data.frame(out$axDat) && nrow(out$axDat) >= 3L) {
    if ("crval" %in% names(out$axDat)) {
      out$axDat[3, "crval"] <- wavelengths[wave_idx[1]]
    }
    if ("crpix" %in% names(out$axDat)) {
      out$axDat[3, "crpix"] <- 1
    }
    if ("cdelt" %in% names(out$axDat) && length(wave_idx) > 1L) {
      out$axDat[3, "cdelt"] <- stats::median(diff(wavelengths[wave_idx]), na.rm = TRUE)
    }
    if ("len" %in% names(out$axDat)) {
      out$axDat[3, "len"] <- length(wave_idx)
    }
  }

  list(
    cubedat = out,
    wave_idx = wave_idx,
    selected_wavelengths = wavelengths[wave_idx]
  )
}

.scale_rows <- function(X, scale_fn, na_to_zero = FALSE) {
  X <- as.matrix(X)

  if (!nrow(X)) {
    return(matrix(numeric(0), nrow = 0L, ncol = ncol(X)))
  }

  scaled <- t(vapply(
    seq_len(nrow(X)),
    function(i) as.numeric(scale_fn(X[i, ])),
    numeric(ncol(X))
  ))

  if (na_to_zero) {
    scaled[!is.finite(scaled)] <- 0
  }

  scaled
}

.compute_cluster_snr <- function(clusters, signal_valid, noise_valid) {
  vapply(sort(unique(clusters)), function(cluster_id) {
    idx <- which(clusters == cluster_id)
    sum(signal_valid[idx]) / sqrt(sum(noise_valid[idx]^2))
  }, numeric(1))
}

.segment_core <- function(input,
                          Ncomp = 5,
                          redshift = 0,
                          scale_fn = median_scale,
                          na_to_zero = FALSE,
                          return_details = FALSE) {
  cubedat <- .as_cubedat(input)
  cube <- cubedat$imDat

  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("`input$imDat` must be a 3D array with dimensions (n_row, n_col, n_wave).")
  }

  n_row <- dim(cube)[1]
  n_col <- dim(cube)[2]
  IFU2D <- cube_to_matrix(cubedat)

  sn <- .compute_signal_noise(IFU2D)
  valid_indices <- which(sn$signal > 0)

  if (!length(valid_indices)) {
    stop("No valid pixels after filtering.")
  }
  if (!is.null(Ncomp) && length(valid_indices) < Ncomp) {
    stop("`Ncomp` is larger than the number of valid pixels.")
  }

  IFU2D_valid <- IFU2D[valid_indices, , drop = FALSE]
  signal_valid <- sn$signal[valid_indices]
  noise_valid <- sn$noise[valid_indices]

  scaled_data <- .scale_rows(IFU2D_valid, scale_fn = scale_fn, na_to_zero = na_to_zero)
  distance_matrix <- torch_dist(scaled_data, p = 1)
  hc <- fastcluster::hclust(distance_matrix, method = "ward.D2")

  if (is.null(Ncomp)) {
    out <- list(
      header = cubedat$hdr,
      axDat = cubedat$axDat,
      original_cube = cubedat
    )
  } else {
    clusters <- stats::cutree(hc, k = Ncomp)
    cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
    cluster_map[valid_indices] <- clusters

    out <- list(
      cluster_map = cluster_map,
      header = cubedat$hdr,
      axDat = cubedat$axDat,
      cluster_snr = .compute_cluster_snr(clusters, signal_valid, noise_valid),
      original_cube = cubedat
    )
  }

  if (return_details) {
    out <- c(out, list(
      hclust = hc,
      valid_indices = valid_indices,
      signal = sn$signal,
      noise = sn$noise
    ))
  }

  out
}
