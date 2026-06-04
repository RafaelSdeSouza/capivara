#' Build an adaptive multi-band/spatial support mask
#'
#' Builds a foreground support from band-wise or wavelength-slice detection
#' evidence rather than from a single white-light image. This is useful for IFU
#' cubes where the source is persistent over many wavelength channels but the
#' background may vary spatially or spectrally.
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array.
#' @param bands Wavelength slices used to build the support. May be \code{NULL}
#'   for all slices, a numeric vector of indices, or a character vector matching
#'   cube dimnames along the third dimension.
#' @param transform Per-slice image transform before robust sky normalization.
#'   Options are \code{"none"}, \code{"asinh"}, \code{"signed_log1p"},
#'   \code{"log1p"}, \code{"copula_uniform"}, and \code{"copula_gaussian"}.
#' @param sky_method Pixels used for robust sky/noise estimation. \code{"border"}
#'   uses image borders; \code{"all"} uses all finite pixels.
#' @param border_fraction Fraction of rows/columns used as image border when
#'   \code{sky_method = "border"}.
#' @param z_threshold Per-slice detection threshold in robust sigma units.
#' @param min_band_persistence Minimum number of selected slices above
#'   \code{z_threshold}. If \code{NULL}, defaults to two slices when possible.
#' @param single_band_z Optional single-slice rescue threshold. Use \code{Inf}
#'   to disable.
#' @param smooth_sigma Optional Gaussian smoothing applied to each significance
#'   map before thresholding. Requires \code{EBImage} or \code{imager}.
#'
#' @return A list with \code{collapsed}, \code{reconstruction}, \code{mask},
#'   \code{evidence}, \code{evidence_count}, \code{z_maps}, \code{band_stats},
#'   and metadata.
#' @export
build_adaptive_support <- function(input,
                                   bands = NULL,
                                   transform = c(
                                     "none",
                                     "asinh",
                                     "signed_log1p",
                                     "log1p",
                                     "copula_uniform",
                                     "copula_gaussian"
                                   ),
                                   sky_method = c("border", "all"),
                                   border_fraction = 0.10,
                                   z_threshold = 3,
                                   min_band_persistence = NULL,
                                   single_band_z = Inf,
                                   smooth_sigma = 0) {
  transform <- match.arg(transform)
  sky_method <- match.arg(sky_method)

  cube <- .as_cubedat(input)$imDat
  stopifnot(is.array(cube), length(dim(cube)) == 3L)

  band_ids <- .resolve_adaptive_bands(cube, bands)
  n_bands <- length(band_ids)

  if (is.null(min_band_persistence)) {
    min_band_persistence <- min(2L, n_bands)
  }
  min_band_persistence <- as.integer(min_band_persistence)
  if (!is.finite(min_band_persistence) || min_band_persistence < 1L) {
    min_band_persistence <- 1L
  }
  min_band_persistence <- min(min_band_persistence, n_bands)

  z_threshold <- as.numeric(z_threshold)
  if (!is.finite(z_threshold)) z_threshold <- 3
  single_band_z <- as.numeric(single_band_z)
  if (!is.finite(single_band_z)) single_band_z <- Inf

  z_maps <- array(NA_real_, dim = c(dim(cube)[1:2], n_bands))
  active <- array(FALSE, dim = dim(z_maps))
  band_stats <- vector("list", n_bands)

  for (i in seq_along(band_ids)) {
    band_index <- band_ids[i]
    image <- .adaptive_transform_image(cube[, , band_index], transform = transform)
    sky_values <- .adaptive_sky_values(
      image,
      method = sky_method,
      border_fraction = border_fraction
    )
    stats <- .adaptive_robust_location_scale(sky_values)

    z <- (image - stats$center) / stats$scale
    z[!is.finite(z)] <- NA_real_
    if (smooth_sigma > 0) {
      z <- .adaptive_smooth_image(z, sigma = smooth_sigma)
    }

    z_maps[, , i] <- z
    active[, , i] <- is.finite(z) & z >= z_threshold
    band_stats[[i]] <- data.frame(
      band_index = band_index,
      center = stats$center,
      scale = stats$scale,
      finite_sky_pixels = stats$n,
      stringsAsFactors = FALSE
    )
  }

  evidence_count <- apply(active, c(1L, 2L), sum, na.rm = TRUE)
  positive_z <- pmax(z_maps, 0)
  positive_z[!is.finite(positive_z)] <- 0
  evidence <- apply(positive_z, c(1L, 2L), mean, na.rm = TRUE)
  max_z <- apply(z_maps, c(1L, 2L), max, na.rm = TRUE)
  max_z[!is.finite(max_z)] <- NA_real_

  mask <- evidence_count >= min_band_persistence
  if (is.finite(single_band_z)) {
    mask <- mask | (is.finite(max_z) & max_z >= single_band_z)
  }
  mask[!is.finite(evidence)] <- FALSE

  list(
    collapsed = evidence,
    decomposition = NULL,
    reconstruction = evidence,
    mask = mask,
    evidence = evidence,
    evidence_count = evidence_count,
    z_maps = z_maps,
    band_stats = do.call(rbind, band_stats),
    bands = band_ids,
    transform = transform,
    sky_method = sky_method,
    border_fraction = border_fraction,
    z_threshold = z_threshold,
    min_band_persistence = min_band_persistence,
    single_band_z = single_band_z,
    smooth_sigma = smooth_sigma
  )
}

.resolve_adaptive_bands <- function(cube, bands) {
  n_wave <- dim(cube)[3]
  if (is.null(bands)) {
    return(seq_len(n_wave))
  }
  if (is.numeric(bands)) {
    bands <- as.integer(bands)
    if (!length(bands) || any(!is.finite(bands)) ||
        any(bands < 1L) || any(bands > n_wave)) {
      stop("Numeric `bands` must be valid indices along the third cube dimension.")
    }
    return(unique(bands))
  }
  if (is.character(bands)) {
    band_names <- dimnames(cube)[[3]]
    if (is.null(band_names)) {
      stop("Character `bands` requires dimnames along the third cube dimension.")
    }
    idx <- match(bands, band_names)
    if (any(is.na(idx))) {
      missing <- paste(bands[is.na(idx)], collapse = ", ")
      stop("Could not match adaptive-support bands: ", missing)
    }
    return(unique(idx))
  }
  stop("`bands` must be NULL, numeric indices, or character band names.")
}

.adaptive_transform_image <- function(image, transform) {
  out <- as.matrix(image)
  if (transform == "none") return(out)
  if (transform == "asinh") return(asinh(out))
  if (transform == "signed_log1p") return(sign(out) * log1p(abs(out)))
  if (transform == "log1p") {
    finite <- out[is.finite(out)]
    floor <- min(finite, na.rm = TRUE)
    if (!is.finite(floor)) floor <- 0
    return(log1p(pmax(out - floor, 0)))
  }
  if (transform %in% c("copula_uniform", "copula_gaussian")) {
    flat <- as.numeric(out)
    ok <- is.finite(flat)
    z <- rep(NA_real_, length(flat))
    if (any(ok)) {
      ranks <- rank(flat[ok], ties.method = "average")
      u <- (ranks - 0.5) / sum(ok)
      if (transform == "copula_gaussian") u <- stats::qnorm(u)
      z[ok] <- u
    }
    return(matrix(z, nrow = nrow(out), ncol = ncol(out)))
  }
  stop("Unknown adaptive-support transform: ", transform)
}

.adaptive_sky_values <- function(image,
                                 method = c("border", "all"),
                                 border_fraction = 0.10) {
  method <- match.arg(method)
  image <- as.matrix(image)
  if (method == "all") return(image[is.finite(image)])

  border_fraction <- as.numeric(border_fraction)
  if (!is.finite(border_fraction) || border_fraction <= 0) {
    border_fraction <- 0.10
  }
  border_fraction <- min(border_fraction, 0.50)

  nr <- nrow(image)
  nc <- ncol(image)
  br <- max(1L, floor(border_fraction * nr))
  bc <- max(1L, floor(border_fraction * nc))

  border <- matrix(FALSE, nrow = nr, ncol = nc)
  border[seq_len(br), ] <- TRUE
  border[(nr - br + 1L):nr, ] <- TRUE
  border[, seq_len(bc)] <- TRUE
  border[, (nc - bc + 1L):nc] <- TRUE

  image[border & is.finite(image)]
}

.adaptive_robust_location_scale <- function(values) {
  values <- values[is.finite(values)]
  n <- length(values)
  if (!n) return(list(center = 0, scale = 1, n = 0L))

  center <- stats::median(values, na.rm = TRUE)
  scale <- stats::mad(values, center = center, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) {
    q <- stats::quantile(values, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    scale <- diff(q) / 1.349
  }
  if (!is.finite(scale) || scale <= 0) scale <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(scale) || scale <= 0) scale <- 1

  list(center = center, scale = scale, n = n)
}

.adaptive_smooth_image <- function(image, sigma) {
  if (requireNamespace("EBImage", quietly = TRUE)) {
    out <- EBImage::gblur(EBImage::Image(image), sigma = sigma)
    return(as.matrix(EBImage::imageData(out)))
  }
  if (requireNamespace("imager", quietly = TRUE)) {
    out <- imager::isoblur(imager::as.cimg(image), sigma = sigma)
    return(as.matrix(out[, , 1, 1]))
  }
  warning("`smooth_sigma` was requested, but neither EBImage nor imager is installed; skipping smoothing.")
  image
}

#' Reject background-like regions after segmentation
#'
#' @param segmentation A result returned by \code{\link{segment}} or
#'   \code{\link{segment_large}}.
#' @param evidence Optional 2-D evidence map. If \code{NULL}, the function tries
#'   to use \code{segmentation$support_info$evidence}.
#' @param evidence_count Optional 2-D persistence map. If \code{NULL}, the
#'   function tries to use \code{segmentation$support_info$evidence_count}.
#' @param min_median_evidence Minimum median evidence required for a region. If
#'   \code{NULL}, uses \code{low_evidence_quantile} of region median evidence.
#' @param min_median_persistence Minimum median persistence required for a
#'   region.
#' @param low_evidence_quantile Quantile used to define low-evidence regions
#'   when \code{min_median_evidence = NULL}.
#' @param reject_edge_low_evidence Logical; reject edge-connected regions with
#'   low median evidence.
#' @param reject_large_low_evidence Logical; reject large regions with low
#'   median evidence.
#' @param large_area_quantile Quantile of region areas used to identify large
#'   regions.
#'
#' @return The input segmentation with an updated \code{cluster_map} and
#'   appended diagnostics.
#' @export
reject_background_regions <- function(segmentation,
                                      evidence = NULL,
                                      evidence_count = NULL,
                                      min_median_evidence = NULL,
                                      min_median_persistence = 1,
                                      low_evidence_quantile = 0.20,
                                      reject_edge_low_evidence = TRUE,
                                      reject_large_low_evidence = TRUE,
                                      large_area_quantile = 0.75) {
  if (!is.list(segmentation) || is.null(segmentation$cluster_map)) {
    stop("`segmentation` must be a Capivara segmentation result with `cluster_map`.")
  }
  labels <- segmentation$cluster_map
  if (is.null(evidence)) evidence <- segmentation$support_info$evidence
  if (is.null(evidence_count)) evidence_count <- segmentation$support_info$evidence_count
  if (is.null(evidence) || is.null(evidence_count)) {
    stop("`evidence` and `evidence_count` are required unless adaptive support was used.")
  }
  if (!identical(dim(labels), dim(evidence)) ||
      !identical(dim(labels), dim(evidence_count))) {
    stop("`evidence` and `evidence_count` must have the same dimensions as `cluster_map`.")
  }

  ids <- sort(unique(stats::na.omit(as.vector(labels))))
  ids <- ids[ids > 0]
  if (!length(ids)) {
    segmentation$background_diagnostics <- data.frame()
    segmentation$background_regions <- integer()
    return(segmentation)
  }

  nr <- nrow(labels)
  nc <- ncol(labels)
  diagnostics <- data.frame(
    region = ids,
    area = vapply(ids, function(id) sum(labels == id, na.rm = TRUE), numeric(1)),
    edge_connected = vapply(ids, function(id) {
      m <- labels == id
      any(m[1, ], na.rm = TRUE) || any(m[nr, ], na.rm = TRUE) ||
        any(m[, 1], na.rm = TRUE) || any(m[, nc], na.rm = TRUE)
    }, logical(1)),
    median_evidence = vapply(ids, function(id) {
      stats::median(evidence[labels == id], na.rm = TRUE)
    }, numeric(1)),
    median_persistence = vapply(ids, function(id) {
      stats::median(evidence_count[labels == id], na.rm = TRUE)
    }, numeric(1))
  )

  if (is.null(min_median_evidence)) {
    low_evidence_quantile <- min(max(as.numeric(low_evidence_quantile), 0), 1)
    if (!is.finite(low_evidence_quantile)) low_evidence_quantile <- 0.20
    min_median_evidence <- stats::quantile(
      diagnostics$median_evidence,
      probs = low_evidence_quantile,
      na.rm = TRUE,
      names = FALSE
    )
  }
  min_median_persistence <- as.numeric(min_median_persistence)
  if (!is.finite(min_median_persistence)) min_median_persistence <- 1

  large_area_quantile <- min(max(as.numeric(large_area_quantile), 0), 1)
  if (!is.finite(large_area_quantile)) large_area_quantile <- 0.75
  area_cut <- stats::quantile(
    diagnostics$area,
    probs = large_area_quantile,
    na.rm = TRUE,
    names = FALSE
  )

  low_evidence <- diagnostics$median_evidence <= min_median_evidence
  low_persistence <- diagnostics$median_persistence < min_median_persistence
  edge_low <- isTRUE(reject_edge_low_evidence) &
    diagnostics$edge_connected & low_evidence
  large_low <- isTRUE(reject_large_low_evidence) &
    diagnostics$area >= area_cut & low_evidence

  diagnostics$background_like <- low_persistence | edge_low | large_low
  background_ids <- diagnostics$region[diagnostics$background_like]

  cleaned <- labels
  cleaned[cleaned %in% background_ids] <- NA_integer_

  segmentation$cluster_map <- cleaned
  segmentation$background_regions <- background_ids
  segmentation$background_diagnostics <- diagnostics
  segmentation
}
