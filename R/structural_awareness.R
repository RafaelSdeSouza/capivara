#' Detect a spatial support mask
#'
#' This is the first structural-awareness layer in \pkg{capivara}. It decides
#' where signal exists before any semantic structure labels are attempted.
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array.
#' @param method Support detector. Currently \code{"starlet"} and
#'   \code{"adaptive"} are supported.
#' @param ... Additional arguments passed to the selected support detector.
#'
#' @return A support object with at least \code{collapsed},
#'   \code{reconstruction}, and \code{mask}.
#' @export
detect_support <- function(input, method = c("starlet", "adaptive"), ...) {
  method <- match.arg(method)

  if (method == "starlet") {
    return(build_starlet_mask(input, ...))
  }

  build_adaptive_support(input, ...)
}

.structure_robust_norm <- function(x, probs = c(0.02, 0.98)) {
  y <- x
  y[!is.finite(y)] <- NA_real_
  q <- stats::quantile(y, probs, na.rm = TRUE)
  if (!all(is.finite(q)) || q[2] <= q[1]) {
    y[is.finite(y)] <- 0
    return(y)
  }
  y <- pmin(pmax((y - q[1]) / (q[2] - q[1]), 0), 1)
  y[!is.finite(y)] <- 0
  y
}

.structure_reflect_pad_vec <- function(v, pad) {
  n <- length(v)
  if (pad <= 0L || n == 0L) return(v)
  pad <- max(0L, min(pad, n - 1L))
  if (pad == 0L) return(v)
  c(rev(v[seq_len(pad)]), v, rev(v[seq.int(n - pad + 1L, n)]))
}

.structure_conv2 <- function(img, kernel) {
  nr <- nrow(img)
  nc <- ncol(img)
  kr <- nrow(kernel)
  kc <- ncol(kernel)
  rr <- floor(kr / 2)
  cc <- floor(kc / 2)
  out <- matrix(0, nr, nc)
  z <- img
  z[!is.finite(z)] <- 0

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      ii <- (i - rr):(i + rr)
      jj <- (j - cc):(j + cc)
      ok_i <- ii >= 1 & ii <= nr
      ok_j <- jj >= 1 & jj <= nc
      patch <- matrix(0, kr, kc)
      patch[ok_i, ok_j] <- z[ii[ok_i], jj[ok_j], drop = FALSE]
      out[i, j] <- sum(patch * kernel)
    }
  }

  out
}

.structure_line_kernel <- function(theta,
                                   radius = 5L,
                                   sigma_long = 4,
                                   sigma_short = 1) {
  g <- seq(-radius, radius)
  xx <- matrix(rep(g, each = length(g)), nrow = length(g))
  yy <- matrix(rep(g, times = length(g)), nrow = length(g))
  u <- xx * cos(theta) + yy * sin(theta)
  v <- -xx * sin(theta) + yy * cos(theta)
  core <- exp(-0.5 * ((u / sigma_long)^2 + (v / sigma_short)^2))
  bg <- exp(-0.5 * ((u / sigma_long)^2 + (v / (2.8 * sigma_short))^2))
  k <- core / sum(core) - bg / sum(bg)
  k - mean(k)
}

.structure_ridge_score <- function(img,
                                   angles = seq(0, pi - pi / 18, length.out = 18),
                                   radii = c(3L, 5L, 7L)) {
  score <- matrix(-Inf, nrow(img), ncol(img))
  angle_map <- matrix(NA_real_, nrow(img), ncol(img))

  for (r in radii) {
    rr <- min(r, max(2L, floor(min(dim(img)) / 7)))
    for (a in angles) {
      resp <- .structure_conv2(
        img,
        .structure_line_kernel(
          a,
          radius = rr,
          sigma_long = max(2, 0.85 * rr),
          sigma_short = max(0.85, 0.22 * rr)
        )
      )
      take <- resp > score
      score[take] <- resp[take]
      angle_map[take] <- a
    }
  }

  score[!is.finite(score)] <- NA_real_
  list(score = score, angle = angle_map)
}

.structure_neighbor_count <- function(mask) {
  out <- matrix(0L, nrow(mask), ncol(mask))
  for (i in seq_len(nrow(mask))) {
    for (j in seq_len(ncol(mask))) {
      out[i, j] <- sum(
        mask[max(1, i - 1):min(nrow(mask), i + 1),
             max(1, j - 1):min(ncol(mask), j + 1)],
        na.rm = TRUE
      )
    }
  }
  out
}

.structure_components <- function(mask) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  lab <- matrix(0L, nr, nc)
  sizes <- integer()
  cur <- 0L

  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      if (!isTRUE(mask[i, j]) || lab[i, j] != 0L) next

      cur <- cur + 1L
      qi <- i
      qj <- j
      lab[i, j] <- cur
      head <- 1L
      sz <- 0L

      while (head <= length(qi)) {
        a <- qi[head]
        b <- qj[head]
        head <- head + 1L
        sz <- sz + 1L

        for (d in list(c(1L, 0L), c(-1L, 0L), c(0L, 1L), c(0L, -1L))) {
          aa <- a + d[1]
          bb <- b + d[2]
          if (aa >= 1L && aa <= nr && bb >= 1L && bb <= nc &&
              isTRUE(mask[aa, bb]) && lab[aa, bb] == 0L) {
            lab[aa, bb] <- cur
            qi <- c(qi, aa)
            qj <- c(qj, bb)
          }
        }
      }

      sizes[cur] <- sz
    }
  }

  list(label = lab, sizes = sizes)
}

.structure_drop_small <- function(mask, min_area = 20L) {
  cc <- .structure_components(mask)
  keep <- which(cc$sizes >= min_area)
  if (!length(keep)) return(mask & FALSE)
  out <- cc$label %in% keep
  dim(out) <- dim(mask)
  out
}

.structure_dilate <- function(mask, radius = 1L) {
  out <- mask & FALSE
  idx <- which(mask, arr.ind = TRUE)
  if (!nrow(idx)) return(out)

  for (k in seq_len(nrow(idx))) {
    out[max(1, idx[k, 1] - radius):min(nrow(mask), idx[k, 1] + radius),
        max(1, idx[k, 2] - radius):min(ncol(mask), idx[k, 2] + radius)] <- TRUE
  }

  out
}

.structure_largest_center_component <- function(mask, center, radius = 5L) {
  cc <- .structure_components(mask)
  if (!length(cc$sizes)) return(mask & FALSE)

  rr <- sqrt((row(mask) - center[1])^2 + (col(mask) - center[2])^2)
  labs <- unique(cc$label[mask & rr <= radius])
  labs <- labs[labs > 0]
  if (!length(labs)) labs <- which.max(cc$sizes)

  keep <- labs[which.max(cc$sizes[labs])]
  out <- cc$label == keep
  dim(out) <- dim(mask)
  out
}

.structure_otsu_threshold <- function(x, bins = 128L) {
  x <- x[is.finite(x)]
  x <- x[x >= 0 & x <= 1]
  if (length(x) < 16L) return(NA_real_)

  h <- hist(x, breaks = seq(0, 1, length.out = bins + 1L), plot = FALSE)
  prob <- h$counts / sum(h$counts)
  omega <- cumsum(prob)
  mu <- cumsum(prob * h$mids)
  mu_t <- sum(prob * h$mids)
  denom <- omega * (1 - omega)
  sigma_b <- (mu_t * omega - mu)^2 / pmax(denom, .Machine$double.eps)
  sigma_b[denom <= 0] <- NA_real_

  h$mids[which.max(sigma_b)]
}

.structure_hysteresis_mask <- function(score,
                                       support_mask,
                                       high_cut,
                                       low_quantile = 0.2) {
  vals <- score[support_mask & is.finite(score)]
  if (length(vals) < 16L || !is.finite(high_cut)) {
    return(support_mask & FALSE)
  }

  low_cut <- stats::quantile(vals, low_quantile, na.rm = TRUE)
  low_cut <- min(low_cut, high_cut, na.rm = TRUE)
  seed <- support_mask & score >= high_cut
  candidate <- support_mask & score >= low_cut
  cc <- .structure_components(candidate)
  seed_labs <- unique(cc$label[seed])
  seed_labs <- seed_labs[seed_labs > 0]
  if (!length(seed_labs)) return(seed)

  out <- cc$label %in% seed_labs
  dim(out) <- dim(score)
  out
}

.structure_fit_component <- function(mask, score, center) {
  pts <- which(mask, arr.ind = TRUE)
  if (nrow(pts) < 8L) return(NULL)

  wt <- pmax(score[pts], 0.05)
  x <- cbind(x = pts[, 2], y = pts[, 1])
  centroid <- colSums(x * wt) / sum(wt)
  xc <- sweep(x, 2, centroid, "-")
  cov <- crossprod(xc, xc * wt) / sum(wt)
  eig <- eigen(cov)

  semi_major <- sqrt(max(eig$values, na.rm = TRUE)) * 2
  semi_minor <- sqrt(min(eig$values, na.rm = TRUE)) * 2
  axis_ratio <- semi_minor / semi_major
  r_pix <- sqrt((pts[, 1] - center[1])^2 + (pts[, 2] - center[2])^2)
  theta <- atan2(pts[, 1] - center[1], pts[, 2] - center[2])
  theta_bin <- floor(((theta + pi) / (2 * pi)) * 24)
  theta_bin <- pmin(pmax(theta_bin, 0L), 23L)

  data.frame(
    area_px = nrow(pts),
    x_centroid = centroid[1],
    y_centroid = centroid[2],
    r_centroid_px = sqrt((centroid[2] - center[1])^2 + (centroid[1] - center[2])^2),
    median_radius_px = stats::median(r_pix, na.rm = TRUE),
    radius_mad_px = stats::mad(r_pix, constant = 1.4826, na.rm = TRUE),
    radius_cv = stats::mad(r_pix, constant = 1.4826, na.rm = TRUE) /
      max(stats::median(r_pix, na.rm = TRUE), 1e-6),
    azimuth_coverage = length(unique(theta_bin)) / 24,
    semi_major_px = semi_major,
    semi_minor_px = semi_minor,
    axis_ratio = axis_ratio,
    ellipticity = 1 - axis_ratio,
    pa_deg = atan2(eig$vectors[2, 1], eig$vectors[1, 1]) * 180 / pi,
    mean_score = mean(score[pts], na.rm = TRUE),
    max_score = max(score[pts], na.rm = TRUE)
  )
}

#' Score structural morphology inside a support mask
#'
#' Computes neutral structure-score maps from a cube. The scores are not final
#' semantic labels; they are evidence maps that can later be thresholded,
#' catalogued, or used as extra segmentation features.
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array.
#' @param support Optional support object from \code{\link{detect_support}}.
#' @param support_method Method passed to \code{\link{detect_support}} when
#'   \code{support} is \code{NULL}.
#' @param structure_scales Starlet scales used for the structural score.
#' @param ridge_weight Weight of the ridge score in the combined score.
#' @param starlet_weight Weight of the starlet score in the combined score.
#' @param starlet_args Extra arguments passed to \code{\link{detect_support}}.
#'
#' @return A \code{capivara_structure_scores} object.
#' @export
score_structures <- function(input,
                             support = NULL,
                             support_method = "starlet",
                             structure_scales = 3:4,
                             ridge_weight = 0.35,
                             starlet_weight = 0.65,
                             starlet_args = list()) {
  if (is.null(support)) {
    support <- do.call(
      detect_support,
      c(list(input = input, method = support_method), starlet_args)
    )
  }

  collapsed <- .structure_robust_norm(support$collapsed)
  center <- which(collapsed == max(collapsed, na.rm = TRUE), arr.ind = TRUE)[1, ]

  if (!is.null(support$decomposition)) {
    starlet_rec <- starlet_reconstruct(
      support$decomposition,
      keep_scales = structure_scales,
      include_coarse = FALSE,
      denoise_k = 0,
      mode = "soft"
    )
  } else {
    starlet_rec <- support$reconstruction
  }

  starlet_score <- .structure_robust_norm(pmax(starlet_rec, 0))
  ridge <- .structure_ridge_score(collapsed)
  ridge_score <- .structure_robust_norm(pmax(ridge$score, 0))

  support_mask <- support$mask & is.finite(collapsed)
  pos <- collapsed[collapsed > 0 & is.finite(collapsed)]
  if (length(pos) >= 16L) {
    raw <- collapsed >= stats::quantile(pos, 0.18, na.rm = TRUE)
    footprint <- .structure_largest_center_component(raw, center, radius = 5L)
    footprint <- .structure_dilate(footprint, 1L) &
      collapsed >= stats::quantile(pos, 0.08, na.rm = TRUE)
    footprint <- .structure_drop_small(
      footprint & .structure_neighbor_count(footprint) >= 3L,
      max(20L, floor(sum(footprint) * 0.004))
    )
    support_mask <- support_mask & footprint
  }

  combined <- .structure_robust_norm(
    starlet_weight * starlet_score + ridge_weight * ridge_score
  )
  combined[!support_mask] <- NA_real_

  smooth_score <- .structure_robust_norm(pmax(collapsed - ridge_score, 0))
  r <- sqrt((row(collapsed) - center[1])^2 + (col(collapsed) - center[2])^2)
  r95 <- stats::quantile(r[support_mask], 0.95, na.rm = TRUE)
  if (!is.finite(r95) || r95 <= 0) r95 <- max(r[support_mask], na.rm = TRUE)
  central_weight <- exp(-0.5 * (r / max(0.22 * r95, 1))^2)
  outer_weight <- pmin(r / max(0.25 * r95, 1), 1)
  annular_weight <- exp(-0.5 * ((r - 0.55 * r95) / max(0.22 * r95, 1))^2)

  central_smooth_score <- .structure_robust_norm(
    central_weight * smooth_score * (1 - ridge_score)
  )
  central_elongated_score <- .structure_robust_norm(
    central_weight * ridge_score * combined
  )
  outer_ridge_score <- .structure_robust_norm(
    outer_weight * ridge_score * combined
  )
  annular_score <- .structure_robust_norm(annular_weight * combined)

  orientation_cos <- 0.5 + 0.5 * cos(2 * ridge$angle)
  orientation_sin <- 0.5 + 0.5 * sin(2 * ridge$angle)
  orientation_cos[!is.finite(orientation_cos)] <- 0
  orientation_sin[!is.finite(orientation_sin)] <- 0

  out <- list(
    support = support,
    center = center,
    support_mask = support_mask,
    maps = list(
      collapsed = collapsed,
      starlet_score = starlet_score,
      ridge_score = ridge_score,
      combined_structure_score = combined,
      smooth_score = smooth_score,
      central_smooth_score = central_smooth_score,
      central_elongated_score = central_elongated_score,
      outer_ridge_score = outer_ridge_score,
      annular_score = annular_score,
      orientation_cos2 = orientation_cos,
      orientation_sin2 = orientation_sin
    ),
    parameters = list(
      structure_scales = structure_scales,
      ridge_weight = ridge_weight,
      starlet_weight = starlet_weight
    )
  )

  class(out) <- c("capivara_structure_scores", "list")
  out
}

#' Threshold structural scores into a structure mask
#'
#' @param scores A \code{capivara_structure_scores} object.
#' @param mode Thresholding mode. \code{"extended"} uses hysteresis to preserve
#'   faint connected arms and rings. \code{"strict"} keeps only high-score
#'   pixels.
#' @param score_map Name of the score map used for thresholding.
#' @param low_quantile Low hysteresis quantile for \code{mode = "extended"}.
#'
#' @return The input \code{scores} object with a \code{masks} entry.
#' @export
threshold_structures <- function(scores,
                                 mode = c("extended", "strict"),
                                 score_map = "combined_structure_score",
                                 low_quantile = 0.2) {
  mode <- match.arg(mode)
  score <- scores$maps[[score_map]]
  if (is.null(score)) stop("Unknown structure score map: ", score_map)

  vals <- score[scores$support_mask & is.finite(score)]
  if (length(vals) < 16L) stop("Too few supported pixels to threshold structures.")

  otsu <- .structure_otsu_threshold(vals)
  robust <- stats::median(vals, na.rm = TRUE) +
    0.35 * stats::mad(vals, constant = 1.4826, na.rm = TRUE)
  floor_cut <- stats::quantile(vals, 0.30, na.rm = TRUE)
  ceiling_cut <- stats::quantile(vals, 0.75, na.rm = TRUE)
  high_cut <- if (is.finite(otsu)) otsu else robust
  high_cut <- max(high_cut, floor_cut, na.rm = TRUE)
  high_cut <- min(high_cut, ceiling_cut, na.rm = TRUE)

  high <- scores$support_mask & score >= high_cut
  high <- .structure_drop_small(
    high & .structure_neighbor_count(high) >= 2L,
    max(10L, floor(sum(scores$support_mask) * 0.004))
  )

  if (mode == "strict") {
    structure <- high
  } else {
    structure <- .structure_hysteresis_mask(
      score,
      scores$support_mask,
      high_cut = high_cut,
      low_quantile = low_quantile
    )
    structure <- .structure_drop_small(
      structure & .structure_neighbor_count(structure) >= 2L,
      max(16L, floor(sum(scores$support_mask) * 0.006))
    )
  }

  intermediate_cut <- stats::quantile(vals, 0.42, na.rm = TRUE)
  intermediate <- scores$support_mask & !high & score >= intermediate_cut
  intermediate <- .structure_drop_small(
    intermediate & .structure_neighbor_count(intermediate) >= 2L,
    max(10L, floor(sum(scores$support_mask) * 0.004))
  )

  diffuse <- scores$support_mask & !structure
  diffuse <- .structure_drop_small(
    diffuse & .structure_neighbor_count(diffuse) >= 3L,
    max(20L, floor(sum(scores$support_mask) * 0.006))
  )

  tier_map <- matrix(NA_character_, nrow(score), ncol(score))
  tier_map[diffuse] <- "diffuse_support"
  tier_map[intermediate] <- "intermediate_structure"
  tier_map[high] <- "high_structure"

  scores$masks <- list(
    structure_mask = structure,
    high_structure = high,
    intermediate_structure = intermediate,
    diffuse_support = diffuse
  )
  scores$tier_map <- tier_map
  scores$threshold <- list(
    mode = mode,
    score_map = score_map,
    high_cut = as.numeric(high_cut),
    low_quantile = if (mode == "extended") low_quantile else NA_real_
  )

  scores
}

#' Build a neutral structural component catalogue
#'
#' This function deliberately avoids assigning astrophysical labels such as
#' bars, rings, spiral arms, or bulges. It records connected structural
#' components and continuous morphology scores. Scientific labels should be
#' assigned in a later, survey- or science-case-specific layer.
#'
#' @param scores A thresholded \code{capivara_structure_scores} object.
#' @param mask_name Name of the mask to catalogue.
#' @param min_area Minimum connected-component area in pixels.
#'
#' @return A list with \code{catalogue}, \code{component_map}, and
#'   \code{morphology_maps}.
#' @export
catalogue_structures <- function(scores,
                                 mask_name = "structure_mask",
                                 min_area = 12L) {
  if (is.null(scores$masks[[mask_name]])) {
    stop("`scores` does not contain mask `", mask_name, "`. Run threshold_structures() first.")
  }

  mask <- scores$masks[[mask_name]]
  combined <- scores$maps$combined_structure_score
  cc <- .structure_components(mask)
  component_map <- matrix(NA_integer_, nrow(mask), ncol(mask))
  rows <- list()
  next_id <- 0L
  support_r <- sqrt((row(mask) - scores$center[1])^2 + (col(mask) - scores$center[2])^2)
  r95 <- stats::quantile(support_r[scores$support_mask], 0.95, na.rm = TRUE)
  if (!is.finite(r95) || r95 <= 0) r95 <- max(support_r[scores$support_mask], na.rm = TRUE)

  for (lab in seq_along(cc$sizes)) {
    if (cc$sizes[lab] < min_area) next
    comp <- cc$label == lab
    fit <- .structure_fit_component(comp, combined, scores$center)
    if (is.null(fit)) next

    next_id <- next_id + 1L
    component_map[comp] <- next_id
    fit$structure_id <- next_id

    comp_idx <- which(comp, arr.ind = TRUE)
    score_mean <- function(name) {
      mean(scores$maps[[name]][comp_idx], na.rm = TRUE)
    }
    neutral_scores <- c(
      central_smooth = score_mean("central_smooth_score"),
      central_elongated = score_mean("central_elongated_score"),
      outer_ridge = score_mean("outer_ridge_score"),
      annular = score_mean("annular_score")
    )
    neutral_scores[!is.finite(neutral_scores)] <- 0

    fit$central_smooth_score <- neutral_scores[["central_smooth"]]
    fit$central_elongated_score <- neutral_scores[["central_elongated"]]
    fit$outer_ridge_score <- neutral_scores[["outer_ridge"]]
    fit$annular_score <- neutral_scores[["annular"]]
    fit$dominant_morphology_score <- names(which.max(neutral_scores))
    fit$dominance <- max(neutral_scores) -
      stats::median(neutral_scores, na.rm = TRUE)

    rows[[length(rows) + 1L]] <- fit
  }

  catalogue <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(catalogue)) {
    catalogue <- catalogue[
      c(
        "structure_id", "dominant_morphology_score", "dominance", "area_px",
        "x_centroid", "y_centroid", "r_centroid_px", "median_radius_px",
        "radius_cv", "azimuth_coverage", "semi_major_px", "semi_minor_px",
        "axis_ratio", "ellipticity", "pa_deg", "mean_score", "max_score",
        "central_smooth_score", "central_elongated_score", "outer_ridge_score",
        "annular_score"
      )
    ]
  }

  morphology_maps <- list(
    central_smooth = scores$maps$central_smooth_score,
    central_elongated = scores$maps$central_elongated_score,
    outer_ridge = scores$maps$outer_ridge_score,
    annular = scores$maps$annular_score,
    combined_structure = scores$maps$combined_structure_score,
    ridge = scores$maps$ridge_score
  )

  out <- list(
    catalogue = catalogue,
    component_map = component_map,
    morphology_maps = morphology_maps,
    scores = scores
  )
  class(out) <- c("capivara_structure_catalogue", "list")
  out
}

# Deprecated internal alias for development notebooks.
classify_structures <- function(scores,
                                mask_name = "structure_mask",
                                min_area = 12L) {
  catalogue_structures(scores, mask_name = mask_name, min_area = min_area)
}

.structure_augment_input <- function(input, scores, feature_maps, repeats) {
  cubedat <- .as_cubedat(input)
  cube <- cubedat$imDat
  maps <- scores$maps[feature_maps]
  maps <- maps[!vapply(maps, is.null, logical(1))]
  if (!length(maps) || repeats <= 0L) return(cubedat)

  extra_n <- length(maps) * repeats
  aug <- array(NA_real_, dim = c(dim(cube)[1], dim(cube)[2], dim(cube)[3] + extra_n))
  aug[, , seq_len(dim(cube)[3])] <- cube

  k <- dim(cube)[3]
  for (m in maps) {
    m[!is.finite(m)] <- 0
    m[!scores$support_mask] <- 0
    for (r in seq_len(repeats)) {
      k <- k + 1L
      aug[, , k] <- m
    }
  }

  cubedat$imDat <- aug
  cubedat
}

#' Segment a thresholded structural mask
#'
#' Runs \code{\link{segment_large}} on a structural mask, optionally appending
#' selected structure-score maps as extra features.
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array.
#' @param scores A thresholded \code{capivara_structure_scores} object.
#' @param Ncomp Number of components.
#' @param mask_name Name of the structural mask to segment.
#' @param feature_maps Structure maps appended as extra feature channels.
#' @param feature_repeats Number of times each structure map is repeated.
#' @param ... Additional arguments passed to \code{\link{segment_large}}.
#'
#' @return A \code{segment_large} result with structural metadata attached.
#' @export
segment_structures <- function(input,
                               scores,
                               Ncomp = 20,
                               mask_name = "structure_mask",
                               feature_maps = c(
                                 "combined_structure_score",
                                 "starlet_score",
                                 "ridge_score"
                               ),
                               feature_repeats = 6L,
                               ...) {
  if (is.null(scores$masks[[mask_name]])) {
    stop("`scores` does not contain mask `", mask_name, "`. Run threshold_structures() first.")
  }

  augmented <- .structure_augment_input(
    input,
    scores = scores,
    feature_maps = feature_maps,
    repeats = as.integer(feature_repeats)
  )

  out <- segment_large(
    augmented,
    Ncomp = Ncomp,
    mask = scores$masks[[mask_name]],
    valid_mode = "finite",
    feature_scale = "robust_col",
    spatial_weight = 0.25,
    ...
  )

  out$structure_info <- list(
    mask_name = mask_name,
    feature_maps = feature_maps,
    feature_repeats = feature_repeats,
    threshold = scores$threshold
  )
  out
}
