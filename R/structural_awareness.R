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
#' @importFrom graphics hist
#' @importFrom utils head tail
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

.structure_keep_touching <- function(mask, seed) {
  cc <- .structure_components(mask)
  if (!length(cc$sizes)) return(mask & FALSE)
  labs <- unique(cc$label[seed & mask])
  labs <- labs[labs > 0]
  if (!length(labs)) labs <- which.max(cc$sizes)
  out <- cc$label %in% labs
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

.structure_angle_diff_deg <- function(a, b) {
  abs(((a - b + 90) %% 180) - 90)
}

.structure_weighted_ellipse <- function(mask, weight, center = NULL) {
  pts <- which(mask, arr.ind = TRUE)
  if (nrow(pts) < 8L) return(NULL)

  wt <- weight[pts]
  wt[!is.finite(wt)] <- 0
  wt <- pmax(wt, 0)
  if (sum(wt) <= 0) wt[] <- 1

  xy <- cbind(x = pts[, 2], y = pts[, 1])
  if (is.null(center)) {
    centroid_xy <- colSums(xy * wt) / sum(wt)
  } else {
    centroid_xy <- c(x = unname(center[2]), y = unname(center[1]))
  }

  xc <- sweep(xy, 2, centroid_xy, "-")
  cov <- crossprod(xc, xc * wt) / sum(wt)
  eig <- eigen(cov)
  ord <- order(eig$values, decreasing = TRUE)
  val <- pmax(eig$values[ord], 0)
  vec <- eig$vectors[, ord, drop = FALSE]

  semi_major <- 2 * sqrt(val[1])
  semi_minor <- 2 * sqrt(val[2])
  axis_ratio <- semi_minor / max(semi_major, .Machine$double.eps)
  pa_deg <- (atan2(vec[1, 1], vec[2, 1]) * 180 / pi) %% 180

  list(
    x_centroid = centroid_xy[["x"]],
    y_centroid = centroid_xy[["y"]],
    semi_major_px = semi_major,
    semi_minor_px = semi_minor,
    axis_ratio = axis_ratio,
    ellipticity = 1 - axis_ratio,
    pa_deg = pa_deg
  )
}

.structure_elliptic_radius <- function(nr, nc, center, pa_deg, axis_ratio) {
  x <- col(matrix(0, nr, nc)) - center[2]
  y <- row(matrix(0, nr, nc)) - center[1]
  theta <- pa_deg * pi / 180
  major <- x * sin(theta) + y * cos(theta)
  minor <- x * cos(theta) - y * sin(theta)
  sqrt(major^2 + (minor / max(axis_ratio, 0.08))^2)
}

.structure_modified_ferrer <- function(r, sigma0, rout, alpha, beta) {
  power <- pmax(2 - beta, 0.15)
  z <- 1 - (pmax(r, 0) / max(rout, .Machine$double.eps))^power
  sigma0 * pmax(z, 0)^pmax(alpha, .Machine$double.eps)
}

.structure_fit_modified_ferrer <- function(radius,
                                           profile,
                                           area = NULL,
                                           radius_min = 1,
                                           radius_max = Inf) {
  ok <- is.finite(radius) & is.finite(profile)
  radius <- radius[ok]
  y <- profile[ok]
  if (!is.null(area)) {
    wt <- area[ok]
  } else {
    wt <- rep(1, length(y))
  }
  wt[!is.finite(wt) | wt <= 0] <- 1

  if (length(y) < 5L || max(y, na.rm = TRUE) <= 0) return(NULL)
  y <- pmax(y, 0)
  ymax <- max(y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) return(NULL)

  radius_max <- min(radius_max, max(radius, na.rm = TRUE) * 1.25, na.rm = TRUE)
  radius_min <- max(radius_min, min(radius[radius > 0], na.rm = TRUE), na.rm = TRUE)
  if (!is.finite(radius_max) || radius_max <= radius_min) return(NULL)

  wt <- sqrt(wt / stats::median(wt, na.rm = TRUE))
  wt[!is.finite(wt)] <- 1

  edge_guess <- max(radius[y >= 0.08 * ymax], na.rm = TRUE)
  if (!is.finite(edge_guess)) edge_guess <- stats::quantile(radius, 0.75, na.rm = TRUE)
  edge_guess <- min(max(edge_guess, radius_min), radius_max)

  starts <- expand.grid(
    sigma0 = c(ymax, 1.15 * ymax),
    rout = unique(pmin(pmax(c(edge_guess, 1.15 * edge_guess, 0.75 * radius_max), radius_min), radius_max)),
    alpha = c(0.5, 1.0, 2.0, 4.0),
    beta = c(0.0, 0.7, 1.2),
    KEEP.OUT.ATTRS = FALSE
  )

  objective <- function(par) {
    pred <- .structure_modified_ferrer(
      radius,
      sigma0 = par[1],
      rout = par[2],
      alpha = par[3],
      beta = par[4]
    )
    sum(wt * (y - pred)^2, na.rm = TRUE)
  }

  lower <- c(0, radius_min, 0.15, 0.0)
  upper <- c(2 * ymax, radius_max, 8.0, 1.85)
  best <- NULL
  best_value <- Inf

  for (i in seq_len(nrow(starts))) {
    start <- as.numeric(starts[i, ])
    start <- pmin(pmax(start, lower + 1e-6), upper - 1e-6)
    fit <- try(
      stats::optim(
        par = start,
        fn = objective,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = 300)
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) next
    if (is.finite(fit$value) && fit$value < best_value) {
      best <- fit
      best_value <- fit$value
    }
  }
  if (is.null(best)) return(NULL)

  pred <- .structure_modified_ferrer(
    radius,
    sigma0 = best$par[1],
    rout = best$par[2],
    alpha = best$par[3],
    beta = best$par[4]
  )
  rmse <- sqrt(mean((y - pred)^2, na.rm = TRUE))
  y_var <- stats::var(y, na.rm = TRUE)
  r2 <- if (is.finite(y_var) && y_var > 0) {
    1 - stats::var(y - pred, na.rm = TRUE) / y_var
  } else {
    NA_real_
  }

  list(
    sigma0 = best$par[1],
    rout = best$par[2],
    alpha = best$par[3],
    beta = best$par[4],
    radius = radius,
    fitted = pred,
    rmse = rmse,
    r2 = r2,
    convergence = best$convergence
  )
}

.structure_bar_empty <- function(scores,
                                 bar_score,
                                 reason,
                                 parameters,
                                 radial_profile = data.frame()) {
  empty <- matrix(FALSE, nrow(bar_score), ncol(bar_score))
  out <- list(
    bar_like = FALSE,
    confidence = 0,
    reason = reason,
    bar_mask = empty,
    candidate_mask = empty,
    bar_score = bar_score,
    diagnostics = data.frame(),
    radial_profile = radial_profile,
    scores = scores,
    parameters = parameters
  )
  class(out) <- c("capivara_bar_detection", "list")
  out
}

.structure_ring_empty <- function(scores,
                                  ring_score,
                                  reason,
                                  parameters,
                                  radial_profile = data.frame()) {
  empty <- matrix(FALSE, nrow(ring_score), ncol(ring_score))
  out <- list(
    ring_like = FALSE,
    confidence = 0,
    reason = reason,
    ring_mask = empty,
    candidate_mask = empty,
    ring_score = ring_score,
    ring_model = empty * NA_real_,
    diagnostics = data.frame(),
    radial_profile = radial_profile,
    scores = scores,
    parameters = parameters
  )
  class(out) <- c("capivara_ring_detection", "list")
  out
}

#' Detect a bar-like central structure
#'
#' @description
#' Experimental bar detector inspired by the modified Ferrer profile commonly
#' used for GALFIT-style bar and lens decompositions. The detector does not run
#' GALFIT; it estimates a central elongated support, fits a modified Ferrer
#' radial envelope to the bar-aligned structural evidence profile, and uses the
#' fitted outer truncation radius to define the bar support.
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array. Used
#'   only when \code{scores} is not supplied.
#' @param scores Optional object from \code{\link{score_structures}}.
#' @param support Optional support object from \code{\link{detect_support}}.
#' @param support_method Support detector used when \code{scores} is not
#'   supplied.
#' @param structure_scales Starlet scales used by \code{\link{score_structures}}.
#' @param max_radius_fraction Maximum tested radius as a fraction of the
#'   95-percent support radius.
#' @param threshold_quantile Quantile of the central bar score used to seed the
#'   candidate mask.
#' @param min_area Minimum candidate area in pixels.
#' @param min_ellipticity Minimum ellipticity required for a confident bar-like
#'   detection.
#' @param max_pa_scatter_deg Maximum allowed position-angle scatter across
#'   elliptical radial annuli.
#' @param min_confidence Minimum confidence for \code{bar_like = TRUE}.
#' @param min_model_fill Minimum fraction of the fitted elongated support that
#'   must be covered by the thresholded evidence mask. This prevents an ellipse
#'   from being accepted when it is supported only by a few disconnected pixels.
#' @param min_evidence_contrast Minimum median score contrast between the fitted
#'   elongated support and a local outer control region.
#' @param radius_method Method used to set the final bar-support radius.
#'   \code{"profile"} fits a modified Ferrer radial envelope along the bar
#'   ellipse and uses the fitted outer truncation radius. \code{"candidate"}
#'   uses the high-score candidate extent directly.
#' @param profile_floor_fraction Fraction of the fitted profile contrast above
#'   the local control level used as the outer bar-support floor when
#'   \code{radius_method = "profile"}.
#' @param support_axis_ratio_floor Minimum axis ratio used for the final filled
#'   support. The evidence fit may be very thin because it follows the brightest
#'   bar ridge; the segmentation support is allowed to be slightly wider.
#' @param support_radius_expand Multiplicative expansion applied to the final
#'   filled support radius. This is applied after the evidence tests, so it
#'   expands a justified support without making weak candidates easier to pass.
#' @param support_low_quantile Low bar-score quantile used to grow the final
#'   support by connected low-evidence pixels around the fitted support. This
#'   helps asymmetric or slightly curved bars whose outer ends are fainter than
#'   the central ridge.
#' @param use_isophote_profile If \code{TRUE}, grow the final support through a
#'   nested sequence of connected bar-score isophotes and keep the largest
#'   connected component whose ellipticity and position angle remain consistent
#'   with the high-evidence bar candidate. This makes the final support follow
#'   the observed bar body rather than a single static ellipse.
#' @param profile_steps Number of nested isophote levels used when
#'   \code{use_isophote_profile = TRUE}.
#' @param bar_body_quantile Quantile of the central-elongated score used to
#'   define a detached high-confidence bar-body support. This compact support is
#'   useful when the broader profile support starts to merge with spiral arms.
#' @param min_body_ellipticity Minimum ellipticity for the detached compact
#'   bar-body fallback. This is intentionally lower than
#'   \code{min_ellipticity}, because compact bar bodies can appear rounder than
#'   the full bar-support envelope.
#' @param solid_support If \code{TRUE}, the final \code{bar_mask} is the solid
#'   fitted elongated support inside the broad galaxy support. The thresholded
#'   score mask is still returned as \code{candidate_mask}. This avoids carving
#'   artificial holes through the centre of a bar.
#' @param keep_unaccepted If \code{TRUE}, keep the proposed support in
#'   \code{bar_mask} even when \code{bar_like = FALSE}. The default is
#'   \code{FALSE}: unaccepted supports are stored as \code{proposed_bar_mask},
#'   while \code{bar_mask} is empty.
#' @param starlet_args Extra arguments passed to \code{\link{score_structures}}
#'   when \code{scores} is not supplied.
#'
#' @return A \code{capivara_bar_detection} object containing \code{bar_mask},
#'   \code{candidate_mask}, \code{bar_score}, \code{diagnostics},
#'   \code{radial_profile}, and \code{ferrer_fit}. The logical field
#'   \code{bar_like} is deliberately conservative.
detect_bar <- function(input = NULL,
                       scores = NULL,
                       support = NULL,
                       support_method = "starlet",
                       structure_scales = 3:4,
                       max_radius_fraction = 0.75,
                       threshold_quantile = 0.58,
                       min_area = 20L,
                       min_ellipticity = 0.35,
                       max_pa_scatter_deg = 20,
                       min_confidence = 0.45,
                       min_model_fill = 0.20,
                       min_evidence_contrast = 0.04,
                       radius_method = c("profile", "candidate"),
                       profile_floor_fraction = 0.18,
                       support_axis_ratio_floor = 0.46,
                       support_radius_expand = 1.18,
                       support_low_quantile = 0.10,
                       use_isophote_profile = TRUE,
                       profile_steps = 10L,
                       bar_body_quantile = 0.65,
                       min_body_ellipticity = 0.10,
                       solid_support = TRUE,
                       keep_unaccepted = FALSE,
                       starlet_args = list()) {
  radius_method <- match.arg(radius_method)

  if (is.null(scores)) {
    if (is.null(input)) stop("Provide `input` or precomputed `scores`.")
    scores <- do.call(
      score_structures,
      list(
        input = input,
        support = support,
        support_method = support_method,
        structure_scales = structure_scales,
        starlet_args = starlet_args
      )
    )
  }

  combined <- scores$maps$combined_structure_score
  ridge <- scores$maps$ridge_score
  central <- scores$maps$central_elongated_score
  if (is.null(combined) || is.null(ridge) || is.null(central)) {
    stop("`scores` must contain combined, ridge, and central-elongated score maps.")
  }

  support_mask <- scores$support_mask & is.finite(combined)
  r_circ <- sqrt((row(combined) - scores$center[1])^2 + (col(combined) - scores$center[2])^2)
  r95 <- stats::quantile(r_circ[support_mask], 0.95, na.rm = TRUE)
  if (!is.finite(r95) || r95 <= 0) r95 <- max(r_circ[support_mask], na.rm = TRUE)

  central_domain <- support_mask & r_circ <= max_radius_fraction * r95
  bar_score <- .structure_robust_norm(0.55 * central + 0.25 * combined + 0.20 * ridge)
  bar_score[!central_domain] <- NA_real_

  parameters <- list(
    support_method = support_method,
    structure_scales = structure_scales,
    max_radius_fraction = max_radius_fraction,
    threshold_quantile = threshold_quantile,
    min_area = min_area,
    min_ellipticity = min_ellipticity,
    max_pa_scatter_deg = max_pa_scatter_deg,
    min_confidence = min_confidence,
    min_model_fill = min_model_fill,
    min_evidence_contrast = min_evidence_contrast,
    radius_method = radius_method,
    profile_floor_fraction = profile_floor_fraction,
    support_axis_ratio_floor = support_axis_ratio_floor,
    support_radius_expand = support_radius_expand,
    support_low_quantile = support_low_quantile,
    use_isophote_profile = use_isophote_profile,
    profile_steps = profile_steps,
    bar_body_quantile = bar_body_quantile,
    min_body_ellipticity = min_body_ellipticity,
    solid_support = solid_support,
    keep_unaccepted = keep_unaccepted
  )

  vals <- bar_score[central_domain & is.finite(bar_score)]
  if (length(vals) < min_area) {
    return(.structure_bar_empty(scores, bar_score, "too_few_central_pixels", parameters))
  }

  body_score <- .structure_robust_norm(0.70 * central + 0.30 * bar_score)
  body_score[!central_domain] <- NA_real_
  body_vals <- body_score[central_domain & is.finite(body_score)]
  body_cut <- stats::quantile(body_vals, bar_body_quantile, na.rm = TRUE)
  bar_body_mask <- central_domain & is.finite(body_score) & body_score >= body_cut
  bar_body_mask <- .structure_largest_center_component(
    bar_body_mask,
    scores$center,
    radius = max(4L, floor(0.08 * r95))
  )
  bar_body_mask <- .structure_drop_small(
    bar_body_mask & .structure_neighbor_count(bar_body_mask) >= 2L,
    min_area
  )
  body_fit <- .structure_weighted_ellipse(bar_body_mask, body_score, center = scores$center)
  body_area <- sum(bar_body_mask, na.rm = TRUE)
  body_like <- !is.null(body_fit) &&
    body_area >= min_area &&
    body_fit$ellipticity >= min_body_ellipticity

  cut <- stats::quantile(vals, threshold_quantile, na.rm = TRUE)
  candidate <- central_domain & bar_score >= cut
  candidate <- .structure_drop_small(
    candidate & .structure_neighbor_count(candidate) >= 2L,
    min_area
  )
  if (sum(candidate, na.rm = TRUE) < min_area) {
    return(.structure_bar_empty(scores, bar_score, "too_few_candidate_pixels", parameters))
  }

  ellipse <- .structure_weighted_ellipse(candidate, bar_score, center = scores$center)
  if (is.null(ellipse)) {
    return(.structure_bar_empty(scores, bar_score, "ellipse_fit_failed", parameters))
  }

  r_ell <- .structure_elliptic_radius(
    nrow(combined),
    ncol(combined),
    center = scores$center,
    pa_deg = ellipse$pa_deg,
    axis_ratio = ellipse$axis_ratio
  )
  candidate_radius <- stats::quantile(r_ell[candidate], 0.95, na.rm = TRUE)
  if (!is.finite(candidate_radius) || candidate_radius <= 0) {
    return(.structure_bar_empty(scores, bar_score, "invalid_bar_radius", parameters))
  }

  profile_radius_max <- max(max_radius_fraction * r95, candidate_radius * 1.35, 1)
  breaks <- seq(0, profile_radius_max, length.out = 13)
  radial_profile <- data.frame(
    r_inner = head(breaks, -1),
    r_outer = tail(breaks, -1),
    r_mid = 0.5 * (head(breaks, -1) + tail(breaks, -1)),
    median_score = NA_real_,
    q75_score = NA_real_,
    area_px = NA_integer_,
    pa_deg = NA_real_
  )
  pa_values <- numeric()
  for (i in seq_len(nrow(radial_profile))) {
    ann <- central_domain & r_ell >= radial_profile$r_inner[i] & r_ell < radial_profile$r_outer[i]
    radial_profile$area_px[i] <- sum(ann, na.rm = TRUE)
    radial_profile$median_score[i] <- stats::median(bar_score[ann], na.rm = TRUE)
    radial_profile$q75_score[i] <- stats::quantile(bar_score[ann], 0.75, na.rm = TRUE)
    ann_candidate <- ann & bar_score >= cut
    if (sum(ann_candidate, na.rm = TRUE) >= 8L) {
      ann_ellipse <- .structure_weighted_ellipse(ann_candidate, bar_score, center = scores$center)
      if (!is.null(ann_ellipse)) {
        radial_profile$pa_deg[i] <- ann_ellipse$pa_deg
        pa_values <- c(pa_values, ann_ellipse$pa_deg)
      }
    }
  }

  pa_scatter <- if (length(pa_values) >= 2L) {
    stats::median(.structure_angle_diff_deg(pa_values, ellipse$pa_deg), na.rm = TRUE)
  } else {
    max_pa_scatter_deg
  }

  peak_score <- max(radial_profile$median_score, na.rm = TRUE)
  outer_score <- stats::median(tail(radial_profile$median_score, 2), na.rm = TRUE)
  edge_to_peak <- outer_score / max(peak_score, .Machine$double.eps)
  truncation_score <- pmin(pmax(1 - edge_to_peak, 0), 1)
  elongation_score <- pmin(pmax((ellipse$ellipticity - min_ellipticity) /
    max(0.7 - min_ellipticity, 1e-6), 0), 1)
  pa_score <- pmin(pmax(1 - pa_scatter / max(max_pa_scatter_deg, 1e-6), 0), 1)
  support_score <- mean(bar_score[candidate], na.rm = TRUE)

  radius_cap <- max(max_radius_fraction * r95, candidate_radius, na.rm = TRUE)
  profile_score <- 0.6 * radial_profile$median_score + 0.4 * radial_profile$q75_score
  profile_score[!is.finite(profile_score) | radial_profile$area_px < 4L] <- NA_real_
  control_level <- stats::median(tail(profile_score, 3), na.rm = TRUE)
  if (!is.finite(control_level)) {
    control_level <- stats::median(profile_score, na.rm = TRUE)
  }
  profile_peak <- max(profile_score, na.rm = TRUE)
  if (!is.finite(profile_peak)) profile_peak <- peak_score
  profile_floor <- control_level +
    profile_floor_fraction * max(profile_peak - control_level, 0)

  ferrer_fit <- NULL
  if (identical(radius_method, "profile")) {
    ferrer_signal <- pmax(profile_score - control_level, 0)
    ferrer_fit <- .structure_fit_modified_ferrer(
      radius = radial_profile$r_mid,
      profile = ferrer_signal,
      area = radial_profile$area_px,
      radius_min = max(1, 0.55 * candidate_radius),
      radius_max = radius_cap
    )
    if (!is.null(ferrer_fit) && is.finite(ferrer_fit$rout)) {
      bar_radius <- max(candidate_radius, ferrer_fit$rout, na.rm = TRUE)
    } else {
      supported_bins <- which(
        is.finite(profile_score) &
          profile_score >= profile_floor &
          radial_profile$r_mid <= max_radius_fraction * r95
      )
      supported_bins <- supported_bins[
        supported_bins <= min(length(profile_score), which.max(profile_score) + 5L)
      ]
      if (length(supported_bins)) {
        bar_radius <- max(radial_profile$r_outer[supported_bins], candidate_radius, na.rm = TRUE)
      } else {
        bar_radius <- candidate_radius
      }
    }
  } else {
    bar_radius <- candidate_radius
  }
  bar_radius <- min(bar_radius, radius_cap, na.rm = TRUE)
  radial_profile$ferrer_signal <- pmax(profile_score - control_level, 0)
  radial_profile$ferrer_model <- NA_real_
  if (!is.null(ferrer_fit)) {
    radial_profile$ferrer_model <- control_level + .structure_modified_ferrer(
      radial_profile$r_mid,
      sigma0 = ferrer_fit$sigma0,
      rout = ferrer_fit$rout,
      alpha = ferrer_fit$alpha,
      beta = ferrer_fit$beta
    )
  }
  ferrer_supported <- !is.null(ferrer_fit) &&
    is.finite(ferrer_fit$rout) &&
    ferrer_fit$rout >= 0.8 * candidate_radius &&
    (
      !is.finite(ferrer_fit$r2) ||
        ferrer_fit$r2 >= -0.20
    )

  support_r_ell <- .structure_elliptic_radius(
    nrow(combined),
    ncol(combined),
    center = scores$center,
    pa_deg = ellipse$pa_deg,
    axis_ratio = max(ellipse$axis_ratio, support_axis_ratio_floor)
  )
  model_mask <- central_domain & r_ell <= bar_radius
  support_bar_radius <- min(
    bar_radius * support_radius_expand,
    max(radius_cap, bar_radius),
    na.rm = TRUE
  )
  support_model_mask <- central_domain & support_r_ell <= support_bar_radius
  bar_envelope <- support_model_mask
  low_cut <- stats::quantile(vals, support_low_quantile, na.rm = TRUE)
  low_support <- central_domain & is.finite(bar_score) & bar_score >= low_cut
  low_support <- .structure_drop_small(
    low_support & .structure_neighbor_count(low_support) >= 2L,
    min_area
  )
  low_cc <- .structure_components(low_support)
  seed_labs <- unique(low_cc$label[support_model_mask | candidate])
  seed_labs <- seed_labs[seed_labs > 0]
  connected_low_support <- if (length(seed_labs)) {
    out <- low_cc$label %in% seed_labs
    dim(out) <- dim(low_support)
    out
  } else {
    support_model_mask
  }

  isophote_profile <- data.frame()
  profile_support_mask <- support_model_mask
  if (isTRUE(use_isophote_profile)) {
    qs <- seq(
      from = threshold_quantile,
      to = min(threshold_quantile, support_low_quantile),
      length.out = max(2L, as.integer(profile_steps))
    )
    profile_rows <- list()
    best_profile <- NULL
    best_area <- 0L

    for (ii in seq_along(qs)) {
      q <- qs[ii]
      q_cut <- stats::quantile(vals, q, na.rm = TRUE)
      level_mask <- central_domain & is.finite(bar_score) & bar_score >= q_cut
      level_mask <- .structure_drop_small(
        level_mask & .structure_neighbor_count(level_mask) >= 2L,
        min_area
      )
      cc <- .structure_components(level_mask)
      labs <- unique(cc$label[candidate])
      labs <- labs[labs > 0]
      if (!length(labs)) next

      comp <- cc$label %in% labs
      dim(comp) <- dim(level_mask)
      comp <- comp & bar_envelope
      comp <- .structure_drop_small(
        comp & .structure_neighbor_count(comp) >= 2L,
        min_area
      )
      if (sum(comp, na.rm = TRUE) < min_area) next

      comp_fit <- .structure_weighted_ellipse(comp, bar_score, center = scores$center)
      if (is.null(comp_fit)) next
      pa_diff <- .structure_angle_diff_deg(comp_fit$pa_deg, ellipse$pa_deg)
      comp_radius <- stats::quantile(r_ell[comp], 0.95, na.rm = TRUE)
      comp_score <- stats::median(bar_score[comp], na.rm = TRUE)
      profile_rows[[length(profile_rows) + 1L]] <- data.frame(
        level = ii,
        quantile = q,
        threshold = as.numeric(q_cut),
        area_px = sum(comp, na.rm = TRUE),
        radius_px = comp_radius,
        ellipticity = comp_fit$ellipticity,
        axis_ratio = comp_fit$axis_ratio,
        pa_deg = comp_fit$pa_deg,
        pa_diff_deg = pa_diff,
        median_score = comp_score,
        accepted = FALSE,
        stringsAsFactors = FALSE
      )

      accepted <- isTRUE(
        comp_fit$ellipticity >= max(0.12, min_ellipticity * 0.55) &&
          pa_diff <= max(25, max_pa_scatter_deg * 1.5) &&
          is.finite(comp_radius) &&
          comp_radius <= radius_cap * 1.25
      )
      if (accepted && sum(comp, na.rm = TRUE) > best_area) {
        best_area <- sum(comp, na.rm = TRUE)
        best_profile <- comp
        profile_rows[[length(profile_rows)]]$accepted <- TRUE
      }
    }

    if (length(profile_rows)) {
      isophote_profile <- do.call(rbind, profile_rows)
    }
    if (!is.null(best_profile)) {
      profile_support_mask <- best_profile
    }
  }
  local_control <- central_domain &
    r_ell > bar_radius &
    r_ell <= min(max_radius_fraction * r95, 1.45 * bar_radius)
  if (sum(local_control, na.rm = TRUE) < min_area) {
    local_control <- central_domain & !model_mask
  }
  model_fill <- sum(candidate & model_mask, na.rm = TRUE) /
    max(sum(model_mask, na.rm = TRUE), 1)
  model_score <- stats::median(bar_score[model_mask], na.rm = TRUE)
  control_score <- stats::median(bar_score[local_control], na.rm = TRUE)
  if (!is.finite(control_score)) control_score <- stats::median(bar_score[central_domain], na.rm = TRUE)
  evidence_contrast <- model_score - control_score
  fill_score <- pmin(pmax(model_fill / max(min_model_fill, 1e-6), 0), 1)
  contrast_score <- pmin(pmax(evidence_contrast / max(min_evidence_contrast, 1e-6), 0), 1)
  model_supported <- isTRUE(
    model_fill >= min_model_fill &&
      evidence_contrast >= min_evidence_contrast
  )
  confidence <- mean(
    c(elongation_score, pa_score, truncation_score, support_score, fill_score, contrast_score),
    na.rm = TRUE
  )

  profile_like <- isTRUE(
    confidence >= min_confidence &&
      ellipse$ellipticity >= min_ellipticity &&
      pa_scatter <= max_pa_scatter_deg &&
      model_supported &&
      (identical(radius_method, "candidate") || ferrer_supported)
  )
  bar_like <- isTRUE(profile_like || body_like)

  if (isTRUE(solid_support)) {
    proposed_bar_mask <- support_model_mask |
      ((connected_low_support | profile_support_mask) & bar_envelope)
  } else {
    proposed_bar_mask <- (support_model_mask |
      ((connected_low_support | profile_support_mask) & bar_envelope)) &
      bar_score >= stats::quantile(bar_score[candidate], 0.25, na.rm = TRUE)
  }
  proposed_bar_mask <- .structure_drop_small(
    proposed_bar_mask & .structure_neighbor_count(proposed_bar_mask) >= 2L,
    min_area
  )
  proposed_bar_mask <- .structure_largest_center_component(
    proposed_bar_mask,
    scores$center,
    radius = max(4L, floor(0.20 * bar_radius))
  )
  proposed_bar_mask <- .structure_drop_small(
    proposed_bar_mask & .structure_neighbor_count(proposed_bar_mask) >= 2L,
    min_area
  )
  bar_mask <- if (isTRUE(profile_like)) {
    proposed_bar_mask
  } else if (isTRUE(body_like)) {
    bar_body_mask
  } else if (isTRUE(keep_unaccepted)) {
    proposed_bar_mask
  } else {
    proposed_bar_mask & FALSE
  }

  diagnostics <- data.frame(
    bar_like = bar_like,
    profile_like = profile_like,
    body_like = body_like,
    confidence = confidence,
    ellipticity = ellipse$ellipticity,
    axis_ratio = ellipse$axis_ratio,
    support_axis_ratio = max(ellipse$axis_ratio, support_axis_ratio_floor),
    pa_deg = ellipse$pa_deg,
    pa_scatter_deg = pa_scatter,
    bar_radius_px = bar_radius,
    support_bar_radius_px = support_bar_radius,
    candidate_radius_px = candidate_radius,
    bar_body_area_px = body_area,
    bar_body_ellipticity = if (!is.null(body_fit)) body_fit$ellipticity else NA_real_,
    bar_body_axis_ratio = if (!is.null(body_fit)) body_fit$axis_ratio else NA_real_,
    bar_body_pa_deg = if (!is.null(body_fit)) body_fit$pa_deg else NA_real_,
    bar_body_threshold = as.numeric(body_cut),
    support_low_threshold = as.numeric(low_cut),
    isophote_profile_area_px = sum(profile_support_mask, na.rm = TRUE),
    ferrer_supported = ferrer_supported,
    ferrer_sigma0 = if (!is.null(ferrer_fit)) ferrer_fit$sigma0 else NA_real_,
    ferrer_rout_px = if (!is.null(ferrer_fit)) ferrer_fit$rout else NA_real_,
    ferrer_alpha = if (!is.null(ferrer_fit)) ferrer_fit$alpha else NA_real_,
    ferrer_beta = if (!is.null(ferrer_fit)) ferrer_fit$beta else NA_real_,
    ferrer_r2 = if (!is.null(ferrer_fit)) ferrer_fit$r2 else NA_real_,
    ferrer_rmse = if (!is.null(ferrer_fit)) ferrer_fit$rmse else NA_real_,
    profile_peak = profile_peak,
    profile_control = control_level,
    profile_floor = profile_floor,
    truncation_score = truncation_score,
    support_score = support_score,
    model_fill = model_fill,
    model_score = model_score,
    control_score = control_score,
    evidence_contrast = evidence_contrast,
    model_supported = model_supported,
    area_px = sum(bar_mask, na.rm = TRUE),
    threshold = as.numeric(cut),
    stringsAsFactors = FALSE
  )

  out <- list(
    bar_like = bar_like,
    confidence = confidence,
    reason = if (bar_like) "central_elongated_support" else "low_confidence_central_elongated_support",
    bar_mask = bar_mask,
    proposed_bar_mask = proposed_bar_mask,
    bar_body_mask = bar_body_mask,
    body_score = body_score,
    candidate_mask = candidate,
    bar_score = bar_score,
    diagnostics = diagnostics,
    radial_profile = radial_profile,
    isophote_profile = isophote_profile,
    ferrer_fit = ferrer_fit,
    scores = scores,
    parameters = parameters
  )
  class(out) <- c("capivara_bar_detection", "list")
  out
}

#' Detect a ring-like annular structure
#'
#' @description
#' Experimental ring detector for structural-aware segmentation. The detector
#' combines the neutral Capivara structure score with a simple elliptical
#' annulus model: it estimates the galaxy geometry from the supported light,
#' builds an elliptical radial profile, finds the strongest non-central
#' annular peak, and thresholds only pixels close to that peak. The central
#' region is explicitly excluded, so a bright bulge is not labelled as a ring.
#'
#' @param input A FITS-like object with \code{imDat}, or a 3-D cube array. Used
#'   only when \code{scores} is not supplied.
#' @param scores Optional object from \code{\link{score_structures}}.
#' @param support Optional support object from \code{\link{detect_support}}.
#' @param support_method Support detector used when \code{scores} is not
#'   supplied.
#' @param structure_scales Starlet scales used by \code{\link{score_structures}}.
#' @param footprint_quantile Quantile of the collapsed positive image used to
#'   build the broad galaxy footprint before searching for annular structure.
#' @param bulge_fraction Elliptical-radius fraction used to exclude the compact
#'   central component from the ring support.
#' @param min_inner_radius_fraction Minimum elliptical radius considered for
#'   the ring, as a fraction of the 95-percent support radius. The default is
#'   intentionally outside the compact bulge.
#' @param max_radius_fraction Maximum elliptical radius considered for the
#'   ring, as a fraction of the 95-percent support radius.
#' @param threshold_quantile Quantile of the ring score used to seed the
#'   candidate mask.
#' @param min_area Minimum candidate area in pixels.
#' @param min_azimuth_coverage Minimum azimuthal coverage, in fractions of
#'   \eqn{2\pi}, required for \code{ring_like = TRUE}.
#' @param min_confidence Minimum confidence for \code{ring_like = TRUE}.
#' @param keep_largest_component If \code{TRUE}, keep only the largest connected
#'   ring component after thresholding. This is useful for suppressing small
#'   starlet islands around the main annular feature.
#' @param ring_model_floor Minimum value of the fitted annular model kept in
#'   the final ring mask. Lower values keep faint/asymmetric ring arcs.
#' @param n_radial_bins Number of bins used in the elliptical radial profile.
#' @param starlet_args Extra arguments passed to \code{\link{score_structures}}
#'   when \code{scores} is not supplied.
#'
#' @return A \code{capivara_ring_detection} object containing \code{ring_mask},
#'   \code{candidate_mask}, \code{ring_score}, \code{diagnostics}, and
#'   \code{radial_profile}. The logical field \code{ring_like} is deliberately
#'   conservative.
detect_ring <- function(input = NULL,
                        scores = NULL,
                        support = NULL,
                        support_method = "starlet",
                        structure_scales = 3:4,
                        footprint_quantile = 0.20,
                        bulge_fraction = 0.25,
                        min_inner_radius_fraction = 0.22,
                        max_radius_fraction = 1.05,
                        threshold_quantile = 0.05,
                        min_area = 24L,
                        min_azimuth_coverage = 0.22,
                        min_confidence = 0.35,
                        keep_largest_component = TRUE,
                        ring_model_floor = 0.08,
                        n_radial_bins = 18L,
                        starlet_args = list()) {
  if (is.null(scores)) {
    if (is.null(input)) stop("Provide `input` or precomputed `scores`.")
    scores <- do.call(
      score_structures,
      list(
        input = input,
        support = support,
        support_method = support_method,
        structure_scales = structure_scales,
        starlet_args = starlet_args
      )
    )
  }

  combined <- scores$maps$combined_structure_score
  collapsed <- scores$maps$collapsed
  if (is.null(combined) || is.null(collapsed)) {
    stop("`scores` must contain combined structure and collapsed maps.")
  }

  pos <- collapsed[is.finite(collapsed) & collapsed > 0]
  if (length(pos) >= 16L) {
    raw <- collapsed > stats::quantile(pos, footprint_quantile, na.rm = TRUE)
    support_mask <- .structure_largest_center_component(raw, scores$center, radius = 5L)
    support_mask <- .structure_dilate(support_mask, 1L) &
      collapsed > stats::quantile(
        pos,
        max(0.05, footprint_quantile * 0.6),
        na.rm = TRUE
      )
    support_mask <- .structure_drop_small(
      support_mask & .structure_neighbor_count(support_mask) >= 3L,
      min_area
    )
  } else {
    support_mask <- if (!is.null(scores$support$mask)) {
      scores$support$mask
    } else {
      scores$support_mask
    }
  }
  support_mask <- support_mask & is.finite(combined)
  parameters <- list(
    support_method = support_method,
    structure_scales = structure_scales,
    footprint_quantile = footprint_quantile,
    bulge_fraction = bulge_fraction,
    min_inner_radius_fraction = min_inner_radius_fraction,
    max_radius_fraction = max_radius_fraction,
    threshold_quantile = threshold_quantile,
    min_area = min_area,
    min_azimuth_coverage = min_azimuth_coverage,
    min_confidence = min_confidence,
    keep_largest_component = keep_largest_component,
    ring_model_floor = ring_model_floor,
    n_radial_bins = n_radial_bins
  )

  if (sum(support_mask, na.rm = TRUE) < min_area) {
    ring_score <- matrix(NA_real_, nrow(combined), ncol(combined))
    return(.structure_ring_empty(scores, ring_score, "too_few_supported_pixels", parameters))
  }

  idx <- which(support_mask & is.finite(collapsed), arr.ind = TRUE)
  if (nrow(idx) < 8L) {
    ring_score <- matrix(NA_real_, nrow(combined), ncol(combined))
    return(.structure_ring_empty(scores, ring_score, "ellipse_fit_failed", parameters))
  }
  ellipse_weight <- pmax(collapsed[idx], 0)
  if (!any(ellipse_weight > 0, na.rm = TRUE)) ellipse_weight[] <- 1
  xy <- cbind(idx[, 2] - scores$center[2], idx[, 1] - scores$center[1])
  cov <- stats::cov.wt(xy, wt = ellipse_weight)$cov
  eig <- eigen(cov)
  theta <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])
  axis_ratio <- sqrt(max(min(eig$values) / max(eig$values), 0.15))
  xx <- col(combined) - scores$center[2]
  yy <- row(combined) - scores$center[1]
  major <- xx * cos(theta) + yy * sin(theta)
  minor <- -xx * sin(theta) + yy * cos(theta)
  r_ell <- sqrt(major^2 + (minor / axis_ratio)^2)
  pa_deg <- (theta * 180 / pi) %% 180
  galaxy_ellipse <- list(
    axis_ratio = axis_ratio,
    ellipticity = 1 - axis_ratio,
    pa_deg = pa_deg
  )
  r95 <- stats::quantile(r_ell[support_mask], 0.95, na.rm = TRUE)
  if (!is.finite(r95) || r95 <= 0) r95 <- max(r_ell[support_mask], na.rm = TRUE)

  r_min <- max(min_inner_radius_fraction, bulge_fraction) * r95
  r_max <- max_radius_fraction * r95
  domain <- support_mask & r_ell >= r_min & r_ell <= r_max
  if (sum(domain, na.rm = TRUE) < min_area) {
    ring_score <- matrix(NA_real_, nrow(combined), ncol(combined))
    return(.structure_ring_empty(scores, ring_score, "too_few_annular_pixels", parameters))
  }

  ring_evidence <- combined
  if (!is.null(scores$support$decomposition$w) &&
      length(scores$support$decomposition$w) >= 4L) {
    w <- scores$support$decomposition$w
    ring_evidence <- .structure_robust_norm(pmax(w[[2]], 0)) +
      .structure_robust_norm(pmax(w[[3]], 0)) +
      0.5 * .structure_robust_norm(pmax(w[[4]], 0))
    ring_evidence <- .structure_robust_norm(ring_evidence)
  }
  ring_evidence[!support_mask] <- NA_real_

  breaks <- seq(r_min, r_max, length.out = max(6L, as.integer(n_radial_bins)) + 1L)
  radial_profile <- data.frame(
    r_inner = head(breaks, -1),
    r_outer = tail(breaks, -1),
    r_mid = 0.5 * (head(breaks, -1) + tail(breaks, -1)),
    median_score = NA_real_,
    q75_score = NA_real_,
    baseline_score = NA_real_,
    bump_score = NA_real_,
    area_px = NA_integer_
  )
  for (i in seq_len(nrow(radial_profile))) {
    ann <- domain & r_ell >= radial_profile$r_inner[i] & r_ell < radial_profile$r_outer[i]
    radial_profile$area_px[i] <- sum(ann, na.rm = TRUE)
    radial_profile$median_score[i] <- stats::median(ring_evidence[ann], na.rm = TRUE)
    radial_profile$q75_score[i] <- stats::quantile(ring_evidence[ann], 0.75, na.rm = TRUE)
  }

  profile_score <- 0.55 * radial_profile$median_score + 0.45 * radial_profile$q75_score
  profile_score[!is.finite(profile_score) | radial_profile$area_px < 4L] <- NA_real_
  if (all(!is.finite(profile_score))) {
    ring_score <- matrix(NA_real_, nrow(combined), ncol(combined))
    return(.structure_ring_empty(scores, ring_score, "radial_profile_failed", parameters, radial_profile))
  }

  baseline <- rep(NA_real_, length(profile_score))
  for (i in seq_along(profile_score)) {
    local <- seq(max(1L, i - 2L), min(length(profile_score), i + 2L))
    far <- setdiff(seq_along(profile_score), local)
    if (length(far) >= 3L) {
      baseline[i] <- stats::median(profile_score[far], na.rm = TRUE)
    }
  }
  if (all(!is.finite(baseline))) {
    baseline[] <- stats::median(profile_score, na.rm = TRUE)
  }
  baseline[!is.finite(baseline)] <- stats::median(profile_score, na.rm = TRUE)
  bump_score <- profile_score - baseline
  radial_profile$baseline_score <- baseline
  radial_profile$bump_score <- bump_score

  if (any(is.finite(bump_score) & bump_score > 0, na.rm = TRUE)) {
    peak_i <- which.max(bump_score)
  } else {
    peak_i <- which.max(profile_score)
  }
  ring_radius <- radial_profile$r_mid[peak_i]
  peak_bump <- bump_score[peak_i]
  if (!is.finite(peak_bump) || peak_bump <= 0) peak_bump <- profile_score[peak_i]
  ok_bins <- which(is.finite(bump_score) & bump_score >= 0.25 * peak_bump)
  ok_bins <- ok_bins[
    ok_bins >= max(1L, peak_i - 4L) &
      ok_bins <= min(nrow(radial_profile), peak_i + 4L)
  ]
  if (length(ok_bins)) {
    ring_width <- 0.5 * diff(range(c(
      radial_profile$r_inner[ok_bins],
      radial_profile$r_outer[ok_bins]
    )))
  } else {
    ring_width <- 0.13 * r95
  }
  ring_width <- max(ring_width, 0.10 * r95, 1)

  ring_score <- .structure_robust_norm(ring_evidence)
  ring_score[!support_mask] <- NA_real_

  vals <- ring_score[is.finite(ring_score) & ring_score > 0]
  if (length(vals) < min_area) {
    return(.structure_ring_empty(scores, ring_score, "too_few_ring_score_pixels", parameters, radial_profile))
  }

  seed_cut <- stats::quantile(vals, threshold_quantile, na.rm = TRUE)
  ring_seed <- is.finite(ring_score) & ring_score >= seed_cut
  base_support <- if (!is.null(scores$support$mask)) {
    scores$support$mask
  } else {
    scores$support_mask
  }
  structure_support <- support_mask & (base_support | ring_seed)
  structure_support <- .structure_drop_small(
    structure_support & .structure_neighbor_count(structure_support) >= 3L,
    min_area
  )
  bulge_mask <- structure_support & r_ell <= bulge_fraction * r95
  bulge_mask <- .structure_dilate(bulge_mask, 1L) & structure_support
  candidate <- structure_support &
    !bulge_mask &
    r_ell <= max_radius_fraction * r95 &
    ring_seed
  candidate <- .structure_drop_small(
    candidate & .structure_neighbor_count(candidate) >= 3L,
    min_area
  )
  if (keep_largest_component && sum(candidate, na.rm = TRUE) >= min_area) {
    cc <- .structure_components(candidate)
    keep <- which.max(cc$sizes)
    candidate <- cc$label == keep
    dim(candidate) <- dim(ring_score)
  }
  model_radius <- stats::median(r_ell[candidate], na.rm = TRUE)
  model_width <- stats::mad(r_ell[candidate], constant = 1.4826, na.rm = TRUE)
  if (!is.finite(model_radius)) model_radius <- ring_radius
  if (!is.finite(model_width) || model_width <= 0) model_width <- ring_width
  model_width <- max(1.5 * model_width, 0.18 * r95, 1)
  ring_model <- exp(-0.5 * ((r_ell - model_radius) / model_width)^2)
  ring_model[!support_mask] <- NA_real_
  candidate <- candidate & ring_model >= ring_model_floor
  candidate <- .structure_drop_small(
    candidate & .structure_neighbor_count(candidate) >= 3L,
    min_area
  )
  if (sum(candidate, na.rm = TRUE) < min_area) {
    return(.structure_ring_empty(scores, ring_score, "too_few_model_supported_pixels", parameters, radial_profile))
  }
  if (sum(candidate, na.rm = TRUE) < min_area) {
    return(.structure_ring_empty(scores, ring_score, "too_few_candidate_pixels", parameters, radial_profile))
  }

  theta <- atan2(row(candidate) - scores$center[1], col(candidate) - scores$center[2])
  theta_bin <- floor(((theta + pi) / (2 * pi)) * 32)
  theta_bin <- pmin(pmax(theta_bin, 0L), 31L)
  azimuth_coverage <- length(unique(theta_bin[candidate])) / 32

  radial_concentration <- 1 - stats::mad(r_ell[candidate], constant = 1.4826, na.rm = TRUE) /
    max(stats::median(r_ell[candidate], na.rm = TRUE), 1)
  radial_concentration <- pmin(pmax(radial_concentration, 0), 1)
  support_score <- mean(ring_score[candidate], na.rm = TRUE)
  central_leakage <- mean(ring_score[support_mask & r_ell < r_min], na.rm = TRUE)
  if (!is.finite(central_leakage)) central_leakage <- 0
  leakage_score <- pmin(pmax(1 - central_leakage / max(support_score, 1e-6), 0), 1)
  coverage_score <- pmin(pmax(azimuth_coverage / max(min_azimuth_coverage, 1e-6), 0), 1)
  confidence <- mean(
    c(radial_concentration, support_score, leakage_score, coverage_score),
    na.rm = TRUE
  )

  ring_like <- isTRUE(
    confidence >= min_confidence &&
      azimuth_coverage >= min_azimuth_coverage
  )

  diagnostics <- data.frame(
    ring_like = ring_like,
    confidence = confidence,
    ring_radius_px = ring_radius,
    ring_width_px = ring_width,
    model_radius_px = model_radius,
    model_width_px = model_width,
    axis_ratio = galaxy_ellipse$axis_ratio,
    ellipticity = galaxy_ellipse$ellipticity,
    pa_deg = galaxy_ellipse$pa_deg,
    azimuth_coverage = azimuth_coverage,
    radial_concentration = radial_concentration,
    support_score = support_score,
    central_leakage = central_leakage,
    area_px = sum(candidate, na.rm = TRUE),
    threshold = as.numeric(seed_cut),
    stringsAsFactors = FALSE
  )

  out <- list(
    ring_like = ring_like,
    confidence = confidence,
    reason = if (ring_like) "ring_like_structure" else "low_confidence_ring_like_structure",
    ring_mask = candidate,
    candidate_mask = candidate,
    ring_score = ring_score,
    ring_model = ring_model,
    diagnostics = diagnostics,
    radial_profile = radial_profile,
    scores = scores,
    parameters = parameters
  )
  class(out) <- c("capivara_ring_detection", "list")
  out
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
