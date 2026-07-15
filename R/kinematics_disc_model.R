#' Axisymmetric arctangent rotating-disc velocity model
#'
#' @param par Named vector/list with `vsys`, `vmax`, and `Rt`.
#' @param R Deprojected radius in spaxels.
#' @param theta Deprojected azimuth.
#' @param inc_rad Inclination in radians.
#' @return Model line-of-sight velocity.
#' @noRd
disc_velocity_model <- function(par, R, theta, inc_rad) {
  vsys <- unname(par[["vsys"]])
  vmax <- unname(par[["vmax"]])
  Rt <- unname(par[["Rt"]])
  if (!is.finite(Rt) || Rt <= 0) {
    return(rep(NA_real_, length(R)))
  }
  vrot <- vmax * (2 / pi) * atan(R / Rt)
  vsys + sin(inc_rad) * vrot * cos(theta)
}

.capivara_piecewise_basis <- function(R, knots) {
  n <- length(knots)
  basis <- matrix(0, nrow = length(R), ncol = n)
  if (n < 2L) {
    basis[, 1] <- 1
    return(basis)
  }
  for (i in seq_along(R)) {
    r <- R[i]
    if (!is.finite(r)) {
      basis[i, ] <- NA_real_
    } else if (r <= knots[1]) {
      basis[i, 1] <- 1
    } else if (r >= knots[n]) {
      basis[i, n] <- 1
    } else {
      hi <- which(knots >= r)[1]
      lo <- hi - 1L
      frac <- (r - knots[lo]) / (knots[hi] - knots[lo])
      basis[i, lo] <- 1 - frac
      basis[i, hi] <- frac
    }
  }
  basis
}

.capivara_second_derivative_penalty <- function(n, lambda) {
  if (!is.finite(lambda) || lambda <= 0 || n < 3L) {
    return(matrix(numeric(0), 0, n))
  }
  out <- matrix(0, nrow = n - 2L, ncol = n)
  for (i in seq_len(n - 2L)) {
    out[i, i:(i + 2L)] <- c(1, -2, 1)
  }
  sqrt(lambda) * out
}

.capivara_fit_linear_velocity_model <- function(y, design, weights,
                                                smooth_groups,
                                                ridge = NULL,
                                                robust = FALSE,
                                                robust_k = 1.345,
                                                robust_maxit = 20L,
                                                robust_tol = 1e-3) {
  ok <- is.finite(y) & rowSums(!is.finite(design)) == 0
  y <- y[ok]
  design <- design[ok, , drop = FALSE]
  weights <- weights[ok]
  weights[!is.finite(weights) | weights <= 0] <- 1

  penalty_rows <- list()
  penalty_rhs <- numeric()
  for (cols in smooth_groups) {
    pen <- .capivara_second_derivative_penalty(length(cols), attr(cols, "lambda"))
    if (nrow(pen)) {
      block <- matrix(0, nrow = nrow(pen), ncol = ncol(design))
      block[, cols] <- pen
      penalty_rows[[length(penalty_rows) + 1L]] <- block
      penalty_rhs <- c(penalty_rhs, rep(0, nrow(block)))
    }
  }
  if (!is.null(ridge) && length(ridge)) {
    ridge <- ridge[is.finite(ridge) & ridge > 0]
    if (length(ridge)) {
      block <- diag(sqrt(ridge), nrow = length(ridge))
      full <- matrix(0, nrow = length(ridge), ncol = ncol(design))
      full[, as.integer(names(ridge))] <- block
      penalty_rows[[length(penalty_rows) + 1L]] <- full
      penalty_rhs <- c(penalty_rhs, rep(0, nrow(full)))
    }
  }

  solve_once <- function(extra_weights) {
    extra_weights[!is.finite(extra_weights) | extra_weights <= 0] <- 0
    w <- sqrt(weights * extra_weights)
    X <- design * w
    z <- y * w
    if (length(penalty_rows)) {
      X <- rbind(X, do.call(rbind, penalty_rows))
      z <- c(z, penalty_rhs)
    }
    as.numeric(qr.solve(X, z))
  }

  if (!isTRUE(robust)) {
    return(solve_once(rep(1, length(y))))
  }

  row_weights <- rep(1, length(y))
  coef <- solve_once(row_weights)
  status <- "robust_not_converged"
  for (iter in seq_len(as.integer(robust_maxit))) {
    resid <- y - as.numeric(design %*% coef)
    scale <- stats::mad(resid, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(scale) || scale <= 0) {
      scale <- stats::median(abs(resid), na.rm = TRUE) / 0.6745
    }
    if (!is.finite(scale) || scale <= 0) {
      status <- "robust_zero_scale"
      break
    }
    cutoff <- robust_k * scale
    new_row_weights <- pmin(1, cutoff / pmax(abs(resid), .Machine$double.eps))
    new_coef <- solve_once(new_row_weights)
    delta <- max(abs(new_coef - coef), na.rm = TRUE) / max(1, max(abs(coef), na.rm = TRUE))
    coef <- new_coef
    row_weights <- new_row_weights
    if (is.finite(delta) && delta < robust_tol) {
      status <- paste0("robust_huber_converged_", iter)
      break
    }
  }
  attr(coef, "robust_status") <- status
  attr(coef, "robust_weights") <- row_weights
  coef
}

#' Estimate the in-plane bar angle for a bisymmetric model
#'
#' @param spaxels Spaxel table after deprojection.
#' @param phi_b_deg Optional supplied in-plane bar angle relative to the major axis.
#' @param white_light Optional collapsed white-light image. When supplied with
#'   `geometry`, its central elongated light distribution provides an automatic
#'   photometric bar-angle prior.
#' @param geometry Disc geometry used to deproject the photometric angle.
#' @return A list with `phi_b_rad`, `phi_b_deg`, and status.
#' @noRd
estimate_bar_geometry <- function(spaxels,
                                  phi_b_deg = NULL,
                                  white_light = NULL,
                                  geometry = NULL) {
  if (!is.null(phi_b_deg) && length(phi_b_deg) == 1L && is.finite(phi_b_deg)) {
    phi <- as.numeric(phi_b_deg) * pi / 180
    return(list(phi_b_rad = phi, phi_b_deg = as.numeric(phi_b_deg), bar_status = "bar_angle_supplied"))
  }

  if (is.matrix(white_light) && is.list(geometry)) {
    axis <- .capivara_white_light_bar_axis(white_light, geometry)
    if (!is.null(axis)) {
      return(axis)
    }
  }

  bar <- spaxels$valid & tolower(spaxels$seg_class) == "bar" &
    is.finite(spaxels$X) & is.finite(spaxels$Yd)
  if (sum(bar) < 5L) {
    stop("Too few Capivara bar spaxels to estimate a bisymmetric bar angle.", call. = FALSE)
  }
  xy <- cbind(spaxels$X[bar], spaxels$Yd[bar])
  xy <- scale(xy, center = TRUE, scale = FALSE)
  eig <- eigen(stats::cov(xy), symmetric = TRUE)
  axis <- eig$vectors[, 1]
  phi <- atan2(axis[2], axis[1])
  phi <- ((phi + pi / 2) %% pi) - pi / 2
  list(phi_b_rad = phi, phi_b_deg = phi * 180 / pi, bar_status = "bar_angle_estimated_from_capivara_bar_support")
}

.capivara_white_light_bar_axis <- function(white_light,
                                           geometry,
                                           inner_radius_fraction = 0.55,
                                           brightness_quantile = 0.65,
                                           min_pixels = 20L) {
  required_geometry <- c("x0", "y0", "pa_rad", "inc_rad", "coordinate_convention")
  if (!all(required_geometry %in% names(geometry)) ||
      !all(is.finite(unlist(geometry[c("x0", "y0", "pa_rad", "inc_rad")]))) ||
      !any(is.finite(white_light))) {
    return(NULL)
  }

  # Match R's column-major matrix vectorisation: row/y changes fastest.
  grid <- expand.grid(
    y = seq_len(nrow(white_light)),
    x = seq_len(ncol(white_light))
  )
  grid <- grid[c("x", "y")]
  light <- as.vector(white_light)
  finite <- is.finite(light)
  if (sum(finite) < min_pixels) {
    return(NULL)
  }

  background <- stats::quantile(light[finite], 0.10, na.rm = TRUE, names = FALSE)
  signal <- pmax(light - background, 0)
  positive <- finite & signal > 0
  if (sum(positive) < min_pixels) {
    return(NULL)
  }

  projected <- deproject_coordinates(grid$x, grid$y, geometry)
  extent_cut <- stats::quantile(signal[positive], 0.20, na.rm = TRUE, names = FALSE)
  extent <- positive & signal >= extent_cut & is.finite(projected$R)
  r95 <- stats::quantile(projected$R[extent], 0.95, na.rm = TRUE, names = FALSE)
  if (!is.finite(r95) || r95 <= 0) {
    return(NULL)
  }

  central <- extent & projected$R <= inner_radius_fraction * r95
  if (sum(central) < min_pixels) {
    return(NULL)
  }
  light_cut <- stats::quantile(signal[central], brightness_quantile, na.rm = TRUE, names = FALSE)
  selected <- central & signal >= light_cut
  if (sum(selected) < min_pixels) {
    return(NULL)
  }

  weights <- signal[selected]
  xy <- cbind(grid$x[selected], grid$y[selected])
  centre <- colSums(xy * weights) / sum(weights)
  centred <- sweep(xy, 2L, centre, "-")
  covariance <- crossprod(centred * sqrt(weights)) / sum(weights)
  eig <- eigen(covariance, symmetric = TRUE)
  if (!all(is.finite(eig$values)) || eig$values[[2]] <= 0) {
    return(NULL)
  }

  projected_axis <- eig$vectors[, 1L]
  endpoint <- deproject_coordinates(
    c(geometry$x0, geometry$x0 + projected_axis[[1L]]),
    c(geometry$y0, geometry$y0 + projected_axis[[2L]]),
    geometry
  )
  phi <- endpoint$theta[[2L]]
  if (!is.finite(phi)) {
    return(NULL)
  }
  phi <- ((phi + pi / 2) %% pi) - pi / 2

  list(
    phi_b_rad = phi,
    phi_b_deg = phi * 180 / pi,
    bar_status = "bar_angle_estimated_from_inner_white_light",
    photometric_axis_ratio = sqrt(eig$values[[1L]] / eig$values[[2L]]),
    photometric_n_pixels = sum(selected),
    photometric_inner_radius = inner_radius_fraction * r95,
    photometric_brightness_quantile = brightness_quantile
  )
}

.capivara_velocity_weights <- function(spaxels, use_errors = TRUE) {
  if (isTRUE(use_errors) && "velocity_error" %in% names(spaxels)) {
    ok_err <- is.finite(spaxels$velocity_error) & spaxels$velocity_error > 0
    weights <- ifelse(ok_err, 1 / (spaxels$velocity_error^2 + 25), 1 / 25)
  } else {
    weights <- rep(1, nrow(spaxels))
  }
  if ("fit_weight" %in% names(spaxels)) {
    soft <- spaxels$fit_weight
    soft[!is.finite(soft) | soft <= 0] <- 1
    weights <- weights * soft
  }
  weights
}

.capivara_velocity_knots <- function(spaxels, n_rings) {
  Rmax <- stats::quantile(spaxels$R[spaxels$valid & is.finite(spaxels$R)], 0.98, na.rm = TRUE)
  seq(0, as.numeric(Rmax), length.out = max(3L, as.integer(n_rings)))
}

.capivara_bar_fit_mask <- function(spaxels) {
  spaxels$valid &
    is.finite(spaxels$velocity) &
    is.finite(spaxels$R) &
    is.finite(spaxels$theta) &
    !tolower(spaxels$seg_class) %in% c("disturbed", "unassigned")
}

#' Fit the projected bisymmetric velocity model used by the bar paper
#'
#' This implements the Section 3 line-of-sight model without PSF convolution:
#' Vlos = Vsys + sin(i) * [Vt cos(theta)
#'   - V2t cos(2(theta - phi_b)) cos(theta)
#'   - V2r sin(2(theta - phi_b)) sin(theta)].
#'
#' @param spaxels Data frame with valid spaxel coordinates, velocities,
#'   uncertainties, deprojected radius/theta, and segmentation class.
#' @param geometry Disc geometry list from `estimate_disc_geometry()`.
#' @param bar_geometry Bar geometry list from `estimate_bar_geometry()`.
#' @param n_rings Number of radial knots for the piecewise velocity profiles.
#' @param use_errors If `TRUE`, weight by velocity errors when available.
#' @param smooth_lambda First-order smoothing penalty for radial profiles.
#' @param second_order_lambda Curvature penalty for second-order bar terms.
#' @param fixed_vsys Optional fixed systemic velocity.
#' @param robust If `TRUE`, iteratively downweight large residuals.
#' @param robust_k Huber tuning constant used by the robust fit.
#' @param robust_maxit Maximum robust reweighting iterations.
#' @param max_v2_fraction Maximum allowed mean second-order amplitude relative
#'   to the circular term before shrinkage.
#' @param max_mean_v2 Absolute cap on mean second-order velocity amplitude.
#' @return A list with fitted spaxels, radial profiles, parameters, and
#'   diagnostics.
#' @noRd
fit_bisymmetric_model <- function(spaxels,
                                  geometry,
                                  bar_geometry,
                                  n_rings = 10L,
                                  use_errors = TRUE,
                                  smooth_lambda = 10,
                                  second_order_lambda = 500,
                                  fixed_vsys = NULL,
                                  robust = FALSE,
                                  robust_k = 1.345,
                                  robust_maxit = 20L,
                                  max_v2_fraction = 0.6,
                                  max_mean_v2 = 350) {
  fit_mask <- .capivara_bar_fit_mask(spaxels)
  dat <- spaxels[fit_mask, , drop = FALSE]
  if (nrow(dat) < max(20L, 3L * n_rings)) {
    stop("Too few valid spaxels for bisymmetric fit: ", nrow(dat), call. = FALSE)
  }

  knots <- .capivara_velocity_knots(dat, n_rings)
  basis <- .capivara_piecewise_basis(dat$R, knots)
  inc <- sin(geometry$inc_rad)
  phi <- bar_geometry$phi_b_rad
  theta <- dat$theta
  vt <- basis * (inc * cos(theta))
  v2t <- basis * (-inc * cos(2 * (theta - phi)) * cos(theta))
  v2r <- basis * (-inc * sin(2 * (theta - phi)) * sin(theta))
  vsys_fixed <- length(fixed_vsys) == 1L && is.finite(fixed_vsys)
  design <- if (vsys_fixed) {
    cbind(vt, v2t, v2r)
  } else {
    cbind(vsys = 1, vt, v2t, v2r)
  }

  n <- ncol(basis)
  offset <- if (vsys_fixed) 0L else 1L
  vt_cols <- (1 + offset):(offset + n)
  v2t_cols <- (1 + offset + n):(offset + 2 * n)
  v2r_cols <- (1 + offset + 2 * n):(offset + 3 * n)
  smooth_groups <- list(vt_cols, v2t_cols, v2r_cols)
  smooth_groups <- lapply(smooth_groups, function(x) {
    attr(x, "lambda") <- smooth_lambda
    x
  })
  ridge <- c(
    stats::setNames(rep(second_order_lambda, length(v2t_cols)), v2t_cols),
    stats::setNames(rep(second_order_lambda, length(v2r_cols)), v2r_cols)
  )
  coef <- .capivara_fit_linear_velocity_model(
    dat$velocity - if (vsys_fixed) fixed_vsys else 0,
    design,
    .capivara_velocity_weights(dat, use_errors),
    smooth_groups,
    ridge = ridge,
    robust = robust,
    robust_k = robust_k,
    robust_maxit = robust_maxit
  )

  full_basis <- .capivara_piecewise_basis(spaxels$R, knots)
  full_theta <- spaxels$theta
  vt_component <- as.numeric((full_basis * (inc * cos(full_theta))) %*% coef[vt_cols])
  v2t_component <- as.numeric((full_basis * (-inc * cos(2 * (full_theta - phi)) * cos(full_theta))) %*% coef[v2t_cols])
  v2r_component <- as.numeric((full_basis * (-inc * sin(2 * (full_theta - phi)) * sin(full_theta))) %*% coef[v2r_cols])
  vsys <- if (vsys_fixed) fixed_vsys else coef[1]
  pred <- vsys + vt_component + v2t_component + v2r_component
  pred[!spaxels$valid] <- NA_real_
  vt_component[!spaxels$valid] <- NA_real_
  v2t_component[!spaxels$valid] <- NA_real_
  v2r_component[!spaxels$valid] <- NA_real_

  spaxels$v_bar_model <- pred
  spaxels$v_t_component <- vt_component
  spaxels$v_2t_component <- v2t_component
  spaxels$v_2r_component <- v2r_component
  spaxels$v_bar_resid <- spaxels$velocity - pred
  spaxels$v_bar_resid[!spaxels$valid] <- NA_real_
  spaxels$v_disc <- pred
  spaxels$v_resid <- spaxels$v_bar_resid

  profile <- data.frame(
    R = knots,
    Vt = coef[vt_cols],
    V2t = coef[v2t_cols],
    V2r = coef[v2r_cols]
  )
  mean_v2 <- mean(sqrt(profile$V2t^2 + profile$V2r^2), na.rm = TRUE)
  vmax <- max(abs(profile$Vt), na.rm = TRUE)
  v2_ratio <- mean_v2 / max(vmax, .Machine$double.eps)
  fit_status <- "ok_bisymmetric_piecewise"
  robust_status <- attr(coef, "robust_status", exact = TRUE)
  if (isTRUE(robust) && !is.null(robust_status)) {
    fit_status <- paste(fit_status, robust_status, sep = ";")
  }
  if (is.finite(max_mean_v2) && is.finite(mean_v2) && mean_v2 > max_mean_v2) {
    fit_status <- paste(fit_status, "unstable_mean_V2", sep = ";")
  }
  if (is.finite(max_v2_fraction) && is.finite(v2_ratio) && v2_ratio > max_v2_fraction) {
    fit_status <- paste(fit_status, "unstable_V2_fraction", sep = ";")
  }
  list(
    spaxels = spaxels,
    parameters = data.frame(
      vsys = vsys,
      vmax = vmax,
      Rt = NA_real_,
      phi_b_deg = bar_geometry$phi_b_deg,
      mean_V2 = mean_v2,
      V2_over_Vt = v2_ratio,
      stringsAsFactors = FALSE
    ),
    profile = profile,
    knots = knots,
    bar_geometry = bar_geometry,
    robust_weights = attr(coef, "robust_weights", exact = TRUE),
    fit_status = fit_status,
    N_fit = nrow(dat)
  )
}

#' Fit a matching axisymmetric piecewise comparison model
#'
#' @param spaxels Data frame with valid spaxel coordinates, velocities,
#'   uncertainties, deprojected radius/theta, and segmentation class.
#' @param geometry Disc geometry list from `estimate_disc_geometry()`.
#' @param n_rings Number of radial knots for the piecewise velocity profile.
#' @param use_errors If `TRUE`, weight by velocity errors when available.
#' @param smooth_lambda First-order smoothing penalty for the radial profile.
#' @param fixed_vsys Optional fixed systemic velocity.
#' @param robust If `TRUE`, iteratively downweight large residuals.
#' @param robust_k Huber tuning constant used by the robust fit.
#' @param robust_maxit Maximum robust reweighting iterations.
#' @return A list with fitted spaxels, radial profile, parameters, and
#'   diagnostics.
#' @noRd
fit_axisymmetric_piecewise_model <- function(spaxels,
                                             geometry,
                                             n_rings = 10L,
                                             use_errors = TRUE,
                                             smooth_lambda = 10,
                                             fixed_vsys = NULL,
                                             robust = FALSE,
                                             robust_k = 1.345,
                                             robust_maxit = 20L) {
  fit_mask <- .capivara_bar_fit_mask(spaxels)
  dat <- spaxels[fit_mask, , drop = FALSE]
  knots <- .capivara_velocity_knots(dat, n_rings)
  basis <- .capivara_piecewise_basis(dat$R, knots)
  vsys_fixed <- length(fixed_vsys) == 1L && is.finite(fixed_vsys)
  design <- if (vsys_fixed) {
    basis * (sin(geometry$inc_rad) * cos(dat$theta))
  } else {
    cbind(vsys = 1, basis * (sin(geometry$inc_rad) * cos(dat$theta)))
  }
  cols <- if (vsys_fixed) seq_len(ncol(design)) else 2:ncol(design)
  attr(cols, "lambda") <- smooth_lambda
  coef <- .capivara_fit_linear_velocity_model(
    dat$velocity - if (vsys_fixed) fixed_vsys else 0,
    design,
    .capivara_velocity_weights(dat, use_errors),
    list(cols),
    ridge = NULL,
    robust = robust,
    robust_k = robust_k,
    robust_maxit = robust_maxit
  )

  full_basis <- .capivara_piecewise_basis(spaxels$R, knots)
  full_design <- full_basis * (sin(geometry$inc_rad) * cos(spaxels$theta))
  vsys <- if (vsys_fixed) fixed_vsys else coef[1]
  rot_coef <- if (vsys_fixed) coef else coef[cols]
  pred <- vsys + as.numeric(full_design %*% rot_coef)
  pred[!spaxels$valid] <- NA_real_
  spaxels$v_axisym_model <- pred
  spaxels$v_axisym_resid <- spaxels$velocity - pred
  spaxels$v_axisym_resid[!spaxels$valid] <- NA_real_

  fit_status <- "ok_axisymmetric_piecewise"
  robust_status <- attr(coef, "robust_status", exact = TRUE)
  if (isTRUE(robust) && !is.null(robust_status)) {
    fit_status <- paste(fit_status, robust_status, sep = ";")
  }

  list(
    spaxels = spaxels,
    parameters = data.frame(
      vsys = vsys,
      vmax = max(abs(rot_coef), na.rm = TRUE),
      Rt = NA_real_,
      stringsAsFactors = FALSE
    ),
    profile = data.frame(R = knots, Vt = rot_coef),
    knots = knots,
    robust_weights = attr(coef, "robust_weights", exact = TRUE),
    fit_status = fit_status,
    N_fit = nrow(dat)
  )
}
