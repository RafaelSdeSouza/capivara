#' Fit an axisymmetric disc model to Capivara-labelled clean disc spaxels
#'
#' Capivara segmentation defines the structural regions of the galaxy. The
#' axisymmetric disc model is fitted only to Capivara-labelled clean disc
#' spaxels. The fitted disc model is then extrapolated over the full IFU field
#' and subtracted from the observed velocity map.
#'
#' @param spaxels Spaxel data frame with velocity, R, theta, valid, and segment
#'   columns.
#' @param geometry Geometry list.
#' @param use_errors Use velocity errors as inverse-variance weights.
#' @param min_disc_spaxels Minimum number of clean disc spaxels.
#' @param fixed Optional named list of fixed parameters among vsys, vmax, Rt.
#' @return A list with updated spaxels, fitted parameters, covariance, and status.
#' @export
fit_disc_model <- function(spaxels,
                           geometry,
                           use_errors = TRUE,
                           min_disc_spaxels = 12L,
                           fixed = list()) {
  fit_mask <- spaxels$valid &
    (spaxels$use_for_disc_fit %in% TRUE) &
    is.finite(spaxels$velocity) &
    is.finite(spaxels$R) &
    is.finite(spaxels$theta)

  disc_spaxels <- spaxels[fit_mask, , drop = FALSE]
  if (nrow(disc_spaxels) < min_disc_spaxels) {
    stop("Too few clean disc spaxels for disc fit: ", nrow(disc_spaxels), call. = FALSE)
  }

  if (isTRUE(use_errors) && "velocity_error" %in% names(disc_spaxels)) {
    ok_err <- is.finite(disc_spaxels$velocity_error) & disc_spaxels$velocity_error > 0
    weights <- ifelse(ok_err, 1 / disc_spaxels$velocity_error^2, 1)
  } else {
    weights <- rep(1, nrow(disc_spaxels))
  }

  start <- list(
    vsys = stats::median(disc_spaxels$velocity, na.rm = TRUE),
    vmax = 0.5 * diff(range(disc_spaxels$velocity, na.rm = TRUE)),
    Rt = max(1, stats::median(disc_spaxels$R, na.rm = TRUE) / 2)
  )
  if (is.finite(geometry$vsys)) {
    start$vsys <- geometry$vsys
  }
  is_fixed_finite <- function(x) length(x) == 1L && is.finite(x)
  for (nm in names(fixed)) {
    if (nm %in% names(start) && is_fixed_finite(fixed[[nm]])) {
      start[[nm]] <- fixed[[nm]]
    }
  }

  fixed_names <- names(fixed)[vapply(fixed, is_fixed_finite, logical(1))]
  free <- setdiff(names(start), fixed_names)
  fit_status <- "ok"
  cov <- NULL

  if (!length(free)) {
    par <- unlist(start)
  } else if (requireNamespace("minpack.lm", quietly = TRUE)) {
    dat <- disc_spaxels
    dat$w <- sqrt(weights)
    formula <- stats::as.formula(
      paste0(
        "w * velocity ~ w * (",
        if ("vsys" %in% free) "vsys" else format(start$vsys, scientific = FALSE),
        " + sin(inc_rad) * (",
        if ("vmax" %in% free) "vmax" else format(start$vmax, scientific = FALSE),
        " * (2 / pi) * atan(R / ",
        if ("Rt" %in% free) "Rt" else format(start$Rt, scientific = FALSE),
        ")) * cos(theta))"
      )
    )
    lower <- c(vsys = -Inf, vmax = -Inf, Rt = 0.1)[free]
    upper <- c(vsys = Inf, vmax = Inf, Rt = Inf)[free]
    fit <- try(
      minpack.lm::nlsLM(
        formula,
        data = transform(dat, inc_rad = geometry$inc_rad),
        start = start[free],
        lower = lower,
        upper = upper,
        control = minpack.lm::nls.lm.control(maxiter = 200)
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      fit_status <- paste("failed:", conditionMessage(attr(fit, "condition")))
      stop(fit_status, call. = FALSE)
    }
    par <- unlist(start)
    par[names(stats::coef(fit))] <- stats::coef(fit)
    cov <- try(stats::vcov(fit), silent = TRUE)
    if (inherits(cov, "try-error")) {
      cov <- NULL
    }
  } else {
    objective <- function(p_free) {
      par_try <- unlist(start)
      par_try[free] <- p_free
      pred <- disc_velocity_model(par_try, disc_spaxels$R, disc_spaxels$theta, geometry$inc_rad)
      sum(weights * (disc_spaxels$velocity - pred)^2, na.rm = TRUE)
    }
    lower <- c(vsys = -Inf, vmax = -Inf, Rt = 0.1)[free]
    upper <- c(vsys = Inf, vmax = Inf, Rt = Inf)[free]
    opt <- try(
      stats::optim(
        par = unlist(start[free]),
        fn = objective,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = 1000)
      ),
      silent = TRUE
    )
    if (inherits(opt, "try-error")) {
      fit_status <- paste("failed:", conditionMessage(attr(opt, "condition")))
      stop(fit_status, call. = FALSE)
    }
    par <- unlist(start)
    par[free] <- opt$par
    if (!identical(opt$convergence, 0L)) {
      fit_status <- paste0("optim_convergence_", opt$convergence)
    } else {
      fit_status <- "ok_optim"
    }
  }

  spaxels$v_disc <- disc_velocity_model(par, spaxels$R, spaxels$theta, geometry$inc_rad)
  spaxels$v_disc[!spaxels$valid] <- NA_real_
  spaxels$v_resid <- spaxels$velocity - spaxels$v_disc
  spaxels$v_resid[!spaxels$valid] <- NA_real_
  if (any(!is.finite(spaxels$v_resid[spaxels$valid]))) {
    fit_status <- paste(fit_status, "nonfinite_residuals", sep = ";")
  }

  list(
    spaxels = spaxels,
    parameters = data.frame(
      vsys = unname(par[["vsys"]]),
      vmax = unname(par[["vmax"]]),
      vmax_abs = abs(unname(par[["vmax"]])),
      Rt = unname(par[["Rt"]]),
      stringsAsFactors = FALSE
    ),
    covariance = cov,
    fit_status = fit_status,
    N_disc_fit = nrow(disc_spaxels)
  )
}
