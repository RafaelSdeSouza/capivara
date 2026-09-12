.capivara_measurement_flags <- function(spaxels, kin, preview = FALSE) {
  index <- cbind(spaxels$y, spaxels$x)
  spaxels$measured_valid <- kin$measured_valid[index]
  spaxels$imputed <- kin$imputed[index]
  spaxels$fit_weight <- ifelse(!is.na(spaxels$imputed) & spaxels$imputed,
                             if (preview) 0.35 else 0, 1)
  if (!preview) spaxels$valid <- spaxels$valid & !is.na(spaxels$measured_valid) & spaxels$measured_valid
  spaxels
}

.capivara_kinematic_fit_qc <- function(fit, analysis_mode, bar_support_available) {
  status <- fit$fit_status
  completed <- is.character(status) && length(status)==1L && startsWith(status,"ok_")
  finite <- is.data.frame(fit$profile) && nrow(fit$profile)>0L &&
    all(is.finite(as.matrix(fit$profile)))
  converged <- completed && finite && !grepl("robust_not_converged|robust_zero_scale",status)
  list(numerical_success=completed && finite, converged=converged,
       model_instability_flag=completed && grepl("unstable_",status),
       scientifically_usable=FALSE,
       reason=if(analysis_mode=="preview") "preview_not_for_scientific_inference" else
         "velocity_uncertainty_and_model_recovery_not_validated",
       bar_support_available=bar_support_available)
}
