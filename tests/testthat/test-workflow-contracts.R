test_that("measurement flags reach the same spaxels as their maps", {
  for (dims in list(c(2,3),c(3,3))) {
    sp <- expand.grid(x=seq_len(dims[2]),y=seq_len(dims[1]))
    sp$valid <- TRUE
    flag <- matrix(FALSE,dims[1],dims[2]); flag[2,1] <- TRUE
    kin <- list(measured_valid=!flag,imputed=flag)
    preview <- .capivara_measurement_flags(sp,kin,preview=TRUE)
    science <- .capivara_measurement_flags(sp,kin)
    expect_equal(which(preview$imputed),which(sp$x==1 & sp$y==2))
    expect_equal(preview$fit_weight[preview$imputed],.35)
    expect_false(science$valid[science$imputed])
    expect_equal(science$fit_weight[science$imputed],0)
  }
})

test_that("native FITS catalogue labels the actual array channels", {
  script <- system.file("extdata","kinematics","native_kinematics_workflow.R",
                        package="capivara",mustWork=TRUE)
  e <- new.env(parent=asNamespace("stats"))
  e$support_starlet <- matrix(1,2,3); e$support <- matrix(2,2,3)
  e$seg <- list(cluster_map=matrix(3,2,3)); e$kin_seg <- list(cluster_map=matrix(10,2,3))
  e$kin <- list(flux=matrix(1e4,2,3),velocity=matrix(5,2,3),sigma=matrix(6,2,3),
                asymmetry=matrix(7,2,3),h3_proxy=matrix(8,2,3),h4_proxy=matrix(9,2,3))
  e$line <- list(slug="halpha")
  for (ex in parse(script)) {
    if (is.call(ex) && identical(ex[[1]],as.name("<-")) &&
        (identical(ex[[2]],as.name("maps")) || startsWith(paste(deparse(ex[[2]]),collapse=""),"names(maps)"))) eval(ex,e)
  }
  expect_equal(names(e$maps),c("starlet_support","kinematic_support","capivara_segment",
    "halpha_log_flux","halpha_velocity_centered","halpha_sigma","halpha_asymmetry",
    "halpha_h3_proxy","halpha_h4_proxy","kinematic_aware_segment"))
  expect_equal(unname(vapply(e$maps,function(x)x[1],numeric(1))),as.numeric(1:10))
})

test_that("explicit checkout resolution outranks installed workflow scripts", {
  root <- tempfile(); dir.create(file.path(root,"inst","extdata","kinematics"),recursive=TRUE)
  f <- file.path(root,"inst","extdata","kinematics","test.R"); file.create(f)
  on.exit(unlink(root,recursive=TRUE))
  expect_equal(.capivara_workflow_file("test.R",root),normalizePath(f))
})

test_that("workflow environment cleanup is reversible on failure", {
  old <- Sys.getenv("CAPIVARA_MODEL_INC_DEG",unset=NA_character_)
  on.exit(.capivara_restore_env("CAPIVARA_MODEL_INC_DEG",old))
  Sys.setenv(CAPIVARA_MODEL_INC_DEG="77")
  run <- function() {
    before <- .capivara_clear_workflow_env("CAPIVARA_MODEL_")
    on.exit(.capivara_restore_env(names(before),before),add=TRUE)
    own <- c(CAPIVARA_MODEL_INC_DEG="35")
    prev <- .capivara_set_env(own)
    on.exit(.capivara_restore_env(names(own),prev),add=TRUE,after=FALSE)
    stop("deliberate failure")
  }
  expect_error(run(),"deliberate failure")
  expect_equal(Sys.getenv("CAPIVARA_MODEL_INC_DEG"),"77")
})

test_that("kinematic completion and robust convergence use actual fitter statuses", {
  fit <- list(fit_status="ok_bisymmetric_piecewise;robust_huber_converged_3",
              profile=data.frame(R=1:3,Vt=100:102))
  qc <- .capivara_kinematic_fit_qc(fit,"preview",TRUE)
  expect_true(qc$numerical_success)
  expect_true(qc$converged)
  expect_false(qc$scientifically_usable)
  fit$fit_status <- "ok_bisymmetric_piecewise;robust_not_converged;unstable_mean_V2"
  qc <- .capivara_kinematic_fit_qc(fit,"science",FALSE)
  expect_true(qc$numerical_success)
  expect_false(qc$converged)
  expect_true(qc$model_instability_flag)
  expect_false(qc$bar_support_available)
})
