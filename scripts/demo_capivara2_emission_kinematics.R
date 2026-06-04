#!/usr/bin/env Rscript

repo_dir <- normalizePath(
  Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara"),
  mustWork = TRUE
)

path_root <- normalizePath(
  Sys.getenv(
    "CAPIVARA2_PATH_ROOT",
    unset = "/Users/rd23aag/Documents/GitHub/spectropath_paper_workspace_2026-05-08/outputs/ifu_gas_path_maps"
  ),
  mustWork = TRUE
)

include_oiii <- tolower(Sys.getenv("CAPIVARA2_DEMO_INCLUDE_OIII", unset = "false")) %in% c("1", "true", "yes")

cases <- data.frame(
  object = c("manga-7443-12703-LOGCUBE", "manga-10224-6104-LOGCUBE"),
  line = c("ifu_halpha_nii_complex", "ifu_halpha_nii_complex"),
  stringsAsFactors = FALSE
)

if (include_oiii) {
  cases <- rbind(
    cases,
    data.frame(
      object = c("manga-7443-12703-LOGCUBE", "manga-10224-6104-LOGCUBE"),
      line = c("ifu_oiii5007", "ifu_oiii5007"),
      stringsAsFactors = FALSE
    )
  )
}

prototype <- file.path(repo_dir, "scripts", "prototype_capivara2_kinematics.R")

for (i in seq_len(nrow(cases))) {
  run_dir <- file.path(path_root, cases$object[i], cases$line[i])
  out_dir <- file.path(repo_dir, "outputs", "capivara2_kinematics_prototype", cases$object[i], cases$line[i])

  message("Running ", cases$object[i], " / ", cases$line[i])
  Sys.setenv(
    CAPIVARA_REPO = repo_dir,
    CAPIVARA2_PATH_RUN = run_dir,
    CAPIVARA2_OUT = out_dir,
    CAPIVARA2_NCOMP = Sys.getenv("CAPIVARA2_NCOMP", unset = "50"),
    CAPIVARA2_KNN = Sys.getenv("CAPIVARA2_KNN", unset = "40"),
    CAPIVARA2_SPATIAL_WEIGHT = Sys.getenv("CAPIVARA2_SPATIAL_WEIGHT", unset = "0.15"),
    CAPIVARA2_FEATURES = Sys.getenv("CAPIVARA2_FEATURES", unset = "p3u,p3F,p4T"),
    CAPIVARA2_SUPPORT_MODE = Sys.getenv("CAPIVARA2_SUPPORT_MODE", unset = "spatial"),
    CAPIVARA2_FEATURE_IMPUTE = Sys.getenv("CAPIVARA2_FEATURE_IMPUTE", unset = "neutral")
  )

  sys.source(prototype, envir = new.env(parent = globalenv()))
}

message("Done. Outputs are under:")
message(file.path(repo_dir, "outputs", "capivara2_kinematics_prototype"))
