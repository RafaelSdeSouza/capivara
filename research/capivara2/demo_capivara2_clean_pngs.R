#!/usr/bin/env Rscript

repo_dir <- normalizePath(
  Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara"),
  mustWork = TRUE
)

Sys.setenv(
  CAPIVARA_REPO = repo_dir,
  CAPIVARA2_MAKE_MOSAICS = Sys.getenv("CAPIVARA2_MAKE_MOSAICS", unset = "false"),
  CAPIVARA2_COPY_TO_PAPER = Sys.getenv("CAPIVARA2_COPY_TO_PAPER", unset = "false"),
  CAPIVARA2_INCLUDE_OIII = Sys.getenv("CAPIVARA2_INCLUDE_OIII", unset = "false"),
  CAPIVARA2_MAKE_STARLET = Sys.getenv("CAPIVARA2_MAKE_STARLET", unset = "false"),
  CAPIVARA2_MAKE_KINEMATICS = Sys.getenv("CAPIVARA2_MAKE_KINEMATICS", unset = "true"),
  CAPIVARA2_MAKE_PPXF = Sys.getenv("CAPIVARA2_MAKE_PPXF", unset = "true")
)

source(file.path(repo_dir, "scripts", "make_capivara2_paper_figures.R"))

message("Clean individual PNG/PDF files are in:")
message(Sys.getenv(
  "CAPIVARA2_FIG_OUT",
  unset = file.path(repo_dir, "outputs", "capivara2_paper_figures", "individual")
))
