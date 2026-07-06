script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "extensions/capivaraKinematics/scripts/run_8078_native_starletfull.R")
repo <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), mustWork = TRUE)
runner <- file.path(repo, "scripts", "run_manga_10218_bar_merger_maps.R")
cube <- "/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/normal_bar/manga-8078-12703-LOGCUBE.fits"
out_dir <- Sys.getenv("CAPIVARA_8078_NATIVE_OUT", unset = "/private/tmp/capivara_native_8078_starletfull")

env <- c(
  CAPIVARA_REDSHIFT = "0.0281",
  CAPIVARA_LINE = "halpha",
  CAPIVARA_OUTPUT_PREFIX = "manga8078_native_bettermask",
  CAPIVARA_KNN = "100",
  CAPIVARA_NCOMP = "25",
  CAPIVARA_PATH_KNN = "100",
  CAPIVARA_PATH_NCOMP = "45",
  CAPIVARA_PATH_SPATIAL_WEIGHT = "0.10",
  CAPIVARA_STARLET_SCALES = "1:5",
  CAPIVARA_STARLET_INCLUDE_COARSE = "true",
  CAPIVARA_LINE_WINDOW_KMS = "600",
  CAPIVARA_PEAK_SEARCH_KMS = "350",
  CAPIVARA_CENTROID_WINDOW_KMS = "220"
)

if (!file.exists(cube)) {
  stop("Missing cube: ", cube, call. = FALSE)
}
if (!file.exists(runner)) {
  stop("Missing native runner: ", runner, call. = FALSE)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
status <- system2("Rscript", c(runner, cube, out_dir), env = paste(names(env), env, sep = "="))
if (!identical(status, 0L)) {
  stop("Native 8078 run failed with status ", status, call. = FALSE)
}

message("Native 8078 products written to: ", out_dir)
