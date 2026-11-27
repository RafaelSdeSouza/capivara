#!/usr/bin/env Rscript
# Generate the native multi-map atlas for a verified MaNGA LOGCUBE.
# Rscript run_atlas_views.R cube.fits plateifu redshift line output_dir [n_regions]
# Uses the existing native kinematics and SpectroPath implementation, not pPXF.
# Run on a local cube; does not download data or assume a disc/bar geometry.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(5L, 6L))
  stop("Usage: Rscript run_atlas_views.R cube.fits plateifu redshift line output_dir [n_regions]")
if (!requireNamespace("jsonlite", quietly = TRUE))
  stop("Install jsonlite with install.packages('jsonlite').")
suppressPackageStartupMessages(library(capivara))
cube <- normalizePath(args[1], mustWork = TRUE)
id <- args[2]
z <- suppressWarnings(as.numeric(args[3]))
line <- args[4]
out <- args[5]
k <- if (length(args) == 6L) suppressWarnings(as.numeric(args[6])) else 25
if (!grepl("^[0-9]+-[0-9]+$", id) ||
    !grepl(paste0("^manga-", id, "-LOGCUBE[.]fits([.]gz)?$"),
           basename(cube), ignore.case = TRUE))
  stop("Use the matching MaNGA LOGCUBE and its exact plate-IFU identifier.")
if (!is.finite(z) || z < 0 || !is.finite(k) || k < 1 ||
    k > .Machine$integer.max || k != floor(k))
  stop("Supply measured redshift and a positive integer region count.")
if (!line %in% list_emission_lines()$name)
  stop("Select a line from list_emission_lines().")
if (file.exists(out) || dir.exists(out)) stop("Choose a new output directory.")
dir.create(out, recursive = TRUE)
out <- normalizePath(out, mustWork = TRUE)
fit <- segment_kinematics(cube_path = cube, object_id = id, redshift = z,
  emission_line = line, segmentation_mode = "path_signature",
  output_dir = file.path(out, "native"), n_segments = as.integer(k),
  n_path_segments = as.integer(k), wavelength_frame = "observed",
  wavelength_medium = "vacuum", support_mode = "starlet", show_plots = FALSE)
n <- fit$native
kin <- n$kinematics
dims <- dim(n$kinematic_aware$cluster_map)
stopifnot(length(dims) == 2L, all(dims > 0), all(dim(kin$velocity) == dims))
valid <- kin$measured_valid & !kin$imputed
stopifnot(is.logical(valid), identical(dim(valid), dims), !anyNA(valid))
make_map <- function(id, title, group, kind, unit, mat, note) {
  stopifnot(is.matrix(mat), is.numeric(mat), identical(dim(mat), dims))
  mat[!is.finite(mat)] <- NA_real_
  if (kind == "categorical") {
    stopifnot(all(is.na(mat) | (mat > 0 & mat == floor(mat))))
  }
  list(id = id, title = title, group = group, kind = kind,
       unit = unit, note = note,
       values = lapply(seq_len(ncol(mat)), function(j) as.list(mat[, j])))
}
measured <- function(m) { m[!valid] <- NA_real_; m }
maps <- list(
  make_map("white_light", "White light", "overview", "continuous", "Native collapsed flux",
    n$star$collapsed, "Native white-light collapse; not an absolute photometric calibration."),
  make_map("line_flux", "Line flux", "emission", "continuous", "Native estimator units",
    measured(kin$flux), "Native single-line flux estimator; no fitted line ratios or covariance."),
  make_map("velocity", "Gas velocity", "kinematics", "diverging", "km/s",
    measured(kin$velocity), "Systemic-relative velocity from the selected line profile."),
  make_map("sigma", "Observed line width", "kinematics", "continuous", "km/s",
    measured(kin$sigma), "Observed profile dispersion; includes instrumental and pixel broadening."),
  make_map("kinematic_regions", "Kinematic segmentation", "kinematics", "categorical", "Region ID",
    n$kinematic_aware$cluster_map, "Regions based on native line flux, velocity, width and profile features."),
  make_map("asymmetry", "Line asymmetry", "morphology", "diverging", "Dimensionless",
    measured(kin$asymmetry), "Native red-blue profile asymmetry."),
  make_map("path_regions", "SpectroPath segmentation", "morphology", "categorical", "Region ID",
    n$path_signature$cluster_map, "Path-signature partition; labels are independent of kinematic region IDs.")
)
for (i in seq_along(n$path_features$features)) {
  feature <- as.character(n$path_features$features[[i]])
  stopifnot(grepl("^[A-Za-z0-9_]+$", feature))
  maps[[length(maps) + 1L]] <- make_map(
    paste0("path_", feature), paste("SpectroPath", feature), "morphology",
    "continuous", "Native descriptor", n$path_features$feature_cube[, , i],
    "Native path-signature coordinate; descriptive feature, not a physical parameter.")
}
bundle <- list(schema = "capivara-galaxy-views-v1", object = id,
  source = "User-supplied MaNGA LOGCUBE; native CAPIVARA measurement workflow",
  nx = dims[1], ny = dims[2], line = line, redshift = z, maps = maps,
  provenance = list(cube = basename(cube), cube_md5 = unname(tools::md5sum(cube)),
    package_version = as.character(utils::packageVersion("capivara")),
    redshift_source = fit$redshift_source, support = n$support_method,
    wavelength_frame = "observed", wavelength_medium = "vacuum",
    coordinate_system = "Same native array grid across bundle maps; no celestial WCS claim",
    requested_kinematic_regions = k, requested_path_regions = k,
    physical_models = "NOT_RUN", uncertainties = "Not calibrated in this workflow",
    width = "Observed profile width; no LSF correction",
    emission_line_segmentation = "Separate segment_emission_lines() analysis not included"))
jsonlite::write_json(bundle, file.path(out, "views.json"), auto_unbox = TRUE,
                    digits = NA, na = "null")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
message("Open the matching galaxy page and import ", file.path(out, "views.json"))
