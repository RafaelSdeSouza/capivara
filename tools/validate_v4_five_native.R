#!/usr/bin/env Rscript

# Read-only production replay of the Experiment-X native availability audit.
# This script reads frozen Experiment-VIII arrays and never constructs a native
# partition. Run it against an isolated installation of CAPIVARA V4.

if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("The read-only native audit requires the suggested package `reticulate`.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The read-only native audit requires `jsonlite` to read the frozen oracle summary.")
}
suppressPackageStartupMessages(library(capivara))

project_root <- normalizePath(getwd(), mustWork = TRUE)
evidence_root <- Sys.getenv(
  "CAPIVARA_FROZEN_EVIDENCE",
  file.path(dirname(project_root), "capivara-hierarchical-segmentation",
            "research", "semantic_representations")
)
experiment8 <- file.path(evidence_root, "experiment8", "native")
experiment10 <- file.path(evidence_root, "experiment10")
if (!dir.exists(experiment8) || !file.exists(file.path(experiment10, "native", "summary.json"))) {
  stop("Frozen Experiment-VIII/X evidence was not found at `", evidence_root, "`.")
}

np <- reticulate::import("numpy", convert = FALSE)
npz <- function(path, keys) {
  z <- np$load(path)
  on.exit(z$close(), add = TRUE)
  stats::setNames(lapply(keys, function(key) reticulate::py_to_r(z$get(key))), keys)
}

profile <- manga_sandra_spectral_shape()
expected_windows <- list(c(5050, 5100), c(5400, 5500),
                         c(6000, 6100), c(6800, 6900))
stopifnot(
  identical(profile$required_support$amplitude_windows, expected_windows),
  identical(profile$coordinate_domain$frame, "rest"),
  identical(profile$coordinate_domain$medium, "vacuum"),
  identical(profile$eligibility$measured_amplitude_snr_min, 30)
)

oracle <- jsonlite::fromJSON(file.path(experiment10, "native", "summary.json"),
                             simplifyVector = FALSE)
oracle_by_id <- stats::setNames(oracle$galaxies,
                               vapply(oracle$galaxies, `[[`, character(1), "galaxy"))
ids <- c("7443-12703", "8443-6102", "10224-6104",
         "8135-12701", "8602-12705")

metric_fields <- c(
  fraction_requested_support_eligible = "spatial_fraction",
  fraction_positive_measured_continuum_retained = "positive_continuum_fraction",
  outer_quartile_fraction = "outer_quartile_fraction",
  r90_retention_ratio = "r90_ratio"
)
rows <- vector("list", length(ids))

for (i in seq_along(ids)) {
  id <- ids[i]
  amplitude <- npz(file.path(experiment8, id, "amplitude_samples.npz"),
                   c("flux", "variance", "usable", "wave_rest"))
  maps <- npz(file.path(experiment8, id, "eligibility_maps.npz"),
              c("A", "sigma", "snr", "complete", "positive", "feature",
                "overlap", "support", "continuum", "radius"))
  spatial_dim <- dim(maps$support)
  n_spaxel <- prod(spatial_dim)
  flux <- matrix(amplitude$flux, nrow = n_spaxel)
  variance <- matrix(amplitude$variance, nrow = n_spaxel)
  usable <- matrix(amplitude$usable, nrow = n_spaxel)
  complete <- rowSums(!usable) == 0L &
    rowSums(!is.finite(flux) | !is.finite(variance) | variance <= 0) == 0L
  measured_A <- rep(NA_real_, n_spaxel)
  measured_sigma <- measured_A
  measured_A[complete] <- rowMeans(flux[complete, , drop = FALSE])
  measured_sigma[complete] <-
    sqrt(rowSums(variance[complete, , drop = FALSE])) / ncol(variance)
  measured_snr <- measured_A / measured_sigma

  frozen_complete <- as.vector(maps$complete)
  stopifnot(
    identical(complete, frozen_complete),
    isTRUE(all.equal(measured_A[complete], as.vector(maps$A)[complete],
                     tolerance = 5e-15, check.attributes = FALSE)),
    isTRUE(all.equal(measured_sigma[complete], as.vector(maps$sigma)[complete],
                     tolerance = 5e-15, check.attributes = FALSE)),
    isTRUE(all.equal(measured_snr[complete], as.vector(maps$snr)[complete],
                     tolerance = 5e-13, check.attributes = FALSE))
  )

  # Exercise the declared public functional itself on deterministic measured
  # spaxels; the vectorized calculation above supplies the exhaustive replay.
  sampled <- head(which(complete & as.vector(maps$support)), 7L)
  for (index in sampled) {
    value <- profile$amplitude_functional$evaluator(
      flux[index, ], variance[index, ], usable[index, ], amplitude$wave_rest
    )
    stopifnot(
      isTRUE(value$complete),
      isTRUE(all.equal(value$value, measured_A[index], tolerance = 5e-15)),
      isTRUE(all.equal(sqrt(value$variance), measured_sigma[index],
                       tolerance = 5e-15))
    )
  }

  eligible <- matrix(
    complete & measured_A > 0 & measured_snr >= 30 &
      as.vector(maps$feature) & as.vector(maps$overlap),
    spatial_dim[1], spatial_dim[2]
  )
  base <- maps$support
  support <- build_capivara_support(
    base, base, base, centre = c(NA_real_, NA_real_),
    construction_method = "frozen Experiment-VIII native support",
    source = "read-only production V4 audit",
    provenance = list(experiment8_commit = "9c23b7d4",
                      experiment10_commit = "5335c878", galaxy = id)
  )
  prepared <- structure(
    list(eligible_map = eligible, measured_continuum = maps$continuum,
         representation = profile),
    class = c("capivara_prepared_representation", "list")
  )
  coverage <- representation_coverage_qc(
    prepared, support, continuum = maps$continuum, radius = maps$radius,
    radius_units = "arcsec"
  )
  expected <- oracle_by_id[[id]]
  stopifnot(
    identical(coverage$requested_support, as.integer(expected$support)),
    identical(coverage$eligible, as.integer(expected$eligible)),
    identical(coverage$applicability, expected$applicability)
  )
  for (field in names(metric_fields)) {
    stopifnot(isTRUE(all.equal(coverage[[field]], expected[[metric_fields[[field]]]],
                               tolerance = 5e-13)))
  }
  rows[[i]] <- data.frame(
    galaxy = id, requested_support = coverage$requested_support,
    eligible = coverage$eligible,
    eligible_fraction = coverage$fraction_requested_support_eligible,
    positive_continuum_fraction =
      coverage$fraction_positive_measured_continuum_retained,
    outer_quartile_fraction = coverage$outer_quartile_fraction,
    r90_ratio = coverage$r90_retention_ratio,
    applicability = coverage$applicability,
    amplitude_replay = "EXACT_WITHIN_BINARY64_TOLERANCE",
    native_partition_generated = FALSE,
    stringsAsFactors = FALSE
  )
}

# Complementary 8602 regression: the recovered measured-light aperture must
# retain the exact Experiment-X availability result at S/N >= 30.
maps8602 <- npz(file.path(experiment8, "8602-12705", "eligibility_maps.npz"),
                c("complete", "positive", "snr", "feature", "overlap",
                  "continuum", "radius"))
accounting8602 <- npz(file.path(experiment8, "8602-12705", "frozen_accounting.npz"),
                      c("retained_measured_light"))
eligible8602 <- maps8602$complete & maps8602$positive & maps8602$snr >= 30 &
  maps8602$feature & maps8602$overlap
recovered <- accounting8602$retained_measured_light
support8602 <- build_capivara_support(
  recovered, recovered, recovered,
  construction_method = "frozen Experiment-VIII recovered-light aperture",
  source = "read-only production V4 audit",
  provenance = list(experiment8_commit = "9c23b7d4", galaxy = "8602-12705")
)
prepared8602 <- structure(
  list(eligible_map = eligible8602, measured_continuum = maps8602$continuum,
       representation = profile),
  class = c("capivara_prepared_representation", "list")
)
recovered_qc <- representation_coverage_qc(
  prepared8602, support8602, maps8602$continuum, maps8602$radius, "arcsec"
)
stopifnot(
  identical(recovered_qc$requested_support,
            as.integer(oracle$recovered_333$support)),
  identical(recovered_qc$eligible, as.integer(oracle$recovered_333$eligible)),
  isTRUE(all.equal(recovered_qc$fraction_positive_measured_continuum_retained,
                   oracle$recovered_333$positive_continuum_fraction,
                   tolerance = 5e-13)),
  isTRUE(all.equal(recovered_qc$r90_retention_ratio,
                   oracle$recovered_333$r90_ratio, tolerance = 5e-13))
)

result <- do.call(rbind, rows)
out_dir <- file.path(project_root, "validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(result, file.path(out_dir, "v4_five_native_coverage.csv"),
                 row.names = FALSE)
writeLines(c(
  paste("package_version:", as.character(utils::packageVersion("capivara"))),
  paste("package_library:", find.package("capivara")),
  paste("evidence_root:", normalizePath(evidence_root)),
  "experiment_viii_commit: 9c23b7d4",
  "experiment_ix_commit: 19d0c60a",
  "experiment_x_commit: 5335c878",
  paste("profile_id:", profile$representation_id),
  "amplitude_replay: EXACT_WITHIN_BINARY64_TOLERANCE",
  "coverage_replay: EXACT_WITHIN_BINARY64_TOLERANCE",
  paste("recovered_8602_eligible:", recovered_qc$eligible),
  paste("recovered_8602_support:", recovered_qc$requested_support),
  "native_partitions_generated: 0",
  "verdict: READ_ONLY_FIVE_OBJECT_AUDIT_PASSED"
), file.path(out_dir, "v4_five_native_audit.txt"))
print(result, row.names = FALSE)
cat("8602 recovered aperture:", recovered_qc$eligible, "/",
    recovered_qc$requested_support, "eligible; positive continuum",
    sprintf("%.6f", recovered_qc$fraction_positive_measured_continuum_retained),
    "; r90 ratio", sprintf("%.6f", recovered_qc$r90_retention_ratio), "\n")
cat("READ_ONLY_FIVE_OBJECT_AUDIT_PASSED\n")
