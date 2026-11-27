#!/usr/bin/env Rscript
# Copy the installed viewer and generate a reproducible synthetic example.
# This script is executed after pkgdown::build_site().
args <- commandArgs(trailingOnly = TRUE)
site <- normalizePath(if (length(args)) args[[1L]] else "docs", mustWork = TRUE)
viewer <- system.file("atlas", "index.html", package = "capivara", mustWork = TRUE)
starter <- system.file("tutorials", "run_prepared_cube.R",
                       package = "capivara", mustWork = TRUE)
exporter <- system.file("tutorials", "export_semantic_atlas.R",
                        package = "capivara", mustWork = TRUE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
dest <- file.path(site, "atlas")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
out <- tempfile("capivara-site-demo-")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, shQuote(c(starter, "--demo", out, "8"))) == 0L)
stopifnot(system2(rscript, shQuote(c(exporter, out))) == 0L)
demo <- jsonlite::read_json(file.path(out, "atlas.json"), simplifyVector = FALSE)
stopifnot(identical(demo$schema, "capivara-semantic-atlas-v1"),
          identical(demo$object, "SYNTHETIC-DEMO"),
          identical(demo$source, "Synthetic demonstration"),
          length(demo$regions) == demo$actual_k)
stopifnot(file.copy(viewer, file.path(dest, "index.html"), overwrite = TRUE),
          file.copy(file.path(out, "atlas.json"), file.path(dest, "demo.json"),
                    overwrite = TRUE))
unlink(out, recursive = TRUE)
cat("Built atlas viewer and synthetic demo in ", dest, "\n", sep = "")
