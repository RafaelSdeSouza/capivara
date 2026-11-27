#!/usr/bin/env Rscript
# Copy the installed viewer and generate a reproducible synthetic example.
# This script is executed after pkgdown::build_site().
args <- commandArgs(trailingOnly = TRUE)
site <- normalizePath(if (length(args)) args[[1L]] else "docs", mustWork = TRUE)
viewer <- system.file("atlas", "index.html", package = "capivara", mustWork = TRUE)
gallery <- system.file("atlas", "gallery.html", package = "capivara", mustWork = TRUE)
starter <- system.file("tutorials", "run_prepared_cube.R",
                       package = "capivara", mustWork = TRUE)
exporter <- system.file("tutorials", "export_semantic_atlas.R",
                        package = "capivara", mustWork = TRUE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install jsonlite.")
# A complete root build is required so Atlas navigation and articles share a root.
stopifnot(file.exists(file.path(site, "articles", "your-own-data.html")),
          file.exists(file.path(site, "articles", "canonical-api.html")))
# GitHub README keeps source-file links. Rewrite only their generated home-page
# equivalents to the built website's actual article destinations.
home <- file.path(site, "index.html")
html <- readLines(home, warn = FALSE)
links <- c("vignettes/your-own-data.Rmd" = "articles/your-own-data.html",
           "examples/README.html" = "articles/examples.html",
           "examples/README.md" = "articles/examples.html")
for (from in names(links)) {
  html <- gsub(paste0('href="', from, '"'),
               paste0('href="', links[[from]], '"'), html, fixed = TRUE)
}
# Quarto's GFM output can drop explicit source heading IDs. Preserve the
# historical About deep link by adding an anchor at the citation section.
if (!any(grepl('id="about"', html, fixed = TRUE))) {
  doc <- paste(html, collapse = "\n")
  pattern <- '(?is)(<h[2-4]\\b[^>]*\\bid="cite-capivara"[^>]*>)'
  if (!grepl(pattern, doc, perl = TRUE))
    stop("Cannot locate the citation heading for the legacy About anchor.")
  doc <- sub(pattern, '<span id="about"></span>\\1', doc, perl = TRUE)
  html <- strsplit(doc, "\n", fixed = TRUE)[[1L]]
}
writeLines(html, home, useBytes = TRUE)
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
stopifnot(file.copy(gallery, file.path(dest, "gallery.html"), overwrite = TRUE),
          file.copy(viewer, file.path(dest, "index.html"), overwrite = TRUE),
          file.copy(file.path(out, "atlas.json"), file.path(dest, "demo.json"),
                    overwrite = TRUE))
unlink(out, recursive = TRUE)
cat("Built atlas viewer and synthetic demo in ", dest, "\n", sep = "")
