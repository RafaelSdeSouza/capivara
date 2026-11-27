#!/usr/bin/env Rscript
# Export a compact atlas from run_prepared_cube.R results.
# Usage: Rscript export_semantic_atlas.R analysis_directory
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript export_semantic_atlas.R analysis_directory")
paths <- capivara::export_atlas(args[1L])
message("Open the Galaxy Atlas and choose ",paths[["atlas"]],".")
