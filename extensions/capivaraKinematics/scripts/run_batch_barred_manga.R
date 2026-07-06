#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "extensions/capivaraKinematics/scripts/run_batch_barred_manga.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "R"))) {
  root <- normalizePath("extensions/capivaraKinematics", mustWork = TRUE)
}
batch_table_file <- if (length(args)) args[[1]] else file.path(root, "configs/example_batch_table.csv")
output_root <- if (length(args) >= 2L) args[[2]] else "outputs"

invisible(lapply(sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)), source))

run_batch_barred_manga(batch_table_file, output_root = output_root)
