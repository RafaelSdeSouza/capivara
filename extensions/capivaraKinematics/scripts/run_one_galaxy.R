#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs()[grep("^--file=", commandArgs())]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "extensions/capivaraKinematics/scripts/run_one_galaxy.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(root, "R"))) {
  root <- normalizePath("extensions/capivaraKinematics", mustWork = TRUE)
}
config_file <- if (length(args)) args[[1]] else file.path(root, "configs/example_one_galaxy.yml")

invisible(lapply(sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)), source))

run_one_galaxy(config_file)
