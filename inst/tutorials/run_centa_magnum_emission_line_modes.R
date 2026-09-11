#!/usr/bin/env Rscript

# Cen A MAGNUM: compare no-mask emission-line segmentation modes.
# This uses redshifted AGN/emission-line windows instead of the full spectrum.

suppressPackageStartupMessages({
  library(FITSio)
  library(ggplot2)
  library(patchwork)
  library(pkgload)
})

capivara_repo <- Sys.getenv("CAPIVARA_REPO", unset = "")
if (nzchar(capivara_repo)) pkgload::load_all(capivara_repo, quiet = TRUE) else library(capivara)

cube_path <- Sys.getenv(
  "CAPIVARA_CUBE_PATH",
  unset = ""
)
if (!nzchar(cube_path)) stop("Set the explicit cube-path environment variable before running this tutorial.")
output_dir <- Sys.getenv(
  "CAPIVARA_OUTPUT_DIR",
  unset = ""
)
if (!nzchar(output_dir)) stop("Set CAPIVARA_OUTPUT_DIR explicitly.")

redshift <- as.numeric(Sys.getenv("CAPIVARA_REDSHIFT", unset = "0.00183"))
n_segments <- as.integer(Sys.getenv("CAPIVARA_NCOMP", unset = "50"))
knn_k <- as.integer(Sys.getenv("CAPIVARA_KNN", unset = "20"))
max_pixels <- as.integer(Sys.getenv("CAPIVARA_MAX_PIXELS", unset = "18000"))
line_window_kms <- as.numeric(Sys.getenv("CAPIVARA_LINE_WINDOW_KMS", unset = "500"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading cube: ", cube_path)
cube <- FITSio::readFITS(cube_path)

run_one <- function(mode) {
  message("Running emission-line mode: ", mode)
  seg <- capivara::segment_emission_lines(
    cube,
    redshift = redshift,
    lines = "agn",
    Ncomp = n_segments,
    knn_k = knn_k,
    feature_mode = mode,
    line_window_kms = line_window_kms,
    spatial_weight = 0.05,
    max_pixels = max_pixels,
    seed = 20260730,
    verbose = TRUE
  )

  prefix <- file.path(output_dir, paste0("centa_emission_", mode, "_n", n_segments))
  seg_to_save <- seg
  seg_to_save$original_cube <- NULL
  saveRDS(seg_to_save, paste0(prefix, "_result.rds"))
  FITSio::writeFITSim(seg$cluster_map * 1, file = paste0(prefix, "_cluster_map.fits"), type = "double")

  p <- capivara::plot_cluster(seg) +
    ggplot2::labs(title = mode) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13))
  ggplot2::ggsave(paste0(prefix, "_segmentation.png"), p, width = 6, height = 6, dpi = 300, bg = "white")
  rm(seg_to_save)
  seg$original_cube <- NULL
  gc()

  list(mode = mode, seg = seg, plot = p, prefix = prefix)
}

modes <- c("windows", "window_derivative", "profile")
results <- lapply(modes, run_one)

panel <- (results[[1]]$plot | results[[2]]$plot | results[[3]]$plot) +
  patchwork::plot_annotation(
    title = "Cen A MAGNUM: no-mask emission-line Capivara segmentation",
    subtitle = sprintf("z = %.5f, AGN line windows, N = %d, k = %d, max_pixels = %d", redshift, n_segments, knn_k, max_pixels)
  )
ggplot2::ggsave(
  file.path(output_dir, "centa_emission_line_mode_comparison.png"),
  panel,
  width = 15,
  height = 5.4,
  dpi = 300,
  bg = "white"
)

summary <- do.call(rbind, lapply(results, function(x) {
  info <- x$seg$backend_info
  data.frame(
    mode = x$mode,
    requested_segments = n_segments,
    actual_segments = x$seg$Ncomp,
    valid_pixels = info$valid_pixels,
    sampled_pixels = if (!is.null(info$sampled_pixels)) info$sampled_pixels else NA_integer_,
    sampled_projection = isTRUE(info$sampled_projection),
    knn_k = info$knn_k,
    redshift = redshift,
    line_window_kms = line_window_kms,
    result_rds = paste0(x$prefix, "_result.rds"),
    cluster_map_fits = paste0(x$prefix, "_cluster_map.fits"),
    segmentation_png = paste0(x$prefix, "_segmentation.png"),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(summary, file.path(output_dir, "centa_emission_line_mode_summary.csv"), row.names = FALSE)

used_lines <- results[[1]]$seg$emission_line_features$lines
utils::write.csv(used_lines, file.path(output_dir, "centa_emission_line_windows_used.csv"), row.names = FALSE)

writeLines(
  c(
    "Cen A MAGNUM emission-line Capivara showcase",
    sprintf("Cube: %s", cube_path),
    sprintf("Redshift: %.5f", redshift),
    "No starlet mask, no supplied support mask, no line fitting, no continuum subtraction.",
    "Recommended default: feature_mode = 'windows'.",
    "Diagnostic modes: 'window_derivative' for sharp profile changes; 'profile' for normalized line-shape grouping.",
    "",
    "Main panel:",
    file.path(output_dir, "centa_emission_line_mode_comparison.png")
  ),
  file.path(output_dir, "README.txt")
)

print(summary)
message("Saved outputs in: ", output_dir)
