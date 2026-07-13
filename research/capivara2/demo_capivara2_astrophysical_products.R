#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

repo_dir <- normalizePath(
  Sys.getenv("CAPIVARA_REPO", unset = "/Users/rd23aag/Documents/GitHub/capivara"),
  mustWork = TRUE
)

out_dir <- Sys.getenv(
  "CAPIVARA2_ASTRO_OUT",
  unset = file.path(repo_dir, "outputs", "capivara2_astrophysical_products", "individual")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_pdf <- tolower(Sys.getenv("CAPIVARA2_WRITE_PDF", unset = "true")) %in% c("1", "true", "yes")
spatial_regularize <- tolower(Sys.getenv("CAPIVARA2_SPATIAL_REGULARIZE", unset = "true")) %in% c("1", "true", "yes")
regularize_iterations <- as.integer(Sys.getenv("CAPIVARA2_REGULARIZE_ITER", unset = "2"))
regularize_lambda <- as.numeric(Sys.getenv("CAPIVARA2_REGULARIZE_LAMBDA", unset = "0.5"))

cases <- list(
  list(
    id = "7443_12703",
    label = "7443-12703",
    seg_rds = file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "manga_7443_12703_capivara_n25.rds"),
    population_csv = file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "manga_7443_12703_capivara_n25_local_ppxf_population_by_bin.csv"),
    population_laplace_csv = file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "local_ppxf_figures", "manga_7443_12703_local_ppxf_laplace_region_values.csv"),
    emission_csv = file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "science", "manga_7443_12703_ppxf_emission_raw_by_bin.csv"),
    emission_laplace_rds = file.path(repo_dir, "outputs", "capivara_ppxf_manga_smoke", "maps", "manga_7443_12703_ppxf_emission_laplace_maps.rds")
  ),
  list(
    id = "8135_12701",
    label = "8135-12701",
    seg_rds = file.path(repo_dir, "outputs", "capivara_ppxf_manga_8135_12701", "manga_8135_12701_capivara_n25.rds"),
        population_csv = file.path(repo_dir, "outputs", "capivara_ppxf_manga_8135_12701", "manga_8135_12701_capivara_n25_local_ppxf_population_by_bin.csv"),
        population_laplace_csv = file.path(repo_dir, "outputs", "capivara_ppxf_manga_8135_12701", "local_ppxf_figures", "manga_8135_12701_local_ppxf_laplace_region_values.csv"),
        emission_csv = file.path(repo_dir, "outputs", "capivara_ppxf_manga_8135_12701", "science", "manga_8135_12701_ppxf_emission_raw_by_bin.csv"),
        emission_laplace_rds = file.path(repo_dir, "outputs", "capivara_ppxf_manga_8135_12701", "maps", "manga_8135_12701_ppxf_emission_laplace_maps.rds")
  )
)

vangogh_base <- c(
  "#0B132B", "#1C4E80", "#2F74B5", "#58A4B0", "#F2D06B",
  "#F2A541", "#C7771F", "#7A3E12", "#EFE6A4", "#244C89",
  "#1E8078", "#B84A28"
)

vangogh_discrete <- function(n) {
  if (n <= length(vangogh_base)) return(vangogh_base[seq_len(n)])
  grDevices::colorRampPalette(vangogh_base, space = "Lab")(n)
}

seq_palette <- grDevices::colorRampPalette(
  c("#0B132B", "#1C4E80", "#58A4B0", "#F2D06B", "#F2A541", "#B84A28"),
  space = "Lab"
)(256)

div_palette <- grDevices::colorRampPalette(
  c("#16345C", "#2C7FB8", "#D7E8E6", "#F2D06B", "#D46A2C", "#7F1D1D"),
  space = "Lab"
)(256)

bpt_palette <- c(
  sf = "#2F74B5",
  star_forming = "#2F74B5",
  composite = "#58A4B0",
  transition = "#58A4B0",
  agn_like = "#C7771F",
  seyfert = "#B84A28",
  liner = "#7A3E12",
  inactive = "#B8B8AA",
  unclassified = "#0B132B"
)

map_df <- function(mat) {
  out <- expand.grid(x = seq_len(nrow(mat)), y = seq_len(ncol(mat)))
  out$value <- as.vector(mat)
  out
}

theme_map <- function(base_size = 8) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = unit(9, "pt"),
      legend.key.width = unit(16, "pt"),
      plot.margin = margin(1, 1, 1, 1)
    )
}

paint_numeric <- function(cluster_map, values, value_col) {
  out <- matrix(NA_real_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
  tab <- values[is.finite(values$bin) & is.finite(values[[value_col]]), c("bin", value_col)]
  val <- stats::setNames(tab[[value_col]], tab$bin)
  idx <- is.finite(cluster_map) & as.character(cluster_map) %in% names(val)
  out[idx] <- val[as.character(cluster_map[idx])]
  out
}

paint_class <- function(cluster_map, values, value_col) {
  out <- matrix(NA_character_, nrow = nrow(cluster_map), ncol = ncol(cluster_map))
  tab <- values[is.finite(values$bin) & !is.na(values[[value_col]]) & nzchar(values[[value_col]]), c("bin", value_col)]
  val <- stats::setNames(as.character(tab[[value_col]]), tab$bin)
  idx <- is.finite(cluster_map) & as.character(cluster_map) %in% names(val)
  out[idx] <- val[as.character(cluster_map[idx])]
  out
}

regularize_numeric_map <- function(mat, support = is.finite(mat), iterations = regularize_iterations, lambda = regularize_lambda) {
  if (!spatial_regularize || iterations < 1 || lambda <= 0) return(mat)
  support <- support & is.finite(mat)
  if (!any(support)) return(mat)

  out <- mat
  for (iter in seq_len(iterations)) {
    previous <- out
    next_map <- previous
    for (i in seq_len(nrow(previous))) {
      for (j in seq_len(ncol(previous))) {
        if (!support[i, j]) next
        vals <- numeric()
        if (i > 1L && support[i - 1L, j]) vals <- c(vals, previous[i - 1L, j])
        if (i < nrow(previous) && support[i + 1L, j]) vals <- c(vals, previous[i + 1L, j])
        if (j > 1L && support[i, j - 1L]) vals <- c(vals, previous[i, j - 1L])
        if (j < ncol(previous) && support[i, j + 1L]) vals <- c(vals, previous[i, j + 1L])
        vals <- vals[is.finite(vals)]
        if (length(vals)) {
          next_map[i, j] <- (previous[i, j] + lambda * sum(vals)) / (1 + lambda * length(vals))
        }
      }
    }
    out <- next_map
  }
  out
}

load_emission_laplace_maps <- function(path) {
  if (!spatial_regularize || !file.exists(path)) return(NULL)
  maps <- readRDS(path)$maps
  if (!is.list(maps)) return(NULL)
  maps
}

finite_limits <- function(mat, probs = c(0.02, 0.98), symmetric = FALSE) {
  x <- mat[is.finite(mat)]
  if (!length(x)) return(NULL)
  if (symmetric) {
    lim <- max(abs(stats::quantile(x, probs, na.rm = TRUE)), na.rm = TRUE)
    return(c(-lim, lim))
  }
  as.numeric(stats::quantile(x, probs, na.rm = TRUE))
}

legend_label <- function(units) {
  if (identical(units, "km s^-1")) return(expression(km~s^{-1}))
  if (identical(units, "M/L_R")) return(expression(M/L[R]))
  units
}

plot_numeric <- function(mat, palette, limits = NULL, legend_title = NULL) {
  ggplot(map_df(mat), aes(x, y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_gradientn(
      colours = palette,
      limits = limits,
      oob = scales::squish,
      na.value = "white",
      name = legend_label(legend_title)
    ) +
    theme_map()
}

plot_bins <- function(cluster_map) {
  df <- map_df(cluster_map)
  df$value <- factor(df$value)
  n <- length(stats::na.omit(unique(df$value)))
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(
      values = vangogh_discrete(max(n, 1)),
      guide = "none",
      na.value = "white",
      na.translate = FALSE
    ) +
    theme_map()
}

plot_category <- function(mat, values = bpt_palette, legend_title = NULL) {
  df <- map_df(mat)
  df$value <- factor(df$value, exclude = NA)
  cols <- values[names(values) %in% levels(df$value)]
  missing <- setdiff(levels(df$value), names(cols))
  if (length(missing)) {
    cols <- c(cols, stats::setNames(vangogh_discrete(length(missing)), missing))
  }
  ggplot(df, aes(x, y, fill = value)) +
    geom_raster(na.rm = TRUE) +
    coord_equal(expand = FALSE) +
    scale_fill_manual(
      values = cols,
      na.value = "white",
      na.translate = FALSE,
      name = legend_title
    ) +
    theme_map()
}

save_plot <- function(plot, stem, width = 2.4, height = 2.4, dpi = 450) {
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  pdf_file <- NA_character_
  if (write_pdf) {
    pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
    ggsave(pdf_file, plot, width = width, height = height, bg = "white")
  }
  data.frame(file = png_file, pdf = pdf_file, stringsAsFactors = FALSE)
}

manifest <- list()

for (case in cases) {
  if (!file.exists(case$seg_rds)) {
    warning("Skipping ", case$label, "; missing segmentation RDS: ", case$seg_rds)
    next
  }

  seg <- readRDS(case$seg_rds)
  cluster_map <- seg$cluster_map

  manifest[[length(manifest) + 1L]] <- cbind(
    object = case$label,
    product = "bins",
    quantity = "capivara_bin",
    units = "categorical",
    save_plot(plot_bins(cluster_map), paste0("capivara2_astro_", case$id, "_bins_vangogh"))
  )

  pop_file <- if (file.exists(case$population_laplace_csv)) case$population_laplace_csv else case$population_csv
  if (file.exists(pop_file)) {
    pop <- utils::read.csv(pop_file)
    pop_products <- list(
      list(col = "mean_log_age", stem = "log_age_yr", units = "log yr", palette = seq_palette, symmetric = FALSE),
      list(col = "mean_metal", stem = "metallicity_mh", units = "[M/H]", palette = div_palette, symmetric = FALSE),
      list(col = "ml_r", stem = "ml_r", units = "M/L_R", palette = seq_palette, symmetric = FALSE),
      list(col = "stellar_vel_median0", stem = "stellar_velocity_kms", units = "km s^-1", palette = div_palette, symmetric = TRUE),
      list(col = "stellar_sigma", stem = "stellar_sigma_kms", units = "km s^-1", palette = seq_palette, symmetric = FALSE)
    )
    for (prod in pop_products) {
      if (!prod$col %in% names(pop)) next
      mat <- paint_numeric(cluster_map, pop, prod$col)
      mat <- regularize_numeric_map(mat)
      lim <- finite_limits(mat, symmetric = prod$symmetric)
      p <- plot_numeric(mat, prod$palette, limits = lim, legend_title = prod$units)
      manifest[[length(manifest) + 1L]] <- cbind(
        object = case$label,
        product = "stellar_population",
        quantity = prod$stem,
        units = prod$units,
        save_plot(p, paste0("capivara2_astro_", case$id, "_", prod$stem))
      )
    }
  } else {
    warning("No stellar-population table for ", case$label)
  }

  if (file.exists(case$emission_csv)) {
    em <- utils::read.csv(case$emission_csv)
    emission_laplace <- load_emission_laplace_maps(case$emission_laplace_rds)
    laplace_name <- c(
      halpha_flux = "halpha_log_flux_laplace",
      gas_vel = "gas_velocity_pattern_laplace",
      log_nii_ha = "log_nii_ha_laplace",
      log_oiii_hb = "log_oiii_hb_laplace"
    )
    emission_products <- list(
      list(col = "halpha_flux", stem = "halpha_log_flux_regularized", units = "log relative flux", palette = seq_palette, transform = "log10_positive", symmetric = FALSE),
      list(col = "gas_vel", stem = "gas_velocity_kms", units = "km s^-1", palette = div_palette, transform = "median0", symmetric = TRUE),
      list(col = "gas_sigma", stem = "gas_sigma_kms", units = "km s^-1", palette = seq_palette, transform = "identity", symmetric = FALSE),
      list(col = "log_nii_ha", stem = "log_nii_ha", units = "dex", palette = div_palette, transform = "identity", symmetric = FALSE),
      list(col = "log_oiii_hb", stem = "log_oiii_hb", units = "dex", palette = div_palette, transform = "identity", symmetric = FALSE)
    )
    for (prod in emission_products) {
      if (!prod$col %in% names(em)) next
      if (!is.null(emission_laplace) &&
          prod$col %in% names(laplace_name) &&
          laplace_name[[prod$col]] %in% names(emission_laplace)) {
        mat <- emission_laplace[[laplace_name[[prod$col]]]]
      } else {
        mat <- paint_numeric(cluster_map, em, prod$col)
        if (prod$transform == "log10_positive") {
          mat[!is.finite(mat) | mat <= 0] <- NA_real_
          mat <- log10(mat)
        }
        if (prod$transform == "median0") mat <- mat - stats::median(mat, na.rm = TRUE)
        mat <- regularize_numeric_map(mat)
      }
      lim <- finite_limits(mat, symmetric = prod$symmetric)
      p <- plot_numeric(mat, prod$palette, limits = lim, legend_title = prod$units)
      manifest[[length(manifest) + 1L]] <- cbind(
        object = case$label,
        product = "emission_line",
        quantity = prod$stem,
        units = prod$units,
        save_plot(p, paste0("capivara2_astro_", case$id, "_", prod$stem))
      )
    }

    if ("bpt_class" %in% names(em)) {
      if (!is.null(emission_laplace) && "bpt_class_laplace" %in% names(emission_laplace)) {
        mat <- emission_laplace$bpt_class_laplace
      } else {
        mat <- paint_class(cluster_map, em, "bpt_class")
      }
      p <- plot_category(mat, legend_title = "BPT")
      manifest[[length(manifest) + 1L]] <- cbind(
        object = case$label,
        product = "emission_line",
        quantity = "bpt_class",
        units = "categorical",
        save_plot(p, paste0("capivara2_astro_", case$id, "_bpt_class"))
      )
    }
  } else {
    warning("No emission-line table for ", case$label)
  }
}

manifest <- if (length(manifest)) do.call(rbind, manifest) else data.frame()
utils::write.csv(manifest, file.path(out_dir, "capivara2_astrophysical_products_manifest.csv"), row.names = FALSE)

message("Wrote astrophysical product maps to:")
message(out_dir)
message("Manifest:")
message(file.path(out_dir, "capivara2_astrophysical_products_manifest.csv"))
