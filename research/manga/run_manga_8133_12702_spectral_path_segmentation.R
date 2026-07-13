#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(capivara)
  library(spectropath)
  library(ggplot2)
  library(FITSio)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

cube_path <- Sys.getenv(
  "CAPIVARA_MANGA8133_CUBE",
  unset = "/Users/rd23aag/Documents/GitHub/capivara_experimental/manga-8133-12702-MEGACUBE.fits"
)
out_dir <- Sys.getenv(
  "CAPIVARA_MANGA8133_PATH_OUT",
  unset = file.path(getwd(), "outputs", "manga_8133_12702_agn", "spectral_window_path")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_combined <- as.integer(Sys.getenv("CAPIVARA_MANGA8133_NCOMP", unset = "40"))
n_oiii <- as.integer(Sys.getenv("CAPIVARA_MANGA8133_NCOMP_OIII", unset = as.character(n_combined)))
n_halpha <- as.integer(Sys.getenv("CAPIVARA_MANGA8133_NCOMP_HALPHA", unset = as.character(n_combined)))
n_path <- as.integer(Sys.getenv("CAPIVARA_MANGA8133_PATH_NCOMP", unset = "40"))

parse_card_value <- function(card) {
  if (nchar(card) < 10 || substr(card, 9, 9) != "=") return(NULL)
  raw <- trimws(strsplit(substr(card, 11, 80), "/", fixed = TRUE)[[1]][1])
  if (!nzchar(raw)) return(NULL)
  if (startsWith(raw, "'")) {
    end <- regexpr("'", substring(raw, 2), fixed = TRUE)[1]
    if (end > 0) return(substr(raw, 2, end))
    return(gsub("^'|'$", "", raw))
  }
  if (raw %in% c("T", "F")) return(raw == "T")
  value <- suppressWarnings(as.numeric(gsub("D", "E", raw, fixed = TRUE)))
  if (is.finite(value)) return(value)
  raw
}

read_fits_index <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  file_size <- file.info(path)$size
  offset <- 0
  hdu <- 0L
  out <- list()

  repeat {
    if (offset >= file_size) break
    cards <- character()
    repeat {
      block <- readBin(con, what = "raw", n = 2880L)
      if (!length(block)) break
      offset <- offset + 2880
      for (i in seq(1, 2880, by = 80)) {
        card <- rawToChar(block[i:(i + 79)])
        cards <- c(cards, card)
        if (startsWith(card, "END")) break
      }
      if (length(cards) && startsWith(cards[length(cards)], "END")) break
    }

    header <- list()
    for (card in cards) {
      key <- trimws(substr(card, 1, 8))
      val <- parse_card_value(card)
      if (nzchar(key) && !is.null(val)) header[[key]] <- val
    }

    hdu <- hdu + 1L
    bitpix <- as.integer(header$BITPIX %||% 0L)
    naxis <- as.integer(header$NAXIS %||% 0L)
    dims <- if (naxis > 0L) {
      as.integer(vapply(seq_len(naxis), function(i) header[[paste0("NAXIS", i)]] %||% 0L, numeric(1)))
    } else {
      integer()
    }
    bytes <- if (naxis > 0L) prod(dims) * abs(bitpix) / 8 else 0
    bytes <- bytes + as.numeric(header$PCOUNT %||% 0)
    bytes <- bytes * as.numeric(header$GCOUNT %||% 1)
    padded <- ceiling(bytes / 2880) * 2880

    out[[hdu]] <- list(
      hdu = hdu,
      extname = trimws(header$EXTNAME %||% if (hdu == 1L) "PRIMARY" else ""),
      header = header,
      cards = cards,
      data_start = offset,
      data_bytes = bytes,
      dims = dims,
      bitpix = bitpix
    )
    seek(con, where = padded, origin = "current")
    offset <- offset + padded
  }
  out
}

read_image_hdu <- function(path, index, extname) {
  extnames <- vapply(index, `[[`, "", "extname")
  rec <- index[[which(extnames == extname)[1]]]
  if (is.null(rec)) stop("Missing HDU: ", extname)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = rec$data_start, origin = "start")
  what <- switch(
    as.character(rec$bitpix),
    "-64" = "double",
    "-32" = "numeric",
    "64" = "integer",
    "32" = "integer",
    "16" = "integer",
    "8" = "raw",
    stop("Unsupported BITPIX: ", rec$bitpix)
  )
  vals <- readBin(con, what = what, n = prod(rec$dims), size = abs(rec$bitpix) / 8, endian = "big")
  if (rec$bitpix == 8) vals <- as.integer(vals)
  dim(vals) <- rec$dims
  vals
}

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else stats::median(x)
}

continuum_subtract_window <- function(cube, wave, main_range, side_blue, side_red) {
  keep <- which(wave >= min(c(main_range, side_blue, side_red)) & wave <= max(c(main_range, side_blue, side_red)))
  main <- which(wave[keep] >= main_range[1] & wave[keep] <= main_range[2])
  blue <- which(wave[keep] >= side_blue[1] & wave[keep] <= side_blue[2])
  red <- which(wave[keep] >= side_red[1] & wave[keep] <= side_red[2])
  sub <- cube[, , keep, drop = FALSE]
  w <- wave[keep]
  nx <- dim(sub)[1]
  ny <- dim(sub)[2]
  out <- array(NA_real_, dim = c(nx, ny, length(main)))
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      spec <- sub[i, j, ]
      if (sum(is.finite(spec)) < 6L) next
      xb <- c(w[blue], w[red])
      yb <- c(spec[blue], spec[red])
      ok <- is.finite(xb) & is.finite(yb)
      if (sum(ok) >= 3L) {
        fit_df <- data.frame(x = xb[ok], y = yb[ok])
        fit <- tryCatch(stats::lm(y ~ x, data = fit_df), error = function(e) NULL)
        cont <- if (is.null(fit)) rep(finite_median(yb), length(w)) else as.numeric(stats::predict(fit, newdata = data.frame(x = w)))
      } else {
        cont <- rep(finite_median(spec), length(w))
      }
      out[i, j, ] <- spec[main] - cont[main]
    }
  }
  list(cube = out, wave = w[main])
}

safe_log10 <- function(x) {
  x[x <= 0 | !is.finite(x)] <- NA_real_
  log10(x)
}

integrate_window <- function(cube, wave, range = range(wave), positive_only = TRUE) {
  idx <- which(wave >= range[1] & wave <= range[2])
  dx <- median(diff(wave[idx]))
  apply(cube[, , idx, drop = FALSE], c(1, 2), function(v) {
    if (positive_only) v <- pmax(v, 0)
    if (!any(is.finite(v))) return(NA_real_)
    sum(v, na.rm = TRUE) * dx
  })
}

robust_support <- function(line_flux_maps) {
  total <- Reduce(`+`, lapply(line_flux_maps, function(x) {
    y <- x
    y[!is.finite(y) | y < 0] <- 0
    y
  }))
  q <- stats::quantile(total[is.finite(total) & total > 0], 0.08, na.rm = TRUE)
  is.finite(total) & total > q
}

build_starlet_support_from_cube <- function(cube, mask2d = NULL) {
  cube[!is.finite(cube)] <- 0
  starlet <- capivara::build_starlet_mask(
    input = list(imDat = cube, hdr = NULL, axDat = NULL),
    starlet_J = 5,
    starlet_scales = 2:5,
    include_coarse = FALSE,
    denoise_k = 0,
    mode = "soft",
    positive_only = TRUE
  )
  support <- starlet$mask
  if (!is.null(mask2d)) {
    support <- support & (is.na(mask2d) | mask2d == 0 | is.finite(mask2d))
  }
  if (!any(support, na.rm = TRUE)) {
    stop("The starlet support is empty; inspect the spectral windows or starlet parameters.")
  }
  list(
    support = support,
    starlet = starlet
  )
}

run_seg <- function(cube, mask, ncomp, spatial_weight = 0.20, label = "") {
  message("Running Capivara spectral-window segmentation: ", label)
  cube[!is.finite(cube)] <- 0
  capivara::segment_large(
    input = list(imDat = cube, hdr = NULL, axDat = NULL),
    Ncomp = ncomp,
    scale_fn = capivara::median_scale,
    feature_scale = "robust_col",
    spatial_weight = spatial_weight,
    mask = mask,
    valid_mode = "finite",
    knn_k = 55,
    auto_k = TRUE,
    max_k = 180,
    verbose = TRUE
  )
}

path_feature_matrix <- function(line_cube, wave, rest_wave, mask, max_abs_velocity = 900) {
  vel <- 299792.458 * (wave / rest_wave - 1)
  keep <- abs(vel) <= max_abs_velocity
  vel <- vel[keep]
  cube <- line_cube[, , keep, drop = FALSE]
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  rows <- vector("list", nx * ny)
  n <- 0L
  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      if (!isTRUE(mask[i, j])) next
      flux <- as.numeric(cube[i, j, ])
      if (sum(is.finite(flux)) < 8L) next
      flux[!is.finite(flux)] <- 0
      if (max(abs(flux), na.rm = TRUE) <= 0) next
      # Shape-normalize for path features, but keep sign and ordering.
      flux <- flux / max(abs(flux), na.rm = TRUE)
      pf <- tryCatch(
        spectropath::path_features(cbind(vel, flux), depth = 4, normalize = TRUE, notation = "paper"),
        error = function(e) NULL
      )
      if (is.null(pf)) next
      n <- n + 1L
      rows[[n]] <- cbind(data.frame(x = i, y = j), pf)
    }
  }
  if (!n) return(data.frame())
  do.call(rbind, rows[seq_len(n)])
}

weighted_velocity <- function(flux, vel, positive_only = TRUE) {
  f <- as.numeric(flux)
  if (positive_only) f <- pmax(f, 0)
  ok <- is.finite(f) & is.finite(vel)
  if (sum(ok) < 3L || sum(abs(f[ok])) <= 0) return(NA_real_)
  sum(vel[ok] * f[ok], na.rm = TRUE) / sum(f[ok], na.rm = TRUE)
}

weighted_sigma <- function(flux, vel, mu, positive_only = TRUE) {
  f <- as.numeric(flux)
  if (positive_only) f <- pmax(f, 0)
  ok <- is.finite(f) & is.finite(vel) & is.finite(mu)
  if (sum(ok) < 3L || sum(abs(f[ok])) <= 0) return(NA_real_)
  sqrt(sum((vel[ok] - mu)^2 * f[ok], na.rm = TRUE) / sum(f[ok], na.rm = TRUE))
}

baseline_subtract <- function(v, y, n_edge = 2L) {
  edge <- c(seq_len(min(n_edge, length(y))), seq.int(max(1L, length(y) - n_edge + 1L), length(y)))
  base <- stats::median(y[edge], na.rm = TRUE)
  y - base
}

nonparam_velocity <- function(v, y) {
  line <- baseline_subtract(v, y)
  pos <- pmax(line, 0)
  flux <- sum(pos, na.rm = TRUE)
  if (!is.finite(flux) || flux <= 0) {
    return(data.frame(
      np_flux = NA_real_,
      np_peak_v = NA_real_,
      np_centroid = NA_real_,
      np_sigma = NA_real_,
      np_w80 = NA_real_
    ))
  }

  centroid <- sum(v * pos, na.rm = TRUE) / flux
  sigma <- sqrt(sum((v - centroid)^2 * pos, na.rm = TRUE) / flux)
  ord <- order(v)
  cum <- cumsum(pos[ord]) / flux
  qv <- stats::approx(cum, v[ord], xout = c(0.1, 0.5, 0.9), ties = "ordered", rule = 2)$y

  data.frame(
    np_flux = flux,
    np_peak_v = v[which.max(line)],
    np_centroid = centroid,
    np_sigma = sigma,
    np_w80 = qv[3] - qv[1]
  )
}

segment_kinematics <- function(seg_map,
                               line_cube,
                               wave,
                               rest_wave,
                               mask,
                               label,
                               max_abs_velocity = 600) {
  vel <- 299792.458 * (wave / rest_wave - 1)
  keep <- abs(vel) <= max_abs_velocity
  vel <- vel[keep]
  line_cube <- line_cube[, , keep, drop = FALSE]
  ids <- sort(unique(as.integer(seg_map[is.finite(seg_map)])))
  rows <- lapply(ids, function(id) {
    pix <- which(seg_map == id & mask, arr.ind = TRUE)
    if (!nrow(pix)) return(NULL)
    spectra <- t(apply(pix, 1, function(rc) line_cube[rc[1], rc[2], ]))
    stack <- colSums(spectra, na.rm = TRUE)
    mu <- weighted_velocity(stack, vel)
    sig <- weighted_sigma(stack, vel, mu)
    pf <- tryCatch(
      spectropath::path_features(cbind(vel, stack / max(abs(stack), na.rm = TRUE)), depth = 4, normalize = TRUE, notation = "paper"),
      error = function(e) data.frame(p2 = NA, p_pm = NA, p3u = NA, p3F = NA, p4F = NA, p4T = NA)
    )
    cbind(data.frame(segment = id, n_pix = nrow(pix), velocity = mu, sigma = sig, line = label), pf)
  })
  out <- do.call(rbind, rows)
  if (!is.null(out) && nrow(out)) {
    v0 <- stats::median(out$velocity[is.finite(out$velocity)], na.rm = TRUE)
    out$velocity_centered <- out$velocity - v0
  }
  out
}

nearest_impute_feature_cube <- function(feature_cube, mask) {
  dims <- dim(feature_cube)
  mat <- matrix(feature_cube, nrow = dims[1] * dims[2], ncol = dims[3])
  xy <- expand.grid(x = seq_len(dims[1]), y = seq_len(dims[2]))
  support <- as.vector(mask)
  complete <- support & apply(mat, 1, function(z) all(is.finite(z)))
  missing <- support & !complete

  if (!any(complete)) {
    stop("No complete path-feature spaxels were available for clustering.")
  }

  if (any(missing)) {
    source_rows <- which(complete)
    source_xy <- as.matrix(xy[source_rows, c("x", "y")])
    for (row in which(missing)) {
      d2 <- (source_xy[, 1] - xy$x[row])^2 + (source_xy[, 2] - xy$y[row])^2
      mat[row, ] <- mat[source_rows[which.min(d2)], ]
    }
  }

  feature_cube[] <- mat
  feature_cube
}

build_path_feature_cube <- function(line_cube,
                                    wave,
                                    rest_wave,
                                    mask,
                                    max_abs_velocity = 600,
                                    feature_names = c("p2", "p3u", "p3F", "p4F", "p4T", "p_pm")) {
  vel <- 299792.458 * (wave / rest_wave - 1)
  keep <- abs(vel) <= max_abs_velocity
  vel <- vel[keep]
  cube <- line_cube[, , keep, drop = FALSE]
  nx <- dim(cube)[1]
  ny <- dim(cube)[2]
  feature_cube <- array(NA_real_, dim = c(nx, ny, length(feature_names)), dimnames = list(NULL, NULL, feature_names))
  rows <- vector("list", sum(mask, na.rm = TRUE))
  n <- 0L

  for (i in seq_len(nx)) {
    for (j in seq_len(ny)) {
      if (!isTRUE(mask[i, j])) next
      flux <- as.numeric(cube[i, j, ])
      if (sum(is.finite(flux)) < 8L) next
      flux[!is.finite(flux)] <- 0
      line <- baseline_subtract(vel, flux)
      amp <- max(abs(line), na.rm = TRUE)
      if (!is.finite(amp) || amp <= 0) next
      profile <- line / amp
      pf <- tryCatch(
        spectropath::path_features(cbind(vel, profile), depth = 4, normalize = TRUE, notation = "paper"),
        error = function(e) NULL
      )
      if (is.null(pf)) next
      vals <- as.numeric(pf[1, feature_names, drop = TRUE])
      feature_cube[i, j, ] <- vals

      np <- nonparam_velocity(vel, line)
      n <- n + 1L
      rows[[n]] <- data.frame(
        x = i,
        y = j,
        line_flux = np$np_flux,
        centroid_kms = np$np_centroid,
        sigma_kms = np$np_sigma,
        w80_kms = np$np_w80,
        pf[, feature_names, drop = FALSE],
        check.names = FALSE
      )
    }
  }

  feature_cube <- nearest_impute_feature_cube(feature_cube, mask)
  table <- if (n) do.call(rbind, rows[seq_len(n)]) else data.frame()

  list(
    feature_cube = feature_cube,
    table = table,
    velocity = vel,
    line_cube = cube,
    features = feature_names
  )
}

run_path_feature_seg <- function(feature_cube, mask, ncomp, label = "") {
  message("Running Capivara path-feature segmentation: ", label)
  capivara::segment_large(
    input = list(imDat = feature_cube, hdr = NULL, axDat = NULL),
    Ncomp = ncomp,
    feature_scale = "robust_col",
    spatial_weight = as.numeric(Sys.getenv("CAPIVARA_MANGA8133_PATH_SPATIAL_WEIGHT", unset = "0.10")),
    mask = mask,
    valid_mode = "finite",
    knn_k = as.integer(Sys.getenv("CAPIVARA_MANGA8133_PATH_KNN", unset = "55")),
    auto_k = TRUE,
    max_k = 180,
    verbose = TRUE
  )
}

map_segment_values <- function(seg_map, summary, value_col) {
  out <- matrix(NA_real_, nrow = nrow(seg_map), ncol = ncol(seg_map))
  vals <- setNames(summary[[value_col]], summary$segment)
  ids <- as.character(as.integer(seg_map))
  out[] <- vals[ids]
  out
}

plot_map <- function(map, path, title = NULL, diverging = FALSE, categorical = FALSE) {
  df <- expand.grid(x = seq_len(nrow(map)), y = seq_len(ncol(map)))
  df$value <- as.vector(map)
  van_gogh <- c("#172033", "#223E6B", "#2F74B5", "#2897A8", "#58B7A6", "#E9D8A6", "#F2A541", "#D9792B", "#B84A28", "#7A2E1F")
  div_cols <- c("#223E6B", "#2F74B5", "#A9C9C8", "#F8F0D0", "#F2A541", "#B84A28", "#7A2E1F")
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.02),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position = "right"
    )
  if (categorical) {
    cols <- grDevices::colorRampPalette(van_gogh)(max(as.integer(df$value), na.rm = TRUE))
    p <- p + scale_fill_gradientn(colours = cols, na.value = "white", guide = "none")
  } else if (diverging) {
    lim <- stats::quantile(abs(df$value), 0.98, na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- max(abs(df$value), na.rm = TRUE)
    p <- p + scale_fill_gradient2(low = div_cols[1], mid = div_cols[4], high = div_cols[7], midpoint = 0, limits = c(-lim, lim), oob = scales::squish, na.value = "white")
  } else {
    qs <- stats::quantile(df$value, c(0.02, 0.98), na.rm = TRUE)
    p <- p + scale_fill_gradientn(colours = van_gogh, limits = qs, oob = scales::squish, na.value = "white")
  }
  ggsave(path, p, width = 5.4, height = 4.8, dpi = 320, bg = "white")
  invisible(p)
}

idx <- read_fits_index(cube_path)
fitspec <- read_image_hdu(cube_path, idx, "FITSPEC")
fitcont <- read_image_hdu(cube_path, idx, "FITCONT")
flux_m <- read_image_hdu(cube_path, idx, "FLUX_M")
mask2d <- read_image_hdu(cube_path, idx, "MASK2D")

wave <- 4401 + seq(0, dim(fitspec)[3] - 1)
resid <- fitspec - fitcont

line_names <- c("hb", "o3_4959", "o3_5007", "he1_5876", "o1_6300", "n2_6548", "ha", "n2_6583", "s2_6716", "s2_6731")
dimnames(flux_m) <- list(NULL, NULL, line_names)

o3_flux <- flux_m[, , "o3_5007"]
hb_flux <- flux_m[, , "hb"]
ha_flux <- flux_m[, , "ha"]
nii_flux <- flux_m[, , "n2_6583"]
sii_flux <- flux_m[, , "s2_6716"] + flux_m[, , "s2_6731"]

o3hb <- continuum_subtract_window(resid, wave, main_range = c(4825, 5035), side_blue = c(4750, 4805), side_red = c(5065, 5120))
ha_nii_sii <- continuum_subtract_window(resid, wave, main_range = c(6500, 6765), side_blue = c(6420, 6480), side_red = c(6785, 6860))
oi_sii <- continuum_subtract_window(resid, wave, main_range = c(6265, 6765), side_blue = c(6200, 6245), side_red = c(6785, 6860))
o3hb$cube[!is.finite(o3hb$cube)] <- 0
ha_nii_sii$cube[!is.finite(ha_nii_sii$cube)] <- 0
oi_sii$cube[!is.finite(oi_sii$cube)] <- 0

combined_cube <- abind::abind(o3hb$cube, ha_nii_sii$cube, along = 3)
support_info <- build_starlet_support_from_cube(fitspec, mask2d = mask2d)
support <- support_info$support
message("Full-spectrum starlet support pixels: ", sum(support, na.rm = TRUE))

combined_seg <- run_seg(combined_cube, support, ncomp = n_combined, spatial_weight = 0.18, label = "OIII/Hbeta + Halpha/NII/SII")
oiii_seg <- run_seg(o3hb$cube, support, ncomp = n_oiii, spatial_weight = 0.15, label = "OIII/Hbeta spectral window")
halpha_seg <- run_seg(ha_nii_sii$cube, support, ncomp = n_halpha, spatial_weight = 0.15, label = "Halpha/NII/SII spectral window")

compute_spaxel_paths <- tolower(Sys.getenv("CAPIVARA_MANGA8133_SPAXEL_PATHS", unset = "false")) %in% c("1", "true", "yes")
if (isTRUE(compute_spaxel_paths)) {
  o3_path <- path_feature_matrix(o3hb$cube, o3hb$wave, rest_wave = 5008.239, mask = support)
  ha_path <- path_feature_matrix(ha_nii_sii$cube, ha_nii_sii$wave, rest_wave = 6564.61, mask = support)
} else {
  o3_path <- data.frame()
  ha_path <- data.frame()
}
write.csv(o3_path, file.path(out_dir, "oiii5007_spaxel_path_features.csv"), row.names = FALSE)
write.csv(ha_path, file.path(out_dir, "halpha_spaxel_path_features.csv"), row.names = FALSE)

o3_kin <- segment_kinematics(
  combined_seg$cluster_map,
  o3hb$cube,
  o3hb$wave,
  rest_wave = 5008.239,
  support,
  "oiii5007",
  max_abs_velocity = 700
)
ha_kin <- segment_kinematics(
  combined_seg$cluster_map,
  ha_nii_sii$cube,
  ha_nii_sii$wave,
  rest_wave = 6564.61,
  support,
  "halpha",
  max_abs_velocity = 500
)
write.csv(o3_kin, file.path(out_dir, "combined_segments_oiii5007_path_kinematics.csv"), row.names = FALSE)
write.csv(ha_kin, file.path(out_dir, "combined_segments_halpha_path_kinematics.csv"), row.names = FALSE)

o3_path_style <- build_path_feature_cube(
  o3hb$cube,
  o3hb$wave,
  rest_wave = 5008.239,
  mask = support,
  max_abs_velocity = 700
)
ha_path_style <- build_path_feature_cube(
  ha_nii_sii$cube,
  ha_nii_sii$wave,
  rest_wave = 6564.61,
  mask = support,
  max_abs_velocity = 500
)
write.csv(o3_path_style$table, file.path(out_dir, "path_signature_oiii5007_spaxel_features.csv"), row.names = FALSE)
write.csv(ha_path_style$table, file.path(out_dir, "path_signature_halpha_spaxel_features.csv"), row.names = FALSE)

o3_path_seg <- run_path_feature_seg(o3_path_style$feature_cube, support, n_path, label = "[OIII] 5007 path features")
ha_path_seg <- run_path_feature_seg(ha_path_style$feature_cube, support, n_path, label = "Halpha path features")
o3_path_kin <- segment_kinematics(
  o3_path_seg$cluster_map,
  o3hb$cube,
  o3hb$wave,
  rest_wave = 5008.239,
  support,
  "oiii5007_path_segments",
  max_abs_velocity = 700
)
ha_path_kin <- segment_kinematics(
  ha_path_seg$cluster_map,
  ha_nii_sii$cube,
  ha_nii_sii$wave,
  rest_wave = 6564.61,
  support,
  "halpha_path_segments",
  max_abs_velocity = 500
)
write.csv(o3_path_kin, file.path(out_dir, "path_signature_oiii5007_segment_kinematics.csv"), row.names = FALSE)
write.csv(ha_path_kin, file.path(out_dir, "path_signature_halpha_segment_kinematics.csv"), row.names = FALSE)

tab <- expand.grid(x = seq_len(dim(resid)[1]), y = seq_len(dim(resid)[2]))
tab$support <- as.vector(support)
tab$combined_segment <- as.vector(combined_seg$cluster_map)
tab$oiii_segment <- as.vector(oiii_seg$cluster_map)
tab$halpha_segment <- as.vector(halpha_seg$cluster_map)
tab$path_oiii_segment <- as.vector(o3_path_seg$cluster_map)
tab$path_halpha_segment <- as.vector(ha_path_seg$cluster_map)
tab$log_oiii_hbeta <- as.vector(safe_log10(o3_flux) - safe_log10(hb_flux))
tab$log_nii_halpha <- as.vector(safe_log10(nii_flux) - safe_log10(ha_flux))
tab$log_sii_halpha <- as.vector(safe_log10(sii_flux) - safe_log10(ha_flux))
write.csv(tab, file.path(out_dir, "spectral_window_segmentation_spaxel_table.csv"), row.names = FALSE)

saveRDS(
  list(
    cube_path = cube_path,
    wave = wave,
    support = support,
    support_info = support_info,
    o3hb_wave = o3hb$wave,
    ha_nii_sii_wave = ha_nii_sii$wave,
    combined_seg = combined_seg,
    oiii_seg = oiii_seg,
    halpha_seg = halpha_seg,
    o3_path_style = o3_path_style,
    ha_path_style = ha_path_style,
    o3_path_seg = o3_path_seg,
    ha_path_seg = ha_path_seg,
    o3_kin = o3_kin,
    ha_kin = ha_kin,
    o3_path_kin = o3_path_kin,
    ha_path_kin = ha_path_kin
  ),
  file.path(out_dir, "spectral_window_path_segmentation_result.rds")
)

plot_map(support + 0, file.path(out_dir, "support_mask.png"), "spectral-window support")
plot_map(support_info$starlet$collapsed, file.path(out_dir, "starlet_collapsed_full_spectra.png"), "full-spectrum white image")
plot_map(support_info$starlet$reconstruction, file.path(out_dir, "starlet_reconstruction_full_spectra.png"), "starlet full-spectrum support image")
plot_map(combined_seg$cluster_map, file.path(out_dir, sprintf("segmentation_combined_lines_n%d.png", n_combined)), "Capivara spectral-window segmentation", categorical = TRUE)
plot_map(oiii_seg$cluster_map, file.path(out_dir, sprintf("segmentation_oiii_hbeta_n%d.png", n_oiii)), "Capivara OIII/Hbeta window", categorical = TRUE)
plot_map(halpha_seg$cluster_map, file.path(out_dir, sprintf("segmentation_halpha_nii_sii_n%d.png", n_halpha)), "Capivara Halpha/NII/SII window", categorical = TRUE)
plot_map(map_segment_values(combined_seg$cluster_map, o3_kin, "velocity_centered"), file.path(out_dir, "path_kinematics_oiii_velocity.png"), "[OIII] path velocity by Capivara segment", diverging = TRUE)
plot_map(map_segment_values(combined_seg$cluster_map, ha_kin, "velocity_centered"), file.path(out_dir, "path_kinematics_halpha_velocity.png"), "Halpha path velocity by Capivara segment", diverging = TRUE)
plot_map(map_segment_values(combined_seg$cluster_map, o3_kin, "sigma"), file.path(out_dir, "path_kinematics_oiii_sigma.png"), "[OIII] path sigma by Capivara segment")
plot_map(o3_path_seg$cluster_map, file.path(out_dir, sprintf("path_signature_oiii5007_segments_n%d.png", n_path)), "[OIII] path-feature segmentation", categorical = TRUE)
plot_map(ha_path_seg$cluster_map, file.path(out_dir, sprintf("path_signature_halpha_segments_n%d.png", n_path)), "Halpha path-feature segmentation", categorical = TRUE)
plot_map(map_segment_values(o3_path_seg$cluster_map, o3_path_kin, "velocity_centered"), file.path(out_dir, "path_signature_oiii5007_velocity.png"), "[OIII] velocity from path segments", diverging = TRUE)
plot_map(map_segment_values(ha_path_seg$cluster_map, ha_path_kin, "velocity_centered"), file.path(out_dir, "path_signature_halpha_velocity.png"), "Halpha velocity from path segments", diverging = TRUE)
plot_map(map_segment_values(o3_path_seg$cluster_map, o3_path_kin, "sigma"), file.path(out_dir, "path_signature_oiii5007_sigma.png"), "[OIII] sigma from path segments")
plot_map(safe_log10(o3_flux) - safe_log10(hb_flux), file.path(out_dir, "map_log_oiii_hbeta_fitflux.png"), "fit flux log [OIII]/Hbeta")
plot_map(safe_log10(nii_flux) - safe_log10(ha_flux), file.path(out_dir, "map_log_nii_halpha_fitflux.png"), "fit flux log [NII]/Halpha")

stack <- array(NA_real_, dim = c(dim(resid)[1], dim(resid)[2], 14))
stack[, , 1] <- support + 0
stack[, , 2] <- combined_seg$cluster_map
stack[, , 3] <- oiii_seg$cluster_map
stack[, , 4] <- halpha_seg$cluster_map
stack[, , 5] <- map_segment_values(combined_seg$cluster_map, o3_kin, "velocity_centered")
stack[, , 6] <- map_segment_values(combined_seg$cluster_map, ha_kin, "velocity_centered")
stack[, , 7] <- map_segment_values(combined_seg$cluster_map, o3_kin, "sigma")
stack[, , 8] <- safe_log10(o3_flux) - safe_log10(hb_flux)
stack[, , 9] <- safe_log10(nii_flux) - safe_log10(ha_flux)
stack[, , 10] <- o3_path_seg$cluster_map
stack[, , 11] <- ha_path_seg$cluster_map
stack[, , 12] <- map_segment_values(o3_path_seg$cluster_map, o3_path_kin, "velocity_centered")
stack[, , 13] <- map_segment_values(ha_path_seg$cluster_map, ha_path_kin, "velocity_centered")
stack[, , 14] <- map_segment_values(o3_path_seg$cluster_map, o3_path_kin, "sigma")
fits_path <- file.path(out_dir, "manga_8133_12702_spectral_path_maps.fits")
if (file.exists(fits_path)) unlink(fits_path)
try(FITSio::writeFITSim(stack, file = fits_path, type = "double"), silent = TRUE)
write.csv(
  data.frame(
    channel = seq_len(dim(stack)[3]),
    name = c(
      "support", "combined_segment", "oiii_segment", "halpha_segment",
      "oiii_spectral_segment_velocity_centered", "halpha_spectral_segment_velocity_centered",
      "oiii_spectral_segment_sigma", "log_oiii_hbeta", "log_nii_halpha",
      "path_oiii_segment", "path_halpha_segment",
      "path_oiii_velocity_centered", "path_halpha_velocity_centered",
      "path_oiii_sigma"
    )
  ),
  file.path(out_dir, "spectral_path_fits_channels.csv"),
  row.names = FALSE
)

message("Wrote spectral-window/path outputs to: ", out_dir)
