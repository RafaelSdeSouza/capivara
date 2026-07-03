#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(capivara)
  library(ggplot2)
  library(FITSio)
})

cube_path <- Sys.getenv(
  "CAPIVARA_MANGA8133_CUBE",
  unset = "/Users/rd23aag/Documents/GitHub/capivara_experimental/manga-8133-12702-MEGACUBE.fits"
)
out_dir <- Sys.getenv(
  "CAPIVARA_MANGA8133_OUT",
  unset = file.path(getwd(), "outputs", "manga_8133_12702_agn")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
  val <- suppressWarnings(as.numeric(gsub("D", "E", raw, fixed = TRUE)))
  if (is.finite(val)) return(val)
  raw
}

read_fits_index <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  size <- file.info(path)$size
  offset <- 0
  hdu <- 0L
  out <- list()

  repeat {
    if (offset >= size) break
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

`%||%` <- function(x, y) if (is.null(x)) y else x

read_image_hdu <- function(path, index, extname) {
  rec <- index[[which(vapply(index, `[[`, "", "extname") == extname)[1]]]
  if (is.na(rec$hdu)) stop("Missing HDU: ", extname)
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
  n <- prod(rec$dims)
  vals <- readBin(con, what = what, n = n, size = abs(rec$bitpix) / 8, endian = "big")
  if (rec$bitpix == 8) vals <- as.integer(vals)
  dim(vals) <- rec$dims
  vals
}

clean <- function(x) {
  x[!is.finite(x)] <- NA_real_
  x
}

line_names <- c(
  "hb", "o3_4959", "o3_5007", "he1_5876", "o1_6300",
  "n2_6548", "ha", "n2_6583", "s2_6716", "s2_6731"
)
param_names <- as.vector(rbind(
  paste0(line_names, "_amp"),
  paste0(line_names, "_vel"),
  paste0(line_names, "_sigma")
))

idx <- read_fits_index(cube_path)
meta <- do.call(rbind, lapply(idx, function(z) {
  data.frame(
    hdu = z$hdu,
    extname = z$extname,
    dims = paste(z$dims, collapse = "x"),
    bitpix = z$bitpix,
    data_start = z$data_start,
    data_bytes = z$data_bytes
  )
}))
write.csv(meta, file.path(out_dir, "fits_hdu_index.csv"), row.names = FALSE)

primary <- idx[[1]]$header
flux_hdr <- idx[[which(vapply(idx, `[[`, "", "extname") == "FLUX")[1]]]$header
config_cards <- idx[[which(vapply(idx, `[[`, "", "extname") == "FITCONFIG")[1]]]$cards
writeLines(idx[[1]]$cards, file.path(out_dir, "primary_header_cards.txt"))
writeLines(config_cards, file.path(out_dir, "fitconfig_header_cards.txt"))

flux_m <- read_image_hdu(cube_path, idx, "FLUX_M")
solution <- read_image_hdu(cube_path, idx, "SOLUTION")
disp <- read_image_hdu(cube_path, idx, "DISP")
mask2d <- read_image_hdu(cube_path, idx, "MASK2D")
gimg <- read_image_hdu(cube_path, idx, "GIMG")
fitspec <- read_image_hdu(cube_path, idx, "FITSPEC")

dimnames(flux_m) <- list(NULL, NULL, line_names)
dimnames(solution) <- list(NULL, NULL, c(param_names, "extra"))
dimnames(disp) <- list(NULL, NULL, param_names)

nx <- dim(flux_m)[1]
ny <- dim(flux_m)[2]

get_line <- function(arr, nm) clean(arr[, , nm])
get_param <- function(line, par) clean(solution[, , paste0(line, "_", par)])

hb <- get_line(flux_m, "hb")
o3 <- get_line(flux_m, "o3_5007")
ha <- get_line(flux_m, "ha")
nii <- get_line(flux_m, "n2_6583")
sii <- get_line(flux_m, "s2_6716") + get_line(flux_m, "s2_6731")
oi <- get_line(flux_m, "o1_6300")

o3_v <- get_param("o3_5007", "vel")
ha_v <- get_param("ha", "vel")
o3_s <- get_param("o3_5007", "sigma")
ha_s <- get_param("ha", "sigma")

safe_log10 <- function(x) {
  x[x <= 0 | !is.finite(x)] <- NA_real_
  log10(x)
}
log_ratio <- function(num, den) safe_log10(num) - safe_log10(den)

log_o3_hb <- log_ratio(o3, hb)
log_nii_ha <- log_ratio(nii, ha)
log_sii_ha <- log_ratio(sii, ha)
log_oi_ha <- log_ratio(oi, ha)

kewley01 <- function(x) 0.61 / (x - 0.47) + 1.19
kauff03 <- function(x) 0.61 / (x - 0.05) + 1.30
bpt_class <- matrix(NA_character_, nx, ny)
valid_bpt <- is.finite(log_o3_hb) & is.finite(log_nii_ha)
bpt_class[valid_bpt & log_o3_hb < kauff03(log_nii_ha)] <- "SF"
bpt_class[valid_bpt & log_o3_hb >= kauff03(log_nii_ha) & log_o3_hb < kewley01(log_nii_ha)] <- "Composite"
bpt_class[valid_bpt & log_o3_hb >= kewley01(log_nii_ha)] <- "AGN/LINER"

finite_line <- is.finite(o3) | is.finite(ha) | is.finite(hb)
positive_line <- (is.finite(o3) & o3 > 0) | (is.finite(ha) & ha > 0)
mask <- finite_line & positive_line & (is.na(mask2d) | mask2d == 0)
if (sum(mask, na.rm = TRUE) < 100) {
  mask <- finite_line & positive_line
}

center_map <- function(x, m = mask) {
  x - stats::median(x[m & is.finite(x)], na.rm = TRUE)
}

feature_cube <- array(NA_real_, dim = c(nx, ny, 12))
feature_names <- c(
  "log_o3_flux", "log_hb_flux", "log_ha_flux", "log_nii_flux", "log_sii_flux",
  "log_o3_hb", "log_nii_ha", "log_sii_ha",
  "o3_velocity_centered", "ha_velocity_centered", "o3_sigma", "ha_sigma"
)
feature_cube[, , 1] <- safe_log10(o3)
feature_cube[, , 2] <- safe_log10(hb)
feature_cube[, , 3] <- safe_log10(ha)
feature_cube[, , 4] <- safe_log10(nii)
feature_cube[, , 5] <- safe_log10(sii)
feature_cube[, , 6] <- log_o3_hb
feature_cube[, , 7] <- log_nii_ha
feature_cube[, , 8] <- log_sii_ha
feature_cube[, , 9] <- center_map(o3_v)
feature_cube[, , 10] <- center_map(ha_v)
feature_cube[, , 11] <- o3_s
feature_cube[, , 12] <- ha_s
dimnames(feature_cube) <- list(NULL, NULL, feature_names)

run_segment_large <- function(cube, ncomp, mask, spatial_weight = 0.35) {
  capivara::segment_large(
    input = list(imDat = cube, hdr = NULL, axDat = NULL),
    Ncomp = ncomp,
    feature_scale = "robust_col",
    spatial_weight = spatial_weight,
    mask = mask,
    valid_mode = "finite",
    knn_k = 50,
    auto_k = TRUE,
    max_k = 160,
    verbose = TRUE
  )
}

message("Running AGN diagnostic segmentation...")
agn_seg <- run_segment_large(feature_cube, ncomp = 30, mask = mask, spatial_weight = 0.35)

fit_wave <- 4401 + seq(0, dim(fitspec)[3] - 1)
oiii_idx <- which(fit_wave >= 4900 & fit_wave <= 5055)
oiii_cube <- fitspec[, , oiii_idx, drop = FALSE]
oiii_cube[!is.finite(oiii_cube)] <- NA_real_
oiii_mask <- mask & is.finite(o3) & o3 > stats::quantile(o3[o3 > 0 & is.finite(o3)], 0.05, na.rm = TRUE)

message("Running OIII-window segmentation...")
oiii_seg <- run_segment_large(oiii_cube, ncomp = 25, mask = oiii_mask, spatial_weight = 0.25)

saveRDS(
  list(
    cube_path = cube_path,
    meta = meta,
    primary_header = primary,
    flux_header = flux_hdr,
    line_names = line_names,
    param_names = param_names,
    feature_names = feature_names,
    mask = mask,
    oiii_mask = oiii_mask,
    agn_seg = agn_seg,
    oiii_seg = oiii_seg
  ),
  file.path(out_dir, "manga_8133_12702_agn_segmentation_result.rds")
)

maps <- list(
  gimg = clean(gimg),
  o3_flux = o3,
  hb_flux = hb,
  ha_flux = ha,
  nii_flux = nii,
  sii_flux = sii,
  log_o3_hb = log_o3_hb,
  log_nii_ha = log_nii_ha,
  log_sii_ha = log_sii_ha,
  log_oi_ha = log_oi_ha,
  o3_velocity = o3_v,
  ha_velocity = ha_v,
  o3_velocity_centered = center_map(o3_v),
  ha_velocity_centered = center_map(ha_v),
  o3_sigma = o3_s,
  ha_sigma = ha_s,
  agn_segment = agn_seg$cluster_map,
  oiii_segment = oiii_seg$cluster_map
)

tab <- expand.grid(x = seq_len(nx), y = seq_len(ny))
for (nm in names(maps)) tab[[nm]] <- as.vector(maps[[nm]])
tab$bpt_class <- as.vector(bpt_class)
tab$mask <- as.vector(mask)
tab$oiii_mask <- as.vector(oiii_mask)
write.csv(tab, file.path(out_dir, "manga_8133_12702_agn_maps.csv"), row.names = FALSE)

write.csv(
  data.frame(
    item = c("plateifu", "mangaid", "objra", "objdec", "redshift_header", "flux_wave_start", "flux_wave_step", "fit_wave_start", "fit_wave_step"),
    value = c(
      primary$PLATEIFU %||% "8133-12702",
      primary$MANGAID %||% NA,
      primary$OBJRA %||% NA,
      primary$OBJDEC %||% NA,
      "not present in MEGACUBE header",
      flux_hdr$CRVAL3 %||% NA,
      flux_hdr$CD3_3 %||% NA,
      4401,
      1
    )
  ),
  file.path(out_dir, "manga_8133_12702_metadata_summary.csv"),
  row.names = FALSE
)

van_gogh <- c(
  "#172033", "#223E6B", "#2F74B5", "#2897A8", "#58B7A6",
  "#E9D8A6", "#F2A541", "#D9792B", "#B84A28", "#7A2E1F"
)
div_cols <- c("#223E6B", "#2F74B5", "#A9C9C8", "#F2E1A6", "#F2A541", "#B84A28", "#7A2E1F")
seg_cols <- grDevices::colorRampPalette(c("#172033", "#2F74B5", "#58B7A6", "#E9D8A6", "#F2A541", "#B84A28", "#7A2E1F"))

plot_map <- function(map, path, title = NULL, palette = van_gogh, midpoint = NULL, categorical = FALSE) {
  df <- expand.grid(x = seq_len(nrow(map)), y = seq_len(ncol(map)))
  df$value <- as.vector(map)
  p <- ggplot(df, aes(x, y, fill = value)) +
    geom_raster() +
    coord_fixed(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.02, size = 12),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.position = "right"
    )
  if (categorical) {
    vals <- sort(unique(as.integer(map[is.finite(map)])))
    cols <- seg_cols(max(vals, na.rm = TRUE))
    p <- p + scale_fill_gradientn(colours = cols, na.value = "white", guide = "none")
  } else if (!is.null(midpoint)) {
    lim <- stats::quantile(abs(df$value), 0.98, na.rm = TRUE)
    if (!is.finite(lim) || lim <= 0) lim <- max(abs(df$value), na.rm = TRUE)
    p <- p + scale_fill_gradient2(
      low = div_cols[1], mid = "#F6F0D0", high = div_cols[length(div_cols)],
      midpoint = midpoint, limits = c(-lim, lim), oob = scales::squish, na.value = "white"
    )
  } else {
    qs <- stats::quantile(df$value, c(0.02, 0.98), na.rm = TRUE)
    p <- p + scale_fill_gradientn(colours = palette, limits = qs, oob = scales::squish, na.value = "white")
  }
  ggsave(path, p, width = 5.4, height = 4.8, dpi = 320, bg = "white")
  invisible(p)
}

plot_map(log_o3_hb, file.path(out_dir, "map_log_oiii_hbeta.png"), "log [OIII]/Hbeta")
plot_map(log_nii_ha, file.path(out_dir, "map_log_nii_halpha.png"), "log [NII]/Halpha")
plot_map(center_map(o3_v), file.path(out_dir, "map_oiii_velocity_centered.png"), "[OIII] velocity, median centered", midpoint = 0)
plot_map(center_map(ha_v), file.path(out_dir, "map_halpha_velocity_centered.png"), "Halpha velocity, median centered", midpoint = 0)
plot_map(o3_s, file.path(out_dir, "map_oiii_sigma.png"), "[OIII] sigma")
plot_map(agn_seg$cluster_map, file.path(out_dir, "segmentation_agn_diagnostic_n30.png"), "AGN diagnostic Capivara segmentation", categorical = TRUE)
plot_map(oiii_seg$cluster_map, file.path(out_dir, "segmentation_oiii_window_n25.png"), "[OIII] window Capivara segmentation", categorical = TRUE)

bpt_df <- expand.grid(x = seq_len(nx), y = seq_len(ny))
bpt_df$class <- factor(as.vector(bpt_class), levels = c("SF", "Composite", "AGN/LINER"))
p_bpt <- ggplot(bpt_df, aes(x, y, fill = class)) +
  geom_raster() +
  coord_fixed(expand = FALSE) +
  scale_fill_manual(
    values = c(SF = "#2F74B5", Composite = "#F2A541", `AGN/LINER` = "#B84A28"),
    na.value = "white",
    drop = FALSE
  ) +
  labs(title = "BPT-like class from fitted line fluxes", x = NULL, y = NULL, fill = NULL) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.02, size = 12),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.position = "right"
  )
ggsave(file.path(out_dir, "map_bpt_like_class.png"), p_bpt, width = 5.4, height = 4.8, dpi = 320, bg = "white")

map_stack <- array(NA_real_, dim = c(nx, ny, length(maps)))
for (i in seq_along(maps)) map_stack[, , i] <- suppressWarnings(as.numeric(maps[[i]]))
try(FITSio::writeFITSim(map_stack, file.path(out_dir, "manga_8133_12702_agn_maps.fits"), overwrite = TRUE), silent = TRUE)
write.csv(data.frame(channel = seq_along(maps), name = names(maps)), file.path(out_dir, "fits_map_channels.csv"), row.names = FALSE)

message("Wrote outputs to: ", out_dir)
