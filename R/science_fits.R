# FITSio closes a supplied connection on success but not on every error path.
# Own exactly this connection and close it on exit in either case.
.capivara_read_fits <- function(path, ...) {
  con <- file(path, "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  out <- FITSio::readFITS(con, ...)
  # Only recognised native DRP spectral HDUs carry this survey guarantee.
  value <- function(key) .fits_header_value(out$hdr, key, "")
  if (identical(value("INSTRUME"), "MaNGA") &&
      value("EXTNAME") %in% c("FLUX", "IVAR", "MASK", "LSFPOST", "LSFPRE") &&
      value("CTYPE3") %in% c("WAVE", "WAVE-LOG")) {
    # Native MaNGA keeps WAVE in the sixth FITSio image HDU, including in
    # the audited MEGACUBE containers. Require its identity and length; never
    # approximate a logarithmic grid with a linear header increment.
    wave_hdu <- .capivara_read_fits(path, hdu = 6)
    if (!identical(.fits_header_value(wave_hdu$hdr, "EXTNAME", ""), "WAVE") ||
        length(wave_hdu$imDat) != dim(out$imDat)[3]) {
      stop("Native MaNGA WAVE HDU does not match the spectral cube.")
    }
    out$wavelength <- as.numeric(wave_hdu$imDat)
    .wavelength_axis(out, dim(out$imDat)[3])
    if (value("EXTNAME") == "FLUX") out$lsf_source <- list(
      path = normalizePath(path), drp_version = value("VERSDRP3"),
      reader = "capivara::read_manga_lsf", representation = "native wavelength-dependent sigma_lambda")
    out$wavelength_frame <- "observed"
    out$wavelength_medium <- "vacuum"
    out$wavelength_frame_source <- "native MaNGA DRP spectral HDU; vacuum heliocentric WAVE"
  }
  out
}

.capivara_read_fits_header <- function(path) {
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, "rb") else file(path, "rb")
  on.exit(close(con), add = TRUE)
  FITSio::parseHdr(FITSio::readFITSheader(con, maxLines = 20000))
}
