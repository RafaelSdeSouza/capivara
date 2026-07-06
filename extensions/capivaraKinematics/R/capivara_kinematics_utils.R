.capivara_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before running the Capivara kinematic pipeline.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.capivara_dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

#' Convert a matrix to long plotting form
#'
#' @param mat Matrix to convert.
#' @param value_name Name for the value column.
#' @return A data frame with x, y, and value columns.
#' @export
matrix_to_long <- function(mat, value_name = "value") {
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  out <- expand.grid(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat))
  )
  idx <- cbind(out$y, out$x)
  out[[value_name]] <- mat[idx]
  out
}

.fits_header_value <- function(hdr, key, default = NA_character_) {
  idx <- which(hdr == key)[1]
  if (!length(idx) || is.na(idx) || idx >= length(hdr)) {
    return(default)
  }
  value <- hdr[idx + 1L]
  if (is.na(value) || !nzchar(value)) default else value
}

.fits_temp_unzip <- function(path) {
  if (!grepl("\\.gz$", path, ignore.case = TRUE)) {
    return(path)
  }
  tmp <- tempfile(fileext = ".fits")
  con <- gzfile(path, "rb")
  on.exit(close(con), add = TRUE)
  raw <- readBin(con, what = "raw", n = file.info(path)$size * 30L)
  writeBin(raw, tmp)
  tmp
}

.fits_hdu_index <- function(path, max_hdu = 80L) {
  rows <- vector("list", max_hdu)
  for (hdu in seq_len(max_hdu)) {
    z <- try(FITSio::readFITS(path, hdu = hdu, maxLines = 20000), silent = TRUE)
    if (inherits(z, "try-error")) {
      rows[[hdu]] <- NULL
      next
    }
    rows[[hdu]] <- data.frame(
      hdu = hdu,
      extname = .fits_header_value(z$hdr, "EXTNAME"),
      dims = paste(dim(z$imDat), collapse = "x"),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.fits_read_ext <- function(path, hdu_index, extname) {
  hit <- which(toupper(hdu_index$extname) == toupper(extname))[1]
  if (!length(hit) || is.na(hit)) {
    stop("Could not find FITS extension '", extname, "'.", call. = FALSE)
  }
  FITSio::readFITS(path, hdu = hdu_index$hdu[hit], maxLines = 20000)
}

.fits_channel_index <- function(hdr, patterns) {
  keys <- grep("^C[0-9]+$", hdr, value = TRUE)
  if (!length(keys)) {
    return(NA_integer_)
  }
  labels <- vapply(keys, function(k) .fits_header_value(hdr, k), character(1))
  hit <- grep(paste(patterns, collapse = "|"), labels, ignore.case = TRUE)[1]
  if (!length(hit) || is.na(hit)) {
    return(NA_integer_)
  }
  as.integer(sub("^C0*", "", keys[hit]))
}

.matrix_from_fits_image <- function(x) {
  mat <- x$imDat
  if (length(dim(mat)) > 2L) {
    stop("Expected a 2-D FITS image but got dimensions: ",
         paste(dim(mat), collapse = "x"), call. = FALSE)
  }
  as.matrix(mat)
}
