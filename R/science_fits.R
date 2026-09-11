# FITSio closes a supplied connection on success but not on every error path.
# Own exactly this connection and close it on exit in either case.
.capivara_read_fits <- function(path, ...) {
  con <- file(path, "rb")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  FITSio::readFITS(con, ...)
}

.capivara_read_fits_header <- function(path) {
  con <- if (grepl("[.]gz$", path, ignore.case = TRUE)) gzfile(path, "rb") else file(path, "rb")
  on.exit(close(con), add = TRUE)
  FITSio::parseHdr(FITSio::readFITSheader(con, maxLines = 20000))
}
