#!/usr/bin/env Rscript
# Convert a native DR17 DRP LOGCUBE to the tutorial's prepared RDS contract.
# Rscript prepare_manga_cube.R input.fits[.gz] host_mask.rds z object_id output.rds
# host_mask.rds MUST contain a logical matrix on the exact native FITSio grid.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=5L) stop("Usage: Rscript prepare_manga_cube.R LOGCUBE host_mask.rds z plateifu output.rds")
src <- normalizePath(args[1L],mustWork=TRUE)
z <- suppressWarnings(as.numeric(args[3L]))
if (length(z)!=1L || !is.finite(z) || z<=-1) stop("Supply measured systemic redshift.")
if (file.exists(args[5L])) stop("Output already exists.")
prepare <- function(src,maskfile,z,id,out) {
  p <- src
  if (grepl("[.]gz$",src,ignore.case=TRUE)) {
    p <- tempfile(fileext=".fits")
    on.exit(unlink(p),add=TRUE)
    con <- gzfile(src,"rb"); dest <- file(p,"wb")
    tryCatch(repeat {
      b <- readBin(con,"raw",n=1024L*1024L)
      if (!length(b)) break
      writeBin(b,dest)
    },finally={close(con);close(dest)})
  }
  header_value <- function(h,key) {
    i <- match(key,h)
    if (is.na(i) || i>=length(h)) return(NA_character_)
    trimws(h[i+1L])
  }
  read_hdu <- function(i,name) {
    con <- file(p,"rb")
    on.exit(try(close(con),silent=TRUE),add=TRUE)
    x <- FITSio::readFITS(con,hdu=i,maxLines=20000)
    if (!identical(header_value(x$hdr,"EXTNAME"),name))
      stop("Expected ",name," at FITSio HDU ",i,"; inspect the actual file layout.")
    x
  }
  f <- read_hdu(1L,"FLUX")
  iv <- read_hdu(2L,"IVAR")$imDat
  dq <- read_hdu(3L,"MASK")$imDat
  wave <- as.numeric(read_hdu(6L,"WAVE")$imDat)
  dims <- dim(f$imDat)
  if (length(dims)!=3L || !identical(dim(iv),dims) || !identical(dim(dq),dims) ||
      length(wave)!=dims[3L] || any(!is.finite(wave)) || any(diff(wave)<=0))
    stop("Native arrays or WAVE axis do not agree.")
  host <- readRDS(maskfile)
  if (!is.logical(host) || !identical(dim(host),dims[1:2]) || anyNA(host))
    stop("Host mask must be a verified logical matrix on the native spatial grid.")
  ok <- is.finite(f$imDat) & is.finite(iv) & iv>0 & is.finite(dq) & dq==0
  v <- array(NA_real_,dims); v[ok] <- 1/iv[ok]
  unit <- header_value(f$hdr,"BUNIT")
  if (is.na(unit) || !nzchar(unit)) stop("Missing FLUX units; inspect the FITS header.")
  saveRDS(list(flux=f$imDat,wavelength=wave,wavelength_frame="observed",
    wavelength_medium="vacuum",flux_unit=unit,variance=v,valid=ok,host_mask=host,
    redshift=z,object_id=id,hdr=f$hdr,axDat=f$axDat,
    provenance=list(source="MaNGA DRP LOGCUBE",file=basename(src),
      sha256=digest::digest(file=src,algo="sha256"),
      host_mask_sha256=digest::digest(file=maskfile,algo="sha256"),
      mask_policy="MASK==0, positive finite IVAR, finite FLUX; conservative",
      wavelength="native WAVE extension; vacuum observed Angstrom",
      lsf="Not ingested by this segmentation-only tutorial")),
    out)
}
prepare(src,args[2L],z,args[4L],args[5L])
