#!/usr/bin/env Rscript
# Build browser-ready CAPIVARA products from the normalized native FITS.
# The original source cube and intermediate arrays stay out of the artifact.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: Rscript run_observed_pilot.R INPUT_DIRECTORY OUTPUT_DIRECTORY")
src <- normalizePath(args[1],mustWork=TRUE)
out <- args[2]
if(file.exists(out)) stop("Choose a new output directory.")
library(capivara)
m <- jsonlite::read_json(file.path(src,"source.json"),simplifyVector=TRUE)
read_hdu <- function(i,name) {
  f <- FITSio::readFITS(file.path(src,"prepared.fits"),hdu=i,maxLines=20000)
  k <- match("EXTNAME",f$hdr)
  stopifnot(!is.na(k),identical(trimws(f$hdr[k+1L]),name))
  f$imDat
}
f <- read_hdu(1,"FLUX")
v <- read_hdu(2,"VARIANCE")
mask <- read_hdu(3,"VALID")
wave <- as.numeric(read_hdu(4,"WAVE"))
host <- read_hdu(5,"SUPPORT")
stopifnot(length(dim(f))==3L,identical(dim(f),dim(v)),identical(dim(f),dim(mask)),
          identical(dim(host),dim(f)[1:2]),length(wave)==dim(f)[3L],
          all(mask%in%c(0,1)),all(host%in%c(0,1)),sum(host)==m$support_count,
          identical(as.integer(dim(f)),rev(as.integer(m$native_grid_shape))))
# This is an explicit measurement-selected footprint, not a host claim.
d <- list(flux=f,variance=v,valid=mask==1,host_mask=host==1,wavelength=wave,
          wavelength_frame=m$wavelength_frame,wavelength_medium=m$wavelength_medium,
          flux_unit=m$flux_unit,redshift=m$redshift,object_id=m$object_id,
          provenance=m)
result <- analyse_galaxy(d,n_regions=16L,representation=full_flux_representation())
labels <- result$labels$cluster_map
stopifnot(all(is.na(labels)|d$host_mask),
          sum(!is.na(labels))==sum(result$regions$n_spaxels))
t <- result$spectra
ok <- t$n_valid_spaxels==t$n_region_spaxels
stopifnot(all(is.na(t$fixed_aperture_sum[!ok])),
          all(is.na(t$diagonal_variance_only[!ok])),
          any(ok),all(is.finite(t$fixed_aperture_sum[ok])))
# Independently verify the sums against the measured voxel data, not the
# segmentation representation. Use complete channels only.
fm <- matrix(f,ncol=length(wave))
for(id in result$regions$region) {
  idx <- which(as.vector(labels)==id)
  q <- t[t$region==id,,drop=FALSE]
  full <- q$n_valid_spaxels==length(idx)
  stopifnot(isTRUE(all.equal(q$fixed_aperture_sum[full],
                             colSums(fm[idx,full,drop=FALSE]),
                             tolerance=1e-10,check.attributes=FALSE)))
}
# The atlas exporter never claims physical fits. Include original source
# metadata in the portable preview and alongside the catalogue.
result$provenance$input <- m
html <- explore_galaxy(result,output_dir=out,open=FALSE)
# Keep the artifact small: omit the large retained input cube.
unlink(file.path(out,"segmentation.rds"))
file.copy(file.path(src,"source.json"),file.path(out,"source.json"))
json <- jsonlite::read_json(file.path(out,"atlas.json"),simplifyVector=FALSE)
stopifnot(json$object==m$object_id,json$actual_k==nrow(result$regions))
# Preserve the verified dataset provenance separately in the HTML.
h <- paste(readLines(html,warn=FALSE,encoding="UTF-8"),collapse="\n")
logo <- "man/figures/logo_capivara.png"
if(file.exists(logo)) {
  enc <- jsonlite::base64_enc(readBin(logo,"raw",n=file.info(logo)$size))
  h <- sub('src="../reference/figures/logo_capivara.png"',
    paste0('src="data:image/png;base64,',enc,'"'),h,fixed=TRUE)
}
h <- sub('href="gallery.html"',
  'href="https://rafaelsdesouza.com.br/capivara/articles/examples.html"',h,fixed=TRUE)
# Update displayed source and insert the full documented archive/support record
# through textContent, never by interpolating arbitrary HTML.
meta <- jsonlite::base64_enc(charToRaw(enc2utf8(jsonlite::toJSON(m,auto_unbox=TRUE))))
addon <- paste0("<script>const archiveMeta=JSON.parse(new TextDecoder().decode(Uint8Array.from(atob('",
  meta,"'),c=>c.charCodeAt(0))));",
  "D.source=archiveMeta.survey+' · observed pilot';",
  "D.provenance.archive=archiveMeta;load(D);",
  "text('status','Observed pilot · no physical fits');</script>")
h <- sub("</body>",paste0(addon,"</body>"),h,fixed=TRUE)
writeLines(enc2utf8(h),html,useBytes=TRUE)
cat("OBSERVED PILOT:",m$object_id,"regions:",nrow(result$regions),
    "assigned pixels:",sum(!is.na(labels)),
    "complete spectral channels:",sum(ok),"\n")
