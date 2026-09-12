# Reproduce the old contract from the preserved V1 installation, never new source.
# Rscript tools/verify_historical_wavelength_defect.R PROJECT OUTPUT
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
project<-normalizePath(args[1]);out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(file.path(project,'results/science_baseline_freeze/library'),.libPaths()))
library(capivara);stopifnot(packageVersion('capivara')=='0.4.0.9000',
  !'feature_wavelength_frame' %in% names(formals(segment_large)))
m<-as.data.frame(arrow::read_parquet(file.path(project,'catalogues/bars_v1/input_manifest.parquet')))
row<-m[which(m$plateifu=='11004-12701'),];path<-file.path(project,row$source_relative_path)
readfits<-get('.capivara_read_fits',asNamespace('capivara'))
cube<-readfits(path,hdu=1);wave<-as.numeric(readfits(path,hdu=6)$imDat)
stopifnot(isTRUE(all.equal(as.numeric(get('.wavelength_axis',asNamespace('capivara'))(cube,length(wave))),wave,tolerance=1e-9)))
x<-cube;x$imDat<-x$imDat[35:38,35:38,,drop=FALSE]
a<-segment_large(x,Ncomp=2,redshift=row$redshift,feature_wavelength_range=c(4800,7400),knn_k=8)
b<-segment(x,Ncomp=2,redshift=row$redshift,feature_wavelength_range=c(4800,7400))
ix<-which(wave>=4800&wave<=7400)
stopifnot(identical(a$feature_wavelength_index,ix),identical(b$feature_wavelength_index,ix))
write.csv(data.frame(channel_index=ix,native_wavelength=wave[ix],rest_wavelength=wave[ix]/(1+row$redshift)),
  file.path(out,'selected_channels.csv'),row.names=FALSE)
saveRDS(list(version=packageVersion('capivara'),library=find.package('capivara'),
  validated_V1_code='2cdce2f54526ecc9b04eb3aa7f2a54b51f1e7100',redshift=row$redshift,
  selected_native_range=range(wave[ix]),actual_rest_range=range(wave[ix])/(1+row$redshift),
  required_observed_range=c(4800,7400)*(1+row$redshift),flux_header=cube$hdr),file.path(out,'proof.rds'))
cat('Historical native:',range(wave[ix]),'; actual rest:',range(wave[ix])/(1+row$redshift),'\n')
