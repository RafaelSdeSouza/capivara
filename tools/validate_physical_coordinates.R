# Seven previously selected objects only; native physical-coordinate audit.
# Rscript tools/validate_physical_coordinates.R PROJECT OUTPUT LIBRARY
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3)
project<-normalizePath(args[1]);output<-normalizePath(args[2]);.libPaths(c(normalizePath(args[3]),.libPaths()))
library(capivara);stopifnot(packageVersion('capivara')=='0.4.2.9000')
getf<-function(n)get(n,asNamespace('capivara'))
manifest<-as.data.frame(arrow::read_parquet(file.path(project,'catalogues/bars_v1/input_manifest.parquet')))
previous<-file.path(project,'results/science_baseline_freeze_v2/validated')
ids<-read.csv(file.path(previous,'comparison_summary.csv'))$plateifu
registry<-emission_lines('vacuum');write.csv(registry,file.path(output,'line_registry.csv'),row.names=FALSE)
write.csv(emission_lines('air'),file.path(output,'line_registry_air.csv'),row.names=FALSE)
rows<-list();windows<-list();lsfrows<-list();samples<-list();details<-list()
for(id in ids) {
 message('PHYSICAL ',id,' ',Sys.time())
 row<-manifest[which(manifest$plateifu==id),];path<-file.path(project,row$source_relative_path);z<-row$redshift
 stopifnot(identical(digest::digest(file=path,algo='sha256'),row$sha256))
 cube<-getf('.capivara_read_fits')(path,hdu=1);wave<-cube$wavelength
 lsf<-read_manga_lsf(path,cube)
 stopifnot(identical(dim(lsf$lsf_sigma_angstrom_pre),dim(cube$imDat)),identical(lsf$wavelength,wave))
 support<-readRDS(file.path(previous,'comparisons',id,'bar_detection.rds'))$scores$support_mask
 for(frame in c('observed','rest')) {
  old<-readRDS(file.path(previous,'comparisons',id,paste0(frame,'_segmentation.rds')))$wavelength_provenance
  current<-getf('.subset_cubedat_wavelength_range')(cube,c(4800,7400),redshift=z,
    wavelength_frame='observed',feature_wavelength_frame=frame)$provenance
  current$systemic_redshift_source <- row$redshift_source
  current$input_sha256 <- row$sha256
  stopifnot(identical(current,old))
  windows[[paste(id,frame)]]<-data.frame(plateifu=id,redshift=z,requested_frame=frame,
    native_min=old$selected_native_wavelength_min,native_max=old$selected_native_wavelength_max,
    rest_min=old$selected_rest_wavelength_min,rest_max=old$selected_rest_wavelength_max,
    first_channel=min(old$selected_channel_indices),last_channel=max(old$selected_channel_indices),n_channels=old$number_of_selected_channels)
 }
 methods<-list(direct=select_manga_lsf(lsf,'direct_profile','point_sampled'),
   template=select_manga_lsf(lsf,'template_convolution','pixel_integrated'))
 stopifnot(methods$direct$provenance$selected_lsf=='post',methods$template$provenance$selected_lsf=='pre')
 for(line in c('halpha','oiii5007')) {
  lab<-registry$rest_wavelength[registry$name==line]
  co<-getf('.systemic_line_coordinate')(wave,lab,z,'observed',line,row$redshift_source,600,'systemic','vacuum')
  ix<-co$wavelength_provenance$selected_channel_indices
  roundtrip<-co$provenance$observed_line_centre*(1+co$velocity/299792.458)
  stopifnot(max(abs(roundtrip-wave))<1e-9)
  airlab<-emission_lines('air')$rest_wavelength[registry$name==line]
  rows[[paste(id,line)]]<-data.frame(plateifu=id,redshift=z,line=line,lab_vacuum=lab,
    predicted_observed=lab*(1+z),requested_velocity_min=-600,requested_velocity_max=600,
    native_min=min(wave[ix]),native_max=max(wave[ix]),actual_velocity_min=min(co$velocity[ix]),actual_velocity_max=max(co$velocity[ix]),
    n_channels=length(ix),roundtrip_max_angstrom=max(abs(roundtrip-wave)),
    air_on_vacuum_bias_kms=299792.458*(lab/airlab-1))
  channel<-which.min(abs(wave-lab*(1+z)))
  for(kind in c('pre','post')) {
   a<-lsf[[paste0('lsf_sigma_angstrom_',kind)]][,,channel][support]
   a<-a[is.finite(a)&a>0]
   stopifnot(length(a)>0,min(a)>.1,max(a)<10)
   q<-quantile(a,c(0,.5,1),names=FALSE)
   lsfrows[[paste(id,line,kind)]]<-data.frame(plateifu=id,drp=lsf$provenance$drp_version,line=line,kind=kind,
     native_extension=lsf$provenance$native_extension[kind],native_channel=channel,wavelength=wave[channel],
     sigma_min=q[1],sigma_median=q[2],sigma_max=q[3],sigma_kms_median=299792.458*q[2]/wave[channel],
     resolving_power_median=wave[channel]/(2*sqrt(2*log(2))*q[2]),n_support_valid=length(a))
  }
 }
 # Asymmetric (x,y) samples allow independent Astropy order verification.
 pos<-arrayInd(which(support),dim(support));pos<-pos[unique(round(seq(1,nrow(pos),length.out=3))),,drop=FALSE]
 profile<-data.frame(wavelength=wave,pre=lsf$lsf_sigma_angstrom_pre[pos[2,1],pos[2,2],],
   post=lsf$lsf_sigma_angstrom_post[pos[2,1],pos[2,2],])
 write.csv(profile,file.path(output,paste0(id,'_lsf_profile.csv')),row.names=FALSE)
 for(j in seq_len(nrow(pos)))for(channel in c(1500L,3000L,5000L)) {
  x<-pos[j,1];y<-pos[j,2]
  samples[[length(samples)+1L]]<-data.frame(plateifu=id,x=x,y=y,channel=channel,wavelength=wave[channel],
    flux=cube$imDat[x,y,channel],pre=lsf$lsf_sigma_angstrom_pre[x,y,channel],post=lsf$lsf_sigma_angstrom_post[x,y,channel])
 }
 valid<-is.finite(lsf$lsf_sigma_angstrom_pre)&is.finite(lsf$lsf_sigma_angstrom_post)
 delta<-lsf$lsf_sigma_angstrom_post[valid]^2-lsf$lsf_sigma_angstrom_pre[valid]^2
 stopifnot(min(delta)>0)
 details[[id]]<-list(provenance=lsf$provenance,pre_range=range(lsf$lsf_sigma_angstrom_pre,na.rm=TRUE),
   post_range=range(lsf$lsf_sigma_angstrom_post,na.rm=TRUE),post_variance_minus_pre_range=range(delta),
   method_contracts=lapply(methods,`[[`,'provenance'))
 rm(cube,lsf,methods,valid,delta);gc()
}
write.csv(do.call(rbind,windows),file.path(output,'exact_windows.csv'),row.names=FALSE)
write.csv(do.call(rbind,rows),file.path(output,'line_coordinate_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,lsfrows),file.path(output,'lsf_line_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,samples),file.path(output,'lsf_orientation_samples.csv'),row.names=FALSE)
saveRDS(details,file.path(output,'lsf_contracts.rds'))
cat('PASS: exact window equivalence, vacuum registry, systemic roundtrips, seven native LSF cubes, both DRP generations\n')
