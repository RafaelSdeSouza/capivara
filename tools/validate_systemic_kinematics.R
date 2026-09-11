# Public kinematic API smoke test; two fixed pilots, Halpha and [O III] 5007.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3)
project<-normalizePath(args[1]);output<-normalizePath(args[2]);.libPaths(c(normalizePath(args[3]),.libPaths()))
library(capivara);stopifnot(packageVersion('capivara')=='0.4.1.9000')
ns<-asNamespace('capivara');getf<-function(n)get(n,ns)
manifest<-as.data.frame(arrow::read_parquet(file.path(project,'catalogues/bars_v1/input_manifest.parquet')))
rows<-list();offset_rows<-list();counter<-0
for(id in c('11004-12701','8602-12705'))for(line in c('halpha','oiii5007')) {
  row<-manifest[which(manifest$plateifu==id),];path<-file.path(project,row$source_relative_path)
  dest<-file.path(output,'kinematic_smoke',id,line);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  message('PUBLIC API ',id,' ',line,' ',Sys.time())
  run<-segment_kinematics(path,redshift=row$redshift,emission_line=line,
    segmentation_mode='path_signature',output_dir=dest,object_id=paste0('v2_',id),
    knn_k=40,n_segments=18,n_path_segments=45,wavelength_frame='observed',
    wavelength_medium='vacuum',profile_centering_mode='systemic',show_plots=FALSE)
  native<-run$native;k<-native$kinematics;pf<-native$path_features
  stopifnot(k$frame_provenance$systemic_redshift==row$redshift,
    k$frame_provenance$systemic_redshift_source=='manual',
    k$frame_provenance$observed_line_centre==k$rest_wave*(1+row$redshift),
    identical(k$velocity,k$velocity_systemic),sum(k$valid)>0,nrow(pf$table)>0,
    pf$frame_provenance$profile_centering_mode=='systemic')
  # A native measured profile demonstrates both representation modes, without
  # resampling the cube or requiring a locally centred extraction window.
  fits<-getf('.capivara_read_fits')(path,hdu=1)
  wave<-as.numeric(getf('.capivara_read_fits')(path,hdu=6)$imDat)
  t<-pf$table[which.max(pf$table$line_flux),]
  spectrum<-fits$imDat[t$x,t$y,];small<-array(spectrum,c(1,1,length(spectrum)))
  local<-getf('build_path_feature_cube')(small,wave,k$lambda0,matrix(TRUE,1,1),
    rest_wave=k$rest_wave,redshift=row$redshift,line_name=line,
    systemic_redshift_source=row$redshift_source,profile_centering_mode='local_centroid')
  systemic<-getf('build_path_feature_cube')(small,wave,k$lambda0,matrix(TRUE,1,1),
    rest_wave=k$rest_wave,redshift=row$redshift,line_name=line,
    systemic_redshift_source=row$redshift_source,profile_centering_mode='systemic')
  stopifnot(nrow(local$table)==1,abs(local$table$represented_centroid_kms)<1e-8,
    identical(local$selected_channel_indices,systemic$selected_channel_indices))
  ix<-local$selected_channel_indices;v<-local$velocity
  flux<-getf('baseline_subtract')(v,spectrum[ix]);amp<-max(abs(flux))
  represented<-getf('.centre_velocity_profile')(v,flux/amp,'local_centroid')
  profiles<-data.frame(channel=ix,native_wavelength=wave[ix],systemic_velocity_kms=v,
    local_velocity_kms=represented$path[,1],native_flux=spectrum[ix],continuum_subtracted_flux=flux)
  write.csv(profiles,file.path(dest,'representative_profiles.csv'),row.names=FALSE)
  conventional<-spectropath::classical_features(cbind(v,flux))
  saveRDS(list(systemic=systemic,local=local,conventional=conventional,
    spaxel=c(row=t$x,column=t$y),flux_units='native F_lambda; coordinate integrals are F_lambda * km/s, not calibrated line flux'),file.path(dest,'representation_smoke.rds'))
  # Same physical offsets, at each actual pilot redshift and each line.
  vv<-c(-300,-100,0,100,300)
  synthetic_wave<-k$rest_wave*(1+row$redshift)*(1+vv/299792.458)
  recovered<-getf('.systemic_line_coordinate')(synthetic_wave,k$rest_wave,row$redshift,
    'observed',line,row$redshift_source,600,'systemic','vacuum')$velocity
  stopifnot(max(abs(recovered-vv))<1e-8)
  counter<-counter+1;offset_rows[[counter]]<-data.frame(plateifu=id,line=line,z=row$redshift,
    expected_velocity=vv,native_wavelength=synthetic_wave,recovered_velocity=recovered)
  old_delta<-NA_real_
  if(line=='halpha') {
    old<-readRDS(file.path(project,'results/science_baseline_freeze/corrected/pilots',id,'kinematics_preview',
      paste0('baseline_',gsub('-','_',id),'_halpha_capivara_kinematic_results.rds')))
    common<-k$valid & old$kinematics$valid
    old_delta<-median(k$velocity[common]-old$kinematics$velocity[common],na.rm=TRUE)
  }
  rows[[counter]]<-data.frame(plateifu=id,line=line,redshift=row$redshift,rest_wavelength=k$rest_wave,
    observed_line_centre=k$lambda0,n_valid=sum(k$valid),n_path_profiles=nrow(pf$table),
    n_kinematic_segments=run$native$kinematic_aware$Ncomp,n_path_segments=run$segmentation$Ncomp,
    median_systemic_velocity_kms=median(k$velocity[k$valid]),
    median_delta_from_historical_kms=old_delta,representative_centroid_kms=systemic$table$represented_centroid_kms,
    local_represented_centroid_kms=local$table$represented_centroid_kms,
    coordinate_max_error_kms=max(abs(recovered-vv)),missing_profiles_imputed=FALSE)
  write.csv(do.call(rbind,rows),file.path(output,'kinematic_smoke_summary.csv'),row.names=FALSE)
  write.csv(do.call(rbind,offset_rows),file.path(output,'physical_velocity_offsets.csv'),row.names=FALSE)
  print(rows[[counter]]);rm(run,native,k,pf,fits,small,local,systemic);gc()
}
