# Exactly the five frozen pilots and two sample redshift extremes. No production run.
# Rscript tools/validate_wavelength_frames.R PROJECT V2_OUTPUT LIBRARY
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
project <- normalizePath(args[1]); output <- normalizePath(args[2]); lib <- normalizePath(args[3])
.libPaths(c(lib,.libPaths()))
library(capivara)
stopifnot(packageVersion('capivara')=='0.4.1.9000')
ns <- asNamespace('capivara'); getf <- function(n)get(n,ns)
manifest <- as.data.frame(arrow::read_parquet(file.path(project,'catalogues/bars_v1/input_manifest.parquet')))
manifest <- manifest[which(!is.na(manifest$plateifu)),]
stopifnot(nrow(manifest)==46, !anyDuplicated(manifest$plateifu))
pilots <- c('11004-12701','11014-3704','8602-12705','8932-3701','9869-9102')
extremes <- manifest$plateifu[c(which.min(manifest$redshift),which.max(manifest$redshift))]
ids <- unique(c(pilots,extremes));stopifnot(length(ids)==7)
write.csv(manifest[,c('plateifu','redshift','redshift_source','sha256')],file.path(output,'sandra46_redshift_inventory.csv'),row.names=FALSE)
settings <- list(Ncomp=18L,feature_wavelength_range=c(4800,7400),knn_k=40L,
  auto_k=TRUE,max_k=100L,feature_scale='robust_col',spatial_weight=.15,
  valid_mode='signal',return_details=FALSE,verbose=TRUE)
starlet <- list(starlet_J=5,starlet_scales=2:5,include_coarse=FALSE,denoise_k=1,mode='soft',positive_only=TRUE)
historical <- file.path(project,'results/science_baseline_freeze/corrected/pilots')
params <- readRDS(file.path(historical,pilots[1],'bar_detection.rds'))$parameters
saveRDS(list(settings=settings,ids=ids,pilots=pilots,extremes=extremes,
  median_sample_redshift=median(manifest$redshift),library=lib,session=sessionInfo()),file.path(output,'comparison_provenance.rds'))
rot <- function(m)t(m[nrow(m):1,,drop=FALSE])
canonical <- function(m)match(as.vector(m),unique(as.vector(m)))
components <- function(mask) {
  # Four-neighbour spatial fragmentation; diagonal touching is not connectivity.
  nr<-nrow(mask);nc<-ncol(mask); seen<-matrix(FALSE,nr,nc); sizes<-integer()
  for(start in which(mask))if(!seen[start]) {
    stack<-integer(sum(mask));stack[1]<-start;seen[start]<-TRUE;head<-1L;tail<-1L
    while(head<=tail){
      k<-stack[head];head<-head+1L;r<-(k-1L)%%nr+1L;c<-(k-1L)%/%nr+1L
      nb<-c(if(r>1)k-1L,if(r<nr)k+1L,if(c>1)k-nr,if(c<nc)k+nr)
      nb<-nb[mask[nb]&!seen[nb]]
      if(length(nb)){seen[nb]<-TRUE;stack[(tail+1L):(tail+length(nb))]<-nb;tail<-tail+length(nb)}
    }
    sizes<-c(sizes,tail)
  }
  sizes
}
metrics <- function(a,b) {
  common<-is.finite(a)&is.finite(b);t<-table(a[common],b[common]);n<-sum(t)
  choose2<-function(x)x*(x-1)/2; ab<-sum(choose2(t));aa<-sum(choose2(rowSums(t)));bb<-sum(choose2(colSums(t)));ex<-aa*bb/choose2(n)
  ari<-(ab-ex)/((aa+bb)/2-ex)
  entropy<-function(p)-sum(p[p>0]*log(p[p>0]));p<-t/n
  vi<-2*entropy(as.vector(p))-entropy(rowSums(p))-entropy(colSums(p))
  grid<-matrix(seq_along(a),nrow(a),ncol(a))
  u<-c(grid[-nrow(grid),],grid[,-ncol(grid)]);v<-c(grid[-1,],grid[,-1]);keep<-common[u]&common[v]
  ea<-a[u[keep]]!=a[v[keep]];eb<-b[u[keep]]!=b[v[keep]]
  list(ARI=ari,VI_nats=vi,boundary_jaccard=sum(ea&eb)/sum(ea|eb),
    boundary_F1=2*sum(ea&eb)/(sum(ea)+sum(eb)),edge_agreement=mean(ea==eb),
    observed_boundary_edges=sum(ea),rest_boundary_edges=sum(eb),common_spaxels=sum(common),
    same_footprint=identical(is.finite(a),is.finite(b)),contingency=t)
}
rows<-list()
for(id in ids) {
  message('START ',id,' ',Sys.time())
  dest<-file.path(output,'comparisons',id);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  row<-manifest[manifest$plateifu==id,];path<-file.path(project,row$source_relative_path)
  stopifnot(identical(digest::digest(file=path,algo='sha256'),row$sha256))
  cube<-getf('.capivara_read_fits')(path,hdu=1)
  cube$imDat[!is.finite(cube$imDat)]<-NA_real_
  wave<-as.numeric(getf('.capivara_read_fits')(path,hdu=6)$imDat)
  stopifnot(cube$wavelength_frame=='observed',cube$wavelength_medium=='vacuum',
    isTRUE(all.equal(as.numeric(getf('.wavelength_axis')(cube,length(wave))),wave,tolerance=1e-9)))
  # Full native flux and the historical detector settings determine one support
  # mask shared by both partitions. No window-dependent or tuned mask.
  bar<-do.call(getf('detect_bar'),c(list(input=cube,starlet_args=starlet),params))
  saveRDS(bar,file.path(dest,'bar_detection.rds'))
  if(id %in% pilots) {
    oldbar<-readRDS(file.path(historical,id,'bar_detection.rds'))
    stopifnot(identical(bar$scores$support_mask,oldbar$scores$support_mask),identical(bar$bar_mask,oldbar$bar_mask))
    rotated<-cube;rotated$imDat<-aperm(cube$imDat[dim(cube$imDat)[1]:1,,,drop=FALSE],c(2,1,3))
    rb<-do.call(getf('detect_bar'),c(list(input=rotated,starlet_args=starlet),params))
    rotation<-list(mask_iou=sum(rb$bar_mask&rot(bar$bar_mask))/sum(rb$bar_mask|rot(bar$bar_mask)),
      pa_error_deg=abs((rb$diagnostics$pa_image_deg[1]-bar$diagnostics$pa_image_deg[1]-90+90)%%180-90),
      mask_area_delta=sum(rb$bar_mask)-sum(bar$bar_mask))
    saveRDS(rotation,file.path(dest,'bar_rotation.rds'));stopifnot(rotation$mask_iou==1,rotation$pa_error_deg<1e-6,rotation$mask_area_delta==0)
    rm(rotated,rb,oldbar);gc()
  }
  results<-list();conservation<-list();region_tables<-list()
  for(frame in c('observed','rest')) {
    message('SEGMENT ',id,' ',frame)
    seg<-do.call(segment_large,c(list(input=cube,redshift=row$redshift,
      wavelength_frame='observed',feature_wavelength_frame=frame,mask=bar$scores$support_mask),settings))
    stopifnot(seg$Ncomp==18)
    p<-seg$wavelength_provenance
    p$systemic_redshift_source<-row$redshift_source
    p$input_sha256<-row$sha256
    seg$wavelength_provenance<-p
    ix<-p$selected_channel_indices
    write.csv(data.frame(channel_index=ix,native_wavelength=wave[ix],rest_wavelength=wave[ix]/(1+row$redshift)),file.path(dest,paste0(frame,'_channels.csv')),row.names=FALSE)
    # The spectrum comparison uses exactly the finite native FLUX contributors
    # used by the historical segmentation. No new IVAR/MASK cut changes support.
    sums<-summarize_cluster_spectra(seg)
    total<-colSums(sums$sum_spectra,na.rm=TRUE)
    mat<-matrix(cube$imDat,ncol=length(wave));assigned<-is.finite(seg$cluster_map)
    direct<-colSums(mat[as.vector(assigned),,drop=FALSE],na.rm=TRUE)
    delta<-total-direct
    conservation[[frame]]<-list(max_abs=max(abs(delta)),relative_L1=sum(abs(delta))/sum(abs(direct)),
      retained_spaxels=sum(assigned),finite_input_policy='native finite FLUX; no new data-quality cut')
    regions<-do.call(rbind,lapply(sums$cluster_ids,function(k){
      mask<-is.finite(seg$cluster_map)&seg$cluster_map==k;cs<-components(mask)
      data.frame(region=k,n_spaxels=sum(mask),n_bar=sum(mask&bar$bar_mask),
        bar_fraction=sum(mask&bar$bar_mask)/sum(mask),n_components=length(cs),
        largest_component_fraction=max(cs)/sum(cs),singleton_components=sum(cs==1))
    }))
    region_tables[[frame]]<-regions
    write.csv(regions,file.path(dest,paste0(frame,'_regions.csv')),row.names=FALSE)
    write.csv(seg$cluster_map,file.path(dest,paste0(frame,'_map.csv')),row.names=FALSE)
    write.csv(data.frame(wavelength=wave,t(sums$sum_spectra),check.names=FALSE),file.path(dest,paste0(frame,'_regional_sums.csv')),row.names=FALSE)
    saveRDS(list(wavelength=wave,cluster_ids=sums$cluster_ids,sum_spectra=sums$sum_spectra,finite_counts=sums$finite_counts,provenance=p),file.path(dest,paste0(frame,'_regional_sums.rds')))
    seg$original_cube<-NULL
    saveRDS(seg,file.path(dest,paste0(frame,'_segmentation.rds')));results[[frame]]<-seg
    rm(seg,sums,mat);gc()
  }
  a<-results$observed$cluster_map;b<-results$rest$cluster_map;mt<-metrics(a,b)
  prior_same<-if(id %in% pilots) identical(canonical(a),canonical(readRDS(file.path(historical,id,'segmentation_validation.rds'))$cluster_map)) else NA
  if(id %in% pilots)stopifnot(prior_same)
  write.csv(as.matrix(mt$contingency),file.path(dest,'contingency.csv'),row.names=FALSE)
  saveRDS(list(metrics=mt,flux_conservation=conservation,regions=region_tables,prior_same=prior_same),file.path(dest,'comparison.rds'))
  r<-data.frame(plateifu=id,redshift=row$redshift,role=if(id %in% pilots)'pilot' else if(id==extremes[1])'lowest_z' else 'highest_z',
    ARI=mt$ARI,VI_nats=mt$VI_nats,boundary_jaccard=mt$boundary_jaccard,boundary_F1=mt$boundary_F1,
    edge_agreement=mt$edge_agreement,same_footprint=mt$same_footprint,common_spaxels=mt$common_spaxels,prior_same=prior_same,
    observed_channels=results$observed$wavelength_provenance$number_of_selected_channels,
    rest_channels=results$rest$wavelength_provenance$number_of_selected_channels,
    rest_native_min=results$rest$wavelength_provenance$selected_native_wavelength_min,
    rest_native_max=results$rest$wavelength_provenance$selected_native_wavelength_max,
    observed_components=sum(region_tables$observed$n_components),rest_components=sum(region_tables$rest$n_components),
    observed_flux_relative_L1=conservation$observed$relative_L1,rest_flux_relative_L1=conservation$rest$relative_L1)
  rows[[id]]<-r;write.csv(do.call(rbind,rows),file.path(output,'comparison_summary.csv'),row.names=FALSE)
  print(r);rm(cube,bar,results);gc();message('DONE ',id,' ',Sys.time())
}
