# Evaluate four physical-coordinate contracts, then hash the retained evidence.
# Rscript tools/finalize_wavelength_freeze.R PROJECT
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
project<-normalizePath(args[1]);setwd(project)
root<-file.path(project,'results/science_baseline_freeze_v2');physical<-file.path(root,'physical')
validated<-file.path(root,'validated');repo<-file.path(project,'worktrees/capivara-science-baseline')
hash<-function(f)setNames(vapply(f,function(x)digest::digest(file=x,algo='sha256'),character(1)),f)
old<-readRDS(file.path(root,'historical_sha256_before.rds'));stopifnot(identical(old,hash(names(old))))
stopifnot(!length(system2('git',c('-C',shQuote(repo),'diff','670a61857e74a840822ec026c659ba304de19c4d','--','docs/SCIENCE_BASELINE_FREEZE.md'),stdout=TRUE)))
tests<-readRDS(file.path(root,'tests_physical.rds'));backend<-readRDS(file.path(root,'backend_physical_tests.rds'))
stopifnot(sum(tests$failed)==0,!any(tests$error),sum(backend$failed)==0,!any(backend$error))
checks<-readLines(file.path(root,'physical_source_install_checks.txt'))
stopifnot(length(checks)==5,all(startsWith(checks,'PASS')))
checkfiles<-file.path(root,'build_physical',c('capivara.Rcheck/00check.log','capivaraPPXF.Rcheck/00check.log'))
stopifnot(all(vapply(checkfiles,function(f)any(grepl('Status: OK',readLines(f),fixed=TRUE)),logical(1))))
pilots<-c('11004-12701','11014-3704','8602-12705','8932-3701','9869-9102')
ids<-c(pilots,'8941-12704','9039-9101')
comparison<-read.csv(file.path(validated,'comparison_summary.csv'))
regional<-read.csv(file.path(validated,'regional_comparison_summary.csv'))
windows<-read.csv(file.path(physical,'exact_windows.csv'))
lines<-read.csv(file.path(physical,'line_coordinate_checks.csv'))
lsf<-readRDS(file.path(physical,'lsf_contracts.rds'))
orientation<-jsonlite::fromJSON(file.path(physical,'orientation_verification.json'))
registry<-read.csv(file.path(physical,'line_registry.csv'))
kin<-read.csv(file.path(physical,'kinematic_smoke_summary.csv'))
stopifnot(setequal(comparison$plateifu,ids),nrow(windows)==14,nrow(lines)==14,nrow(kin)==4,
  setequal(names(lsf),ids),nrow(registry)==14,all(registry$wavelength_medium=='vacuum'),
  all(nzchar(registry$reference)),registry$rest_wavelength[registry$name=='halpha']==6564.608,
  registry$rest_wavelength[registry$name=='oiii5007']==5008.240,
  orientation$passed,orientation$samples==63,orientation$max_absolute_difference==0)
# Windows are native integer Angstrom samples, not relabelled rest coordinates.
stopifnot(all(abs(windows$rest_min-windows$native_min/(1+windows$redshift))<1e-9),
  all(abs(windows$rest_max-windows$native_max/(1+windows$redshift))<1e-9),
  all(windows$n_channels==windows$last_channel-windows$first_channel+1),
  all(windows$native_min[windows$requested_frame=='observed']==4800),
  all(windows$native_max[windows$requested_frame=='observed']==7400),
  all(windows$rest_min[windows$requested_frame=='rest']>=4800),
  all(windows$rest_max[windows$requested_frame=='rest']<=7400))
stopifnot(all(abs(lines$predicted_observed-lines$lab_vacuum*(1+lines$redshift))<1e-9),
  all(lines$roundtrip_max_angstrom<1e-9),all(lines$air_on_vacuum_bias_kms>80),
  all(lines$actual_velocity_min>=-600),all(lines$actual_velocity_max<=600))
stopifnot(setequal(vapply(lsf,function(x)x$provenance$drp_version,character(1)),c('v2_7_1','v3_1_1')))
for(a in lsf) {
 p<-a$provenance
 stopifnot(p$quantity=='Gaussian sigma_lambda',p$units=='Angstrom',!p$scalar_fallback,
   p$wavelength_medium=='vacuum',p$wavelength_frame=='observed',p$dimensions[3]==6732,
   min(a$pre_range)>.1,max(a$post_range)<10,min(a$post_variance_minus_pre_range)>0,
   a$method_contracts$direct$selected_lsf=='post',a$method_contracts$template$selected_lsf=='pre')
}
bar<-lapply(pilots,function(id)readRDS(file.path(validated,'comparisons',id,'bar_rotation.rds')))
bar_ready<-all(vapply(bar,function(x)x$mask_iou==1&&x$pa_error_deg<1e-6&&x$mask_area_delta==0,logical(1)))
required<-c('line_rest_wavelength','wavelength_medium','input_wavelength_medium','systemic_redshift','observed_line_centre','line_reference')
kin_details<-list()
for(i in seq_len(nrow(kin))) {
 id<-kin$plateifu[i];line<-kin$line[i]
 path<-file.path(physical,'kinematic_smoke',id,line,paste0('v2_',gsub('-','_',id),'_',line,'_capivara_kinematic_results.rds'))
 k<-readRDS(path);p<-k$kinematics$frame_provenance
 stopifnot(all(required %in% names(p)),p$wavelength_medium=='vacuum',p$input_wavelength_medium=='vacuum',
  p$line_rest_wavelength==registry$rest_wavelength[registry$name==line],
  p$lsf$quantity=='Gaussian sigma_lambda',!p$lsf$correction_applied,!p$lsf$scalar_fallback,
  identical(k$path_features$frame_provenance$lsf,p$lsf),
  dim(k$kinematics$lsf_sigma_angstrom_pre)[3]==length(k$kinematics$line_idx),
  identical(dim(k$kinematics$lsf_sigma_angstrom_pre),dim(k$kinematics$lsf_sigma_angstrom_post)))
 kin_details[[i]]<-data.frame(plateifu=id,line=line,measured=sum(k$kinematics$valid),
   assigned=sum(is.finite(k$kinematic_aware$cluster_map)),path_profiles=nrow(k$path_features$table))
}
kin_details<-do.call(rbind,kin_details);write.csv(kin_details,file.path(physical,'kinematic_counts.csv'),row.names=FALSE)
spectral_ready<-all(comparison$same_footprint)&&all(comparison$prior_same[comparison$role=='pilot'])&&
  all(comparison$observed_flux_relative_L1<1e-12&comparison$rest_flux_relative_L1<1e-12)
kinematic_ready<-all(kin$n_path_profiles>0)&all(kin$coordinate_max_error_kms<1e-8)&
  all(abs(kin$local_represented_centroid_kms)<1e-8)&!any(kin$missing_profiles_imputed)
stopifnot(spectral_ready,kinematic_ready,bar_ready)
readiness<-c('READY_FOR_RESTFRAME_SPECTRAL_SEGMENTATION','READY_FOR_SYSTEMIC_FRAME_KINEMATIC_SEGMENTATION','READY_FOR_SANDRA_SPATIAL_ATLAS_BASELINE')
paths<-c(capivara=repo,capivaraPPXF=file.path(project,'worktrees/capivaraPPXF-science-baseline'),spectropath='/Users/rd23aag/Documents/GitHub/spectropath')
pins<-setNames(vapply(names(paths),function(p)strsplit(checks[grepl(paste0('^PASS ',p,' '),checks)][1],' +')[[1]][5],character(1)),names(paths))
source_files<-unlist(lapply(paths,function(p)c(file.path(p,c('DESCRIPTION','NAMESPACE')),list.files(file.path(p,'R'),full.names=TRUE),
  list.files(file.path(p,'inst'),pattern='[.](R|py)$',recursive=TRUE,full.names=TRUE))))
source_files<-c(source_files,list.files(file.path(repo,'tools'),pattern='(wavelength|physical|systemic|lsf|frame)',full.names=TRUE))
artifacts<-c(list.files(physical,recursive=TRUE,full.names=TRUE),list.files(validated,recursive=TRUE,full.names=TRUE),checkfiles,
 file.path(root,c('tests_physical.rds','backend_physical_tests.rds','physical_source_install_checks.txt','lsf_header_inventory.json')),
 file.path(project,c('capivara_0.4.2.9000.tar.gz','capivaraPPXF_0.0.2.9000.tar.gz')))
artifacts<-artifacts[!dir.exists(artifacts)]
manifest<-list(created_utc=format(Sys.time(),tz='UTC'),validated_code_commits=pins,
 report_commit=system2('git',c('-C',shQuote(repo),'rev-parse','HEAD'),stdout=TRUE),
 report_sha256=hash(file.path(repo,'docs/SCIENCE_BASELINE_FREEZE_V2.md')),
 spectral_partition_commit='a243fdfe551a66c32e90ee39e04ce113918fb9aa',
 spectral_equivalence='Final native inputs have identical SHA256 and exact selected channels; clustering implementation unchanged',
 package_versions=vapply(names(paths),function(p)as.character(packageVersion(p)),character(1)),
 historical_files_verified=length(old),historical_sha256=old,source_sha256=hash(source_files),artifact_sha256=hash(artifacts),
 tests=tests,backend_tests=backend,source_install_checks=checks,comparison=comparison,regional_comparison=regional,
 windows=windows,line_coordinates=lines,lsf=lsf,orientation=orientation,kinematic_smoke=kin,bar_rotation=bar,readiness=readiness,
 scope='Physical coordinates, fixed spectral partitions and observed kinematic/path descriptors. Intrinsic dispersions and legacy native-LSF pPXF fitting are not certified.',session=sessionInfo())
saveRDS(manifest,file.path(root,'freeze_manifest_v2.rds'))
writeLines(c(readiness,paste('Historical files unchanged:',length(old)),paste('CAPIVARA expectations passed:',sum(tests$passed)),
 paste('Artifacts hashed:',length(artifacts)),manifest$scope),file.path(root,'freeze_summary.txt'))
cat(readiness,sep='\n')
