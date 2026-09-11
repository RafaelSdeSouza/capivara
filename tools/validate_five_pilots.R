# Fixed-five correctness validation, never the Sandra-46 runner.
# Usage: Rscript tools/validate_five_pilots.R PROJECT OUTPUT LIBRARY SPS_FILE [REFERENCE_RDS]
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)%in%c(4L,5L)) stop("Supply PROJECT OUTPUT LIBRARY SPS_FILE [REFERENCE_RDS].")
project <- normalizePath(args[1],mustWork=TRUE)
output <- args[2]
.libPaths(c(normalizePath(args[3],mustWork=TRUE),.libPaths()))
sps_file <- normalizePath(args[4],mustWork=TRUE)
stopifnot(packageVersion("capivara")=="0.4.0.9000",
          packageVersion("capivaraPPXF")=="0.0.1.9000")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
ids <- c("11004-12701","11014-3704","8602-12705","8932-3701","9869-9102")
manifest <- arrow::read_parquet(file.path(project,"catalogues/bars_v1/input_manifest.parquet"))
manifest <- manifest[match(ids,manifest$plateifu),]
stopifnot(identical(as.character(manifest$plateifu),ids),
          all(manifest$science_input_policy=="MANGA_DRP_HDU_SUBSET_IN_MEGACUBE"),
          !any(manifest$megacube_derived_products_used))
ns <- asNamespace("capivara")
read_fits <- get(".capivara_read_fits",ns)
detect_bar <- get("detect_bar",ns)
rot <- function(m)t(m[nrow(m):1,,drop=FALSE])
angle_diff <- function(a,b)abs(((a-b+90)%%180)-90)
iou <- function(a,b)if (sum(a|b)>0) sum(a&b)/sum(a|b) else NA_real_
partition <- function(x) {
  x[!is.finite(x)|x<=0] <- NA
  match(as.vector(x),unique(as.vector(x)))
}
collapse_ids <- function(x)if(length(x))paste(x,collapse=",") else "none"
fit_ok <- function(x) !is.na(x) & x
starlet_args <- list(starlet_J=5,starlet_scales=2:5,include_coarse=FALSE,denoise_k=1,
                     mode="soft",positive_only=TRUE)
# Snapshot old products before any validation writes.
baseline <- if(length(args)==5L)readRDS(args[5]) else lapply(ids,function(id) {
  base <- file.path(project,"results/bars_v1")
  seg <- readRDS(file.path(base,"segmentation",id,"segmentation.rds"))
  popfile <- file.path(base,"stellar_populations",id,"stellar_populations.parquet")
  list(bar=readRDS(file.path(base,"bar_detection",id,"bar_detection.rds")),
       cluster_map=seg$cluster_map,
       segmentation_metadata=seg$ifun_metadata,
       populations=if(file.exists(popfile))arrow::read_parquet(popfile) else NULL)
})
names(baseline) <- ids
saveRDS(baseline,file.path(output,"pre_fix_reference.rds"))
writeLines(capture.output(sessionInfo()),file.path(output,"sessionInfo.txt"))
saveRDS(list(manifest=manifest,sps_md5=tools::md5sum(sps_file),
             packages=lapply(c("capivara","capivaraPPXF"),function(p)
               list(package=p,version=as.character(packageVersion(p)),
                    installed_path=find.package(p))),
             starlet_args=starlet_args,created_utc=format(Sys.time(),tz="UTC")),
        file.path(output,"run_provenance.rds"))
qc_plot <- function(id,bar,map,file) {
  grDevices::png(file,width=1500,height=650,res=150)
  on.exit(grDevices::dev.off(),add=TRUE)
  graphics::par(mfrow=c(1,2),mar=c(3,3,2,1))
  for (kind in c("light","segments")) {
    m <- if(kind=="light")bar$scores$support$collapsed else map
    if(is.null(m)) m <- if(kind=="light")bar$bar_score else map
    if(kind=="light")m <- asinh(pmax(m,0)/max(stats::median(m[m>0],na.rm=TRUE),1e-9))
    graphics::image(seq_len(ncol(m)),seq_len(nrow(m)),t(m),
                    col=if(kind=="light")grDevices::hcl.colors(128,"Grays") else grDevices::hcl.colors(20,"Dark 3"),
                    asp=1,xlab="column",ylab="row",main=paste(id,kind),useRaster=TRUE)
    if(any(bar$bar_mask))graphics::contour(seq_len(ncol(m)),seq_len(nrow(m)),
                       t(bar$bar_mask*1),levels=.5,add=TRUE,drawlabels=FALSE,col="#D55E00",lwd=1.5)
    # Only draw the actual fitted image PA; never plot the deprojected phi.
    if(nrow(bar$diagnostics)) {
      pa <- bar$diagnostics$pa_image_deg[1]*pi/180
      centre <- bar$scores$center
      if(is.null(centre))centre <- (dim(m)+1)/2
      radius <- bar$diagnostics$bar_radius_px[1]
      graphics::segments(centre[2]-radius*sin(pa),centre[1]-radius*cos(pa),
                         centre[2]+radius*sin(pa),centre[1]+radius*cos(pa),col="#0072B2",lwd=1.5)
    }
  }
}
results <- list()
for(id in ids) {
  message("\n===== Fixed pilot ",id," =====")
  run_one <- function() {
    row <- manifest[manifest$plateifu==id,]
    path <- file.path(project,row$source_relative_path)
    actual_sha256 <- digest::digest(file=path,algo="sha256")
    stopifnot(identical(actual_sha256,as.character(row$sha256)))
    dest <- file.path(output,"pilots",id)
    dir.create(dest,recursive=TRUE,showWarnings=FALSE)
    cube <- read_fits(path,hdu=1)
    cube$imDat[!is.finite(cube$imDat)] <- NA_real_
    wave <- as.numeric(read_fits(path,hdu=6)$imDat)
    stopifnot(length(wave)==dim(cube$imDat)[3],all(diff(wave)>0))
    # Preserve segmentation's original header-based wavelength path, while
    # checking it against the explicitly identified native WAVE extension.
    oldwave <- get(".wavelength_axis",ns)(cube,dim(cube$imDat)[3])
    stopifnot(isTRUE(all.equal(as.numeric(oldwave),wave,tolerance=1e-9)))
    params <- baseline[[id]]$bar$parameters
    detect_args <- c(list(input=cube,starlet_args=starlet_args),params)
    message("Automatic bar detection")
    bar <- do.call(detect_bar,detect_args)
    saveRDS(bar,file.path(dest,"bar_detection.rds"))
    message("90-degree raw-cube rotation test")
    rotated <- cube
    rotated$imDat <- aperm(cube$imDat[dim(cube$imDat)[1]:1,,,drop=FALSE],c(2,1,3))
    rb <- do.call(detect_bar,c(list(input=rotated,starlet_args=starlet_args),params))
    rotation <- list(pa_error_deg=angle_diff(rb$diagnostics$pa_image_deg[1],
                                            bar$diagnostics$pa_image_deg[1]+90),
                     mask_iou=iou(rb$bar_mask,rot(bar$bar_mask)),
                     mask_area_delta=sum(rb$bar_mask)-sum(bar$bar_mask),
                     accepted_same=identical(rb$bar_like,bar$bar_like),
                     rotated_diagnostics=rb$diagnostics)
    saveRDS(rotation,file.path(dest,"rotation_qc.rds"))
    rm(rotated,rb); gc()
    message("Full-spectrum sparse Ward segmentation")
    seg <- capivara::segment_large(cube,Ncomp=18L,redshift=row$redshift,
                feature_wavelength_range=c(4800,7400),knn_k=40L,auto_k=TRUE,
                max_k=100L,feature_scale="robust_col",spatial_weight=.15,
                mask=bar$scores$support_mask,valid_mode="signal",
                return_details=FALSE,verbose=TRUE)
    same_partition <- identical(partition(seg$cluster_map),partition(baseline[[id]]$cluster_map))
    same_footprint <- identical(is.finite(seg$cluster_map)&seg$cluster_map>0,
                               is.finite(baseline[[id]]$cluster_map)&baseline[[id]]$cluster_map>0)
    saveRDS(list(cluster_map=seg$cluster_map,Ncomp=seg$Ncomp,backend=seg$backend,
                 feature_wavelengths=seg$feature_wavelengths,
                 same_partition=same_partition,same_footprint=same_footprint),
            file.path(dest,"segmentation_validation.rds"))
    labels <- sort(unique(seg$cluster_map[is.finite(seg$cluster_map)&seg$cluster_map>0]))
    overlap <- do.call(rbind,lapply(labels,function(k) {
      support <- seg$cluster_map==k & is.finite(seg$cluster_map)
      data.frame(bin=k,n_spaxels=sum(support),n_bar=sum(support&bar$bar_mask),
                 bar_fraction=sum(support&bar$bar_mask)/sum(support),
                 mask_source="automatic_CAPIVARA_candidate_not_independent_truth")
    }))
    saveRDS(overlap,file.path(dest,"segment_bar_overlap.rds"))
    qc_plot(id,bar,seg$cluster_map,file.path(dest,"bar_segmentation_qc.png"))
    message("Native IVAR and MASK; summed region spectra")
    ivar <- read_fits(path,hdu=2)$imDat
    mask <- read_fits(path,hdu=3)$imDat
    stopifnot(identical(dim(ivar),dim(cube$imDat)),identical(dim(mask),dim(ivar)))
    valid <- is.finite(ivar)&ivar>0&is.finite(mask)&mask==0&is.finite(cube$imDat)
    variance <- array(NA_real_,dim(ivar)); variance[valid] <- 1/ivar[valid]
    seg$original_cube$imDat[!valid] <- NA_real_
    seg$original_cube$wavelength <- wave
    input <- capivaraPPXF::as_ppxf_input(seg,var_cube=variance,spectrum="sum",
          redshift=row$redshift,metadata=list(input_hdus=c("FLUX","IVAR","MASK","WAVE"),
          input_sha256=row$sha256,mask_policy="MASK==0 and finite positive IVAR",
          variance_model="diagonal native IVAR; spatial covariance not included",
          segmentation_changed=!same_partition,megacube_derived_products_used=FALSE))
    saveRDS(input,file.path(dest,"ppxf_input.rds"))
    rm(ivar,mask,valid,variance,seg,cube);gc()
    message("Fit all 18 regions with paired log grids and propagated variance")
    fit <- capivaraPPXF::fit_ppxf_population(input,sps_file=sps_file,
             redshift=row$redshift,lam_range_rest=c(4800,7400),fwhm_gal=2.76,
             mdegree=8L,quiet=TRUE)
    saveRDS(fit,file.path(dest,"population_fits.rds"))
    message("Explicit preview kinematics consuming the automatic mask and image PA")
    candidate <- list(vetted=FALSE,source="automatic CAPIVARA detector candidate",
                pa_image_deg=bar$diagnostics$pa_image_deg[1],
                axis_ratio=bar$diagnostics$axis_ratio[1],bar_mask=bar$bar_mask)
    kin <- capivara::run_kinematic_analysis(path,redshift=row$redshift,
          model="bisymmetric_bar",emission_line="halpha",object_id=paste0("baseline_",id),
          output_dir=file.path(dest,"kinematics_preview"),knn_k=40L,n_segments=18L,
          model_control=list(analysis_mode="preview",bar_geometry=candidate),
          show_plots=FALSE)
    kr <- kin$model_result
    mask_shared <- identical(kr$bar_geometry$bar_mask,bar$bar_mask)
    pa_shared <- isTRUE(all.equal(kr$bar_geometry$pa_image_deg,candidate$pa_image_deg))
    kin_summary <- list(geometry=kr$geometry,bar_geometry=kr$bar_geometry,
                        qc=kr$qc,config=kr$config,diagnostics=kr$diagnostics,
                        mask_shared=mask_shared,pa_shared=pa_shared)
    saveRDS(kin_summary,file.path(dest,"kinematic_validation.rds"))
    oldbar <- baseline[[id]]$bar
    oldfit <- baseline[[id]]$populations
    summary <- data.frame(plateifu=id,status="completed",n_segments=length(labels),
       segmentation_success=length(labels)==18L,same_partition=same_partition,
       same_footprint=same_footprint,bar_accepted=isTRUE(bar$bar_like),
       pa_image_deg=bar$diagnostics$pa_image_deg[1],
       axis_ratio=bar$diagnostics$axis_ratio[1],bar_radius_px=bar$diagnostics$bar_radius_px[1],
       bar_area=sum(bar$bar_mask),old_mask_iou=iou(bar$bar_mask,oldbar$bar_mask),
       old_pa_delta_deg=angle_diff(bar$diagnostics$pa_image_deg[1],oldbar$diagnostics$pa_deg[1]),
       rotation_pa_error_deg=rotation$pa_error_deg,rotation_mask_iou=rotation$mask_iou,
       rotation_mask_area_delta=rotation$mask_area_delta,
       pa_mask_rotation_pass=rotation$pa_error_deg<1e-6 && rotation$mask_iou>0.999999,
       n_fits=nrow(fit),n_numerical_success=sum(fit_ok(fit$numerical_success)),
       n_converged=sum(fit_ok(fit$converged)),n_bound=sum(fit_ok(fit$parameter_on_bound)),
       n_fit_level_eligible=sum(fit_ok(fit$scientifically_usable)),
       failed_bins=collapse_ids(fit$bin[!fit_ok(fit$numerical_success)]),
       bound_bins=collapse_ids(fit$bin[fit_ok(fit$parameter_on_bound)]),
       old_n_fit_ok=if(is.null(oldfit))NA_integer_ else sum(fit_ok(oldfit$fit_ok)),
       median_chi2=stats::median(fit$chi2,na.rm=TRUE),
       kinematic_numerical_success=isTRUE(kr$qc$numerical_success),
       kinematic_converged=isTRUE(kr$qc$converged),
       kinematic_instability=isTRUE(kr$qc$model_instability_flag),
       kinematic_scientifically_usable=isTRUE(kr$qc$scientifically_usable),
       disc_pa_image_deg=kr$geometry$pa_image_deg,
       preview_phi_bar_disc_deg=kr$bar_geometry$phi_bar_disc_deg,
       kinematic_mask_shared=mask_shared,kinematic_pa_shared=pa_shared,
       error=NA_character_)
    saveRDS(summary,file.path(dest,"summary.rds"))
    print(summary)
    summary
  }
  pilot_warnings <- character()
  results[[id]] <- tryCatch(withCallingHandlers(run_one(),warning=function(w) {
    pilot_warnings <<- c(pilot_warnings,conditionMessage(w))
    invokeRestart("muffleWarning")
  }),error=function(e) {
    message("PILOT FAILED: ",conditionMessage(e))
    data.frame(plateifu=id,status="failed",error=conditionMessage(e))
  })
  saveRDS(unique(pilot_warnings),file.path(output,paste0(id,"_warnings.rds")))
  saveRDS(results,file.path(output,"five_pilot_results.rds"))
  gc()
}
if(any(vapply(results,function(x)x$status!="completed",logical(1))))quit(status=1L)
if(any(vapply(results,function(x)!isTRUE(x$segmentation_success && x$pa_mask_rotation_pass &&
      x$kinematic_mask_shared && x$kinematic_pa_shared),logical(1))))quit(status=2L)
