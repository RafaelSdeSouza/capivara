#!/usr/bin/env Rscript
# One-galaxy observational pilot; no pPXF, no production website deployment.
# Rscript tools/run_manga_atlas_pilot.R LOGCUBE.fits preflight_directory output_directory
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Expected LOGCUBE, preflight directory, new output directory.")
suppressPackageStartupMessages(library(capivara))
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite required.")
cube <- normalizePath(args[1], mustWork = TRUE)
pre <- normalizePath(args[2], mustWork = TRUE)
out <- args[3]
if (file.exists(out) || dir.exists(out)) stop("Output already exists.")
dir.create(out, recursive = TRUE)
out <- normalizePath(out, mustWork = TRUE)
meta <- jsonlite::read_json(file.path(pre, "preflight.json"), simplifyVector = TRUE)
stopifnot(identical(meta$object, "8602-12705"),
          identical(meta$release, "SDSS DR17"),
          is.finite(meta$redshift), meta$redshift > 0,
          meta$wcs_max_separation_arcsec <= 1e-3,
          identical(meta$cube_sha256, digest::digest(file=cube, algo="sha256")))
id <- meta$object
z <- meta$redshift
dims <- as.integer(c(meta$nx, meta$ny))
aperture <- read.csv(file.path(pre, "aperture.csv"))
stopifnot(nrow(aperture) == prod(dims), !anyDuplicated(aperture[c("axis1","axis2")]),
          all(aperture$axis1 %in% seq_len(dims[1])),
          all(aperture$axis2 %in% seq_len(dims[2])),
          all(aperture$included %in% c(0L, 1L)))
host <- matrix(FALSE, dims[1], dims[2])
host[as.matrix(aperture[c("axis1","axis2")])] <- aperture$included == 1L
stopifnot(sum(host) == meta$aperture_pixels)
saveRDS(host, file.path(out, "aperture.rds"))
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, argv) {
  status <- system2(rscript, shQuote(c(script, argv)))
  if (status != 0L) stop("Subprocess failed: ", basename(script), " (", status, ").")
}
prep <- system.file("tutorials", "prepare_manga_cube.R", package="capivara", mustWork=TRUE)
prepared_path <- file.path(out, "prepared.rds")
run(prep, c(cube, file.path(out, "aperture.rds"), format(z, digits=17), id, prepared_path))
d <- readRDS(prepared_path)
stopifnot(identical(dim(d$flux)[1:2], dims), length(d$wavelength) == meta$n_wave)
w <- read.csv(file.path(pre,"witness.csv"))
ix <- as.matrix(w[c("axis1","axis2","channel")])
same <- function(x,y) isTRUE(all.equal(as.numeric(x),as.numeric(y), tolerance=1e-7))
stopifnot(nrow(w)>=6L, all(d$valid[ix]), same(d$flux[ix],w$flux),
          same(1/d$variance[ix],w$ivar), same(d$wavelength[w$channel],w$wavelength))
d$provenance$aperture <- meta$aperture
d$provenance$drp_dap_wcs_max_separation_arcsec <- meta$wcs_max_separation_arcsec
d$provenance$redshift_source <- meta$redshift_source
d$provenance$dap_reference_sha256 <- meta$maps_sha256
# A new controlled LOGCUBE comparison, not a claim to reproduce frozen MEGACUBE labels.
# Keep complete native spectra in extraction; restrict both segmentation representations
# to the same rest-frame optical feature interval.
shape <- manga_sandra_spectral_shape()
flux <- full_flux_representation(coordinate_domain = shape$coordinate_domain)
k <- 16L
maps <- list()
map_record <- function(key,title,group,kind,unit,m,note) {
  stopifnot(is.matrix(m), is.numeric(m), identical(dim(m),dims))
  m[!is.finite(m)] <- NA_real_
  if (kind=="categorical") stopifnot(all(is.na(m)|(m>0&m==floor(m))))
  list(id=key,title=title,group=group,kind=kind,unit=unit,note=note,
       values=lapply(seq_len(ncol(m)),function(j)as.list(m[,j])))
}
status <- list(preflight="PASS", full_flux="NOT_RUN", spectral_shape="NOT_RUN",
               native_kinematics="NOT_RUN", emission_segmentation="NOT_RUN",
               physical_models="NOT_RUN")
write_status <- function() jsonlite::write_json(status,file.path(out,"status.json"),
                                                auto_unbox=TRUE,pretty=TRUE)
write_status()
for (key in c("full_flux","spectral_shape")) {
  result <- analyse_galaxy(d, n_regions=k,
    representation=if(key=="full_flux") flux else shape)
  labels <- result$labels$cluster_map
  stopifnot(identical(dim(labels),dims),all(is.na(labels[!host])),
            length(unique(labels[!is.na(labels)]))>=k)
  dest <- file.path(out,key)
  exported <- export_atlas(result,dest)
  stopifnot(file.exists(exported[["atlas"]]))
  maps[[length(maps)+1L]] <- map_record(key,
    if(key=="full_flux") "Full-flux segmentation" else "Spectral-shape segmentation",
    "spectral","categorical","Region ID",labels,
    if(key=="full_flux") "CAPIVARA: absolute flux on native rest-frame 4800–7400 Å channels."
    else "CAPIVARA: declared four-window amplitude; measured amplitude S/N ≥ 30.")
  status[[key]] <- list(state="PASS",regions=result$labels$actual_k,
                        assigned_pixels=sum(!is.na(labels)),
                        aperture_pixels=sum(host))
  write_status()
  rm(result); gc(verbose=FALSE)
}
# Carry the independently fitted survey maps alongside CAPIVARA partitions.
reference <- jsonlite::read_json(file.path(pre,"dap_reference.json"), simplifyVector=FALSE)
stopifnot(length(reference)>0L,all(vapply(reference,function(m)
  startsWith(m$id,"dap_") && length(m$values)==dims[2] &&
  all(lengths(m$values)==dims[1]),logical(1))))
maps <- c(maps,reference)
bundle <- list(schema="capivara-galaxy-views-v1",object=id,
  source="SDSS DR17 · CAPIVARA segmentation + DAP survey reference",
  nx=dims[1],ny=dims[2],line="Multiple DAP lines",redshift=z,maps=maps,
  provenance=c(meta,list(segmentation="New LOGCUBE run; not the frozen MEGACUBE replay",
    requested_regions=k,feature_interval_rest_angstrom=c(4800,7400),
    variance="Diagonal only; no spatial/spectral covariance calibration",
    dap_reference="Survey fits, not CAPIVARA fits. No independent gas velocities claimed.",
    physical_models="NOT_RUN")))
write_bundle <- function() {
  ids <- vapply(bundle$maps,function(m)m$id,character(1))
  stopifnot(!anyDuplicated(ids),length(ids)<=100L)
  jsonlite::write_json(bundle,file.path(out,"views.json"),auto_unbox=TRUE,
                      digits=NA,na="null")
}
write_bundle()
# Native single-line kinematics and SpectroPath are separate exploratory products.
# Their native line-support selection is NOT the DAP aperture or a variance-calibrated fit.
tryCatch({
  script <- system.file("tutorials","run_atlas_views.R",package="capivara",mustWork=TRUE)
  run(script,c(cube,id,format(z,digits=17),"oiii5007",file.path(out,"native_views"),as.character(k)))
  nv <- jsonlite::read_json(file.path(out,"native_views","views.json"),simplifyVector=FALSE)
  stopifnot(identical(nv$object,id),nv$nx==dims[1],nv$ny==dims[2],
            identical(nv$provenance$cube_md5,unname(tools::md5sum(cube))))
  native_maps <- lapply(nv$maps,function(m) {
    m$id <- paste0("native_",m$id)
    m$title <- paste("CAPIVARA",m$title)
    m$note <- paste(m$note,"Exploratory native profile measurement; not a DAP fit.")
    m
  })
  bundle$maps <- c(bundle$maps,native_maps)
  bundle$line <- "Native [O III] + DAP lines"
  bundle$provenance$native <- nv$provenance
  write_bundle()
  status$native_kinematics <- list(state="EXPLORATORY",maps=length(native_maps))
},error=function(e) {status$native_kinematics <<- list(state="FAILED",message=conditionMessage(e))})
write_status()
# Line-sensitive segmentation uses strictly quality-masked samples and a declared aperture.
# Disable sampled projection to avoid presenting a different sampled geometry.
tryCatch({
  valid <- d$valid
  for(j in seq_len(dim(valid)[3])) valid[,,j] <- valid[,,j] & host
  cube_input <- list(imDat=d$flux,wavelength=d$wavelength,
    wavelength_frame=d$wavelength_frame,wavelength_medium=d$wavelength_medium,
    hdr=d$hdr,axDat=d$axDat)
  cube_input$imDat[!valid] <- NA_real_
  emission <- segment_emission_lines(cube_input,redshift=z,lines=c("hbeta","oiii5007","halpha","nii6583","sii6716","sii6731"),
    Ncomp=k,feature_mode="windows",continuum_subtract=TRUE,spatial_weight=0,
    valid_mode="finite",max_pixels=Inf,wavelength_frame="observed",wavelength_medium="vacuum",
    systemic_redshift_source=meta$redshift_source)
  m <- emission$cluster_map
  stopifnot(identical(dim(m),dims),all(is.na(m[!host])),any(!is.na(m)))
  bundle$maps[[length(bundle$maps)+1L]] <- map_record("emission_regions",
    "Emission-line segmentation","emission","categorical","Region ID",m,
    "CAPIVARA: continuum-subtracted native line windows; strict complete features; zero spatial weight.")
  bundle$provenance$emission_line_configuration <- emission$emission_line_features
  write_bundle()
  saveRDS(list(cluster_map=m,feature_configuration=emission$emission_line_features,
               backend=emission$backend_info),file.path(out,"emission_regions.rds"))
  status$emission_segmentation <- list(state="EXPLORATORY",regions=length(unique(m[!is.na(m)])))
},error=function(e) {status$emission_segmentation <<- list(state="FAILED",message=conditionMessage(e))})
write_status()
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
# Only compact web-ready products in the handoff. Preserve local research artifacts.
web <- file.path(out,"web")
dir.create(web)
stopifnot(file.copy(file.path(out,c("views.json","status.json","sessionInfo.txt")),web),
          file.copy(file.path(out,"full_flux","atlas.json"),file.path(web,"full_flux.json")),
          file.copy(file.path(out,"spectral_shape","atlas.json"),file.path(web,"spectral_shape.json")),
          file.copy(file.path(pre,"preflight.json"),web))
print(status)
