#!/usr/bin/env Rscript
# A reproducible starter for the existing CAPIVARA segmentation API.
# Rscript run_prepared_cube.R --demo output_directory [Ncomp]
# Rscript run_prepared_cube.R prepared_cube.rds output_directory [Ncomp]
# The RDS format is documented in the "Analyse your own data" vignette.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop("Usage: Rscript run_prepared_cube.R --demo|prepared_cube.rds output_dir [Ncomp]")
}
k <- if (length(args) == 3L) suppressWarnings(as.numeric(args[3L])) else 8
if (length(k) != 1L || !is.finite(k) || k < 1 || k != floor(k) ||
    k > .Machine$integer.max) stop("Ncomp must be a positive integer.")
out <- args[2L]
if (file.exists(out) || dir.exists(out)) stop("Choose a new output directory; existing results are not overwritten.")
suppressPackageStartupMessages(library(capivara))
if (identical(args[1L], "--demo")) {
  set.seed(42)
  wave <- seq(4800, 6800, length.out = 96)
  x <- row(matrix(0, 20, 20)) - 10.5
  y <- col(matrix(0, 20, 20)) - 10.5
  host <- (x/8)^2 + (y/6)^2 <= 1
  g <- function(l, w) exp(-0.5 * ((wave-l)/w)^2)
  a <- 1 + exp(-(x^2+y^2)/20)
  b <- exp(-((x-4)^2+(y+1)^2)/8)
  flux <- array(NA_real_, c(20L,20L,length(wave)))
  for (i in seq_along(wave)) {
    flux[,,i] <- a*(1-.2*g(5175,35)[i]) +
      b*(.3*g(5007,24)[i]+.7*g(6563,24)[i]) +
      matrix(rnorm(400,sd=.01),20,20)
    flux[,,i][!host] <- NA_real_
  }
  d <- list(flux=flux, wavelength=wave, wavelength_frame="rest",
    wavelength_medium="vacuum", flux_unit="arbitrary simulation units",
    variance=array(.01^2,dim(flux)), valid=is.finite(flux),
    host_mask=host, redshift=0, object_id="SYNTHETIC-DEMO",
    provenance=list(source="fixed-seed simulated cube, not observational data"))
} else {
  d <- readRDS(args[1L])
}
if (!is.list(d)) stop("Prepared input must be a list.")
required <- c("flux","wavelength","wavelength_frame","wavelength_medium",
              "flux_unit","variance","valid","host_mask","redshift","object_id","provenance")
if (!all(required %in% names(d))) stop("Missing prepared input fields: ",
                                      paste(setdiff(required,names(d)),collapse=", "))
dims <- dim(d$flux)
if (!is.numeric(d$flux) || length(dims)!=3L || any(dims<1L) ||
    !is.numeric(d$wavelength) || length(d$wavelength)!=dims[3L] ||
    any(!is.finite(d$wavelength)) || any(diff(d$wavelength)<=0))
  stop("Require numeric [spatial1, spatial2, wavelength] flux and increasing wavelengths in Angstrom.")
if (length(d$wavelength_frame)!=1L || !d$wavelength_frame %in% c("observed","rest") ||
    length(d$wavelength_medium)!=1L || !d$wavelength_medium %in% c("air","vacuum") ||
    !is.numeric(d$redshift) || length(d$redshift)!=1L ||
    !is.finite(d$redshift) || d$redshift<=-1)
  stop("Explicit wavelength frame, medium and systemic redshift required.")
if (!is.numeric(d$variance) || !identical(dim(d$variance),dims) ||
    !is.logical(d$valid) || !identical(dim(d$valid),dims) || anyNA(d$valid) ||
    !is.logical(d$host_mask) || !identical(dim(d$host_mask),dims[1:2]) || anyNA(d$host_mask))
  stop("Require aligned variance, logical voxel eligibility and logical host mask.")
if (length(d$object_id)!=1L || is.na(d$object_id) || !nzchar(d$object_id) ||
    length(d$flux_unit)!=1L || is.na(d$flux_unit) || !nzchar(d$flux_unit) ||
    !is.list(d$provenance)) stop("Require object ID, flux unit and provenance list.")
valid <- d$valid & is.finite(d$flux) & is.finite(d$variance) & d$variance>0
for (j in seq_len(dims[3L])) valid[,,j] <- valid[,,j] & d$host_mask
if (!any(valid)) stop("No eligible science samples.")
cube <- list(imDat=d$flux, wavelength=d$wavelength,
             wavelength_frame=d$wavelength_frame, wavelength_medium=d$wavelength_medium,
             hdr=d$hdr, axDat=d$axDat)
# Explicit native-channel full-flux geometry; never zero-fill missing samples.
# A fixed K is an exploratory hierarchy cut, not a physical population count.
seg <- segment_large(input=cube, Ncomp=as.integer(k), redshift=d$redshift,
                      var_cube=d$variance, sample_validity=valid,
                      representation=full_flux_representation())
# Extract from unchanged calibrated flux, with separate raw finite sum and
# quality-masked fixed-aperture products. Diagonal errors are not covariance.
raw <- summarize_cluster_spectra(seg)
quality_seg <- seg
quality_seg$original_cube$imDat[!valid] <- NA_real_
p <- summarize_cluster_spectra(quality_seg,var_cube=d$variance)
complete <- sweep(p$finite_counts,1L,unname(p$n_spaxels),"==")
fixed <- p$sum_spectra
fixed[!complete] <- NA_real_
diagonal_variance <- p$sum_variance
diagonal_variance[!complete] <- NA_real_
ids <- p$cluster_ids
tab <- data.frame(region=rep(ids,each=length(p$wavelength)),
                  wavelength_angstrom=rep(p$wavelength,times=length(ids)),
                  raw_finite_sum=as.vector(t(raw$sum_spectra)),
                  quality_valid_sum=as.vector(t(p$sum_spectra)),
                  fixed_aperture_sum=as.vector(t(fixed)),
                  diagonal_variance_only=as.vector(t(diagonal_variance)),
                  n_valid_spaxels=as.vector(t(p$finite_counts)),
                  n_region_spaxels=rep(unname(p$n_spaxels),each=length(p$wavelength)))
stopifnot(all(is.na(tab$fixed_aperture_sum) |
                tab$n_valid_spaxels==tab$n_region_spaxels))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if (!dir.exists(out)) stop("Could not create output directory.")
saveRDS(seg,file.path(out,"segmentation.rds"))
saveRDS(list(cluster_map=seg$cluster_map, requested_k=k,
             actual_k=length(ids), backend=seg$backend),
        file.path(out,"labels_native.rds"))
write.csv(tab,file.path(out,"regional_spectra.csv"),row.names=FALSE)
write.csv(data.frame(region=ids,n_spaxels=unname(p$n_spaxels)),
          file.path(out,"regions.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
saveRDS(list(object_id=d$object_id, flux_unit=d$flux_unit,
             wavelength_frame=d$wavelength_frame, wavelength_medium=d$wavelength_medium,
             redshift=d$redshift, input=d$provenance, requested_k=k, actual_k=length(ids),
             representation=seg$representation, backend=seg$backend,
             uncertainty="diagonal variance only; spatial/spectral covariance not propagated",
             physical_fits="NOT_RUN", sky_wcs_export="NOT_CERTIFIED"),
        file.path(out,"provenance.rds"))
# Native pixel coordinates only; no claim of celestial WCS.
map <- seg$cluster_map
grid <- expand.grid(spatial1=seq_len(nrow(map)),spatial2=seq_len(ncol(map)))
grid$region <- factor(as.vector(map))
fig <- ggplot2::ggplot(grid,ggplot2::aes(spatial1,spatial2,fill=region))+
  ggplot2::geom_raster()+ggplot2::coord_equal()+
  scale_capivara_identity(name="Region")+theme_capivara("paper")+
  ggplot2::labs(title=d$object_id,x="Native spatial axis 1 [pixel]",y="Native spatial axis 2 [pixel]")
ggplot2::ggsave(file.path(out,"segmentation.png"),fig,width=6,height=5,dpi=180)
message("Completed segmentation and regional extraction: ",length(ids)," regions. ",
        "No physical fitting or calibrated anomaly detection performed.")
