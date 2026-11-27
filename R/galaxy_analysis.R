#' Analyse a prepared galaxy cube
#'
#' Define spectral regions and retain raw, quality-filtered and fixed-aperture
#' spectral products. This workflow does not perform physical fitting.
#' @param prepared Prepared input list with flux, variance, valid, host_mask,
#'   wavelength, wavelength_frame, wavelength_medium, flux_unit, redshift,
#'   object_id and provenance. See the Analyse your own data vignette.
#' @param n_regions Positive integer requested hierarchy cut.
#' @param representation Explicit scientific representation declaration.
#' @return A capivara_galaxy_analysis object with segmentation, regional
#'   spectra, measurements and provenance. Incomplete fixed-aperture channels
#'   are missing; uncertainties are diagonal only.
#' @export
analyse_galaxy <- function(prepared, n_regions, representation) {
  d <- prepared
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
  seg <- segment_spectra(cube, n_regions = n_regions,
                           representation = representation, variance = d$variance,
                           sample_validity = valid, redshift = d$redshift)
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
  provenance <- list(object_id=d$object_id, flux_unit=d$flux_unit,
    wavelength_frame=d$wavelength_frame, wavelength_medium=d$wavelength_medium,
    redshift=d$redshift, input=d$provenance, requested_k=as.integer(n_regions),
    actual_k=length(ids), representation=seg$representation, backend=seg$backend,
    uncertainty="diagonal variance only; spatial/spectral covariance not propagated",
    physical_fits="NOT_RUN", sky_wcs_export="NOT_CERTIFIED")
  structure(list(segmentation=seg, raw_spectra=raw, quality_spectra=p,
    spectra=tab, regions=data.frame(region=ids,n_spaxels=unname(p$n_spaxels)),
    labels=list(cluster_map=seg$cluster_map,requested_k=as.integer(n_regions),
                actual_k=length(ids),backend=seg$backend),
    provenance=provenance), class="capivara_galaxy_analysis")
}

# Write the same seven analysis products as the prepared-cube tutorial.
.capivara_write_analysis <- function(analysis, output_dir) {
  if (!inherits(analysis, "capivara_galaxy_analysis"))
    stop("Supply a result from analyse_galaxy().",call.=FALSE)
  if (!is.character(output_dir) || length(output_dir)!=1L ||
      is.na(output_dir) || !nzchar(output_dir) ||
      file.exists(output_dir) || dir.exists(output_dir))
    stop("Choose a new output directory; existing results are not overwritten.",call.=FALSE)
  if (!dir.create(output_dir, recursive=TRUE, showWarnings=FALSE))
    stop("Could not create output directory.",call.=FALSE)
  out <- normalizePath(output_dir,mustWork=TRUE)
  saveRDS(analysis$segmentation,file.path(out,"segmentation.rds"))
  saveRDS(analysis$labels,file.path(out,"labels_native.rds"))
  utils::write.csv(analysis$spectra,file.path(out,"regional_spectra.csv"),row.names=FALSE)
  utils::write.csv(analysis$regions,file.path(out,"regions.csv"),row.names=FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()),file.path(out,"sessionInfo.txt"))
  saveRDS(analysis$provenance,file.path(out,"provenance.rds"))
  map <- analysis$labels$cluster_map
  grid <- expand.grid(spatial1=seq_len(nrow(map)),spatial2=seq_len(ncol(map)))
  grid$region <- factor(as.vector(map))
  spatial1 <- spatial2 <- region <- NULL
  fig <- ggplot2::ggplot(grid,ggplot2::aes(spatial1,spatial2,fill=region))+
    ggplot2::geom_raster()+ggplot2::coord_equal()+
    scale_capivara_identity(name="Region")+theme_capivara("paper")+
    ggplot2::labs(title=analysis$provenance$object_id,
                  x="Native spatial axis 1 [pixel]",y="Native spatial axis 2 [pixel]")
  ggplot2::ggsave(file.path(out,"segmentation.png"),fig,width=6,height=5,dpi=180)
  out
}
