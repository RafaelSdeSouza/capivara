test_that("regional API retains native sums, covariance and failed regions", {
  skip_if_not_installed("reticulate"); skip_if_not_installed("Matrix")
  skip_if_not(reticulate::py_module_available("scipy"))
  w <- 4800:4927; dims <- c(2L,1L,length(w)); f <- array(2, dims)
  lsf <- list(wavelength = as.numeric(w), lsf_sigma_angstrom_pre = array(1.1,dims),
    lsf_sigma_angstrom_post = array(rep(c(1.2,1.5),length(w)),dims),
    provenance = list(quantity="Gaussian sigma_lambda",units="Angstrom",wavelength_medium="vacuum",
      wavelength_frame="observed",drp_version="v3_1_1",native_extension=c(pre="LSFPRE",post="LSFPOST")))
  seg <- list(cluster_map=matrix(1,2,1))
  r <- prepare_segment_spectra(seg,f,array(1,dims),lsf,array(TRUE,dims),"fixture","fixed",.03)
  expect_s3_class(r,"capivara_regional_spectra")
  a <- r$regions[[1]]
  expect_equal(as.numeric(a$native_flux),rep(4,length(w)))
  expect_equal(as.numeric(a$native_variance),rep(2,length(w)))
  expect_equal(a$selected_lsf,"post")
  expect_equal(a$segmentation_id,"fixed")
  expect_true(inherits(a$covariance,"sparseMatrix"))
  tmp <- tempfile();saveRDS(r,tmp);expect_equal(readRDS(tmp)$regions[[1]]$native_flux,a$native_flux)
  bad <- prepare_segment_spectra(seg,f,array(1,dims),lsf,array(FALSE,dims),"fixture","fixed",.03)
  expect_length(bad$regions,1)
  expect_equal(bad$regions[[1]]$qc,"INSUFFICIENT_WAVELENGTH")
  expect_error(prepare_segment_spectra(seg,f,array(1,dims),lsf,array(1,dims),"x","y",.03),"logical mask")
})
