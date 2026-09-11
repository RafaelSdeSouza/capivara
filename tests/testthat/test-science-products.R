test_that("regional variances follow the contributing measured flux", {
  z <- list(original_cube=list(imDat=array(c(10,NA),c(2,1,1))),
            cluster_map=matrix(1L,2,1))
  s <- summarize_cluster_spectra(z,var_cube=array(1,c(2,1,1)))
  expect_equal(unname(s$weighted_mean_spectra[1,1]),10)
  expect_equal(unname(s$weighted_mean_variance[1,1]),1)
  expect_equal(unname(s$sum_variance[1,1]),1)
  z$original_cube$imDat[] <- NA
  s <- summarize_cluster_spectra(z,var_cube=array(NA_real_,c(2,1,1)))
  expect_true(all(is.na(s$sum_spectra)))
  expect_true(all(is.na(s$sum_variance)))
  z$original_cube$imDat[] <- c(10,20)
  s <- summarize_cluster_spectra(z,var_cube=array(c(1,NA),c(2,1,1)))
  expect_equal(unname(s$sum_spectra[1,1]),30)
  expect_true(is.na(s$sum_variance[1,1]))
  expect_equal(unname(s$weighted_mean_spectra[1,1]),10)
  expect_error(summarize_cluster_spectra(z,variance_inflation=-1),"positive")
})

test_that("starlet reflection supports small fields without changing ordinary padding", {
  expect_equal(.reflect_pad_vec(1:5,2),c(2,1,1:5,5,4))
  expect_equal(.reflect_pad_vec(7,3),rep(7,7))
  set.seed(31)
  d <- starlet_mask(matrix(runif(24*19),24,19),J=5)
  expect_true(all(vapply(c(d$w,list(d$cJ)),function(x)all(is.finite(x)),logical(1))))
})

test_that("structure features never become wavelength samples", {
  x <- list(imDat=array(seq_len(4*5*8),c(4,5,8)),wavelength=seq(5000,5007))
  scores <- list(maps=list(test=matrix(.5,4,5)),support_mask=matrix(TRUE,4,5),
                 masks=list(structure_mask=matrix(TRUE,4,5)),threshold=list())
  z <- segment_structures(x,scores,Ncomp=2,feature_maps="test",feature_repeats=3,knn_k=5)
  expect_equal(z$original_cube,x)
  expect_equal(ncol(summarize_cluster_spectra(z)$sum_spectra),8)
})

test_that("emission segmentation excludes absent spectra and retains wavelength metadata", {
  wave <- seq(6500,6620,length.out=201)
  a <- array(1,c(2,3,length(wave)))
  for(i in 1:2)for(j in 1:3)a[i,j,] <- 1+(i+j)*exp(-.5*((wave-6562.8)/2)^2)
  a[1,1,] <- NA_real_
  x <- list(imDat=a,axDat=data.frame(ctype=c("x","y","WAVE"),crpix=c(1,1,1),
           crval=c(1,1,6500),cdelt=c(1,1,.6),len=c(2,3,201)))
  for (limit in c(Inf,4)) {
    z <- segment_emission_lines(x,0,lines="halpha",Ncomp=2,knn_k=3,max_pixels=limit)
    expect_true(is.na(z$cluster_map[1,1]))
    expect_equal(z$axDat,x$axDat)
  }
})

test_that("FITS connections close on success and failure", {
  path <- tempfile(fileext=".fits")
  FITSio::writeFITSim(array(1:24,c(2,3,4)),path)
  on.exit(unlink(path))
  before <- getAllConnections()
  expect_equal(dim(.capivara_read_fits(path)$imDat),c(2L,3L,4L))
  for (i in 1:4) expect_error(suppressWarnings(.capivara_read_fits(path,hdu=99)))
  expect_equal(getAllConnections(),before)
  expect_true("NAXIS" %in% .capivara_read_fits_header(path))
})

test_that("absent bar support is undefined, not a measured zero", {
  s <- data.frame(valid=TRUE,v_resid=c(1,2,3),seg_class="disc")
  d <- compute_residual_diagnostics(s)$summary
  expect_true(is.na(d$N_bar))
  expect_true(is.na(d$f_bar))
  expect_true(is.na(d$Q_kin))
  expect_equal(d$bar_diagnostic_qc,"bar_support_not_supplied")
  d <- compute_residual_diagnostics(s,bar_support_available=TRUE)$summary
  expect_equal(d$N_bar,0)
  expect_true(is.na(d$f_bar))
  expect_equal(d$bar_diagnostic_qc,"no_valid_spaxels_in_supplied_bar_support")
})

test_that("scientific models fail before I/O if geometry is missing", {
  expect_error(run_kinematic_analysis("does-not-exist.fits",show_plots=FALSE),
               "requires model_control")
  expect_error(run_manga_bar_model("does-not-exist.fits",show_plots=FALSE,
                                  model_control=list(disc_inc_deg=40)),
               "explicit phi_bar_disc_deg")
  expect_equal(.capivara_model_control(list(analysis_mode="preview"))$analysis_mode,"preview")
  expect_warning(z <- .capivara_model_control(list(disc_pa_deg=30)),"deprecated")
  expect_equal(z$disc_pa_image_deg,30)
})
