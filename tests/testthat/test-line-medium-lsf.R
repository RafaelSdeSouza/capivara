test_that("registry wavelengths and explicit air conversion preserve line identity", {
  vac <- emission_lines("vacuum"); air <- emission_lines("air")
  expect_false(anyDuplicated(vac$name) > 0)
  expect_true(all(nzchar(vac$reference)))
  expect_equal(.capivara_match_emission_lines(vac$label)$name, vac$name)
  expect_equal(vac$rest_wavelength[vac$name == "halpha"], 6564.608)
  expect_equal(vac$rest_wavelength[vac$name == "oiii5007"], 5008.240)
  # SDSS Classic tabulated pairs differ by up to 0.003 Angstrom from its
  # published conversion approximation; do not treat them as exact inverses.
  classic_vac <- c(5008.239, 6564.614)
  expect_lt(max(abs(convert_wavelength_medium(classic_vac, "vacuum", "air") -
               c(5006.843, 6562.801))), .004)
  expect_equal(convert_wavelength_medium(air$rest_wavelength, "air", "vacuum"),
               vac$rest_wavelength, tolerance = 1e-12)
  expect_error(convert_wavelength_medium(c(NA, 5000), "air", "vacuum"), "finite")
  expect_error(.capivara_match_emission_lines("oii3727"), "Ambiguous")
  for (line in c("halpha", "oiii5007")) for (z in c(.02, .1253)) {
    lab <- vac$rest_wavelength[vac$name == line]
    airlab <- air$rest_wavelength[air$name == line]
    v <- c(-300, 0, 300); wave <- lab*(1+z)*(1+v/299792.458)
    coordinate <- function(rest, medium, input_medium = "vacuum")
      .systemic_line_coordinate(wave, rest, z, "observed", line, "truth", 600,
                                "systemic", medium, input_medium)
    expect_error(coordinate(airlab, "air"), "medium differ")
    expect_error(coordinate(airlab, "vacuum"), "registry")
    expect_gt(abs(299792.458*(wave[2]/(airlab*(1+z))-1)), 80)
    airwave <- convert_wavelength_medium(wave,"vacuum","air")
    aircoordinate <- .systemic_line_coordinate(airwave,airlab,z,"observed",line,"truth",600,"systemic","air","air")
    expect_equal(aircoordinate$velocity,v,tolerance=1e-8)
    p <- coordinate(lab, "vacuum")
    expect_equal(p$velocity, v, tolerance = 1e-8)
    expect_equal(p$provenance$input_wavelength_medium, "vacuum")
    expect_equal(p$provenance$wavelength_medium, "vacuum")
    expect_equal(p$provenance$observed_line_centre*(1+p$velocity/299792.458), wave, tolerance=1e-12)
  }
})

lsf_fixture <- function(generation) {
  dims <- c(2L, 3L, 41L); wave <- seq(5000, 5020, length.out=dims[3])
  sigma <- array(rep(seq(1, 1.2, length.out=dims[3]), each=6), dims)
  list(flux=list(imDat=array(1,dims),hdr=c("EXTNAME", "FLUX")),
    pre=list(imDat=sigma,hdr=c("BUNIT", "Angstrom")),
    post=list(imDat=sqrt(sigma^2+.5^2/12),hdr=c("BUNIT", "Angstrom")), wavelength=wave,
    drp_version=generation,path="fixture",native_names=if(generation=="v2_7_1")
      c(pre="PREDISP",post="DISP") else c(pre="LSFPRE",post="LSFPOST"))
}

test_that("both LSF generations retain sigma units, sampling and fitting semantics", {
  reference <- NULL
  for (generation in c("v2_7_1", "v3_1_1")) {
    fixture <- lsf_fixture(generation)
    lsf <- do.call(.canonical_manga_lsf, fixture)
    direct <- select_manga_lsf(lsf, "direct_profile", "point_sampled")
    template <- select_manga_lsf(lsf, "template_convolution", "pixel_integrated")
    expect_identical(direct$sigma_angstrom, fixture$post$imDat)
    expect_identical(template$sigma_angstrom, fixture$pre$imDat)
    expect_identical(template$wavelength, fixture$wavelength)
    expect_identical(template$provenance$selected_lsf, "pre")
    expect_identical(direct$provenance$selected_lsf, "post")
    expect_identical(direct$provenance$drp_version, generation)
    expect_false(direct$provenance$scalar_fallback)
    expect_error(select_manga_lsf(lsf, "template_convolution"), "missing")
    expect_equal(select_manga_lsf(lsf,"template_convolution","point_sampled")$sigma_angstrom,
                 direct$sigma_angstrom)
    # Independent forward model: Gaussian convolution plus detector top-hat.
    grid <- seq(-8,8,by=.5); true_sigma <- .4; inst <- template$sigma_angstrom[1,1,1]
    integrated <- function(s) diff(stats::pnorm(c(grid-.25,tail(grid,1)+.25), sd=s))
    signal <- integrated(sqrt(inst^2+true_sigma^2))
    fit <- optimize(function(s) sum((signal-integrated(sqrt(inst^2+s^2)))^2),c(0,2),tol=1e-10)
    expect_equal(fit$minimum,true_sigma,tolerance=1e-7)
    wrong <- optimize(function(s) sum((signal-integrated(sqrt(direct$sigma_angstrom[1,1,1]^2+s^2)))^2),c(0,2))
    expect_gt(abs(wrong$minimum-true_sigma), .02)
    if (is.null(reference)) reference <- direct$sigma_angstrom
    expect_equal(direct$sigma_angstrom,reference)
    bad <- fixture; bad$pre$imDat <- aperm(bad$pre$imDat,c(3,2,1))
    expect_error(do.call(.canonical_manga_lsf,bad),"dimensions/orientation")
    bad <- fixture; bad$pre$hdr <- c("BUNIT", "km/s")
    expect_error(do.call(.canonical_manga_lsf,bad),"units")
    bad <- fixture; bad$pre$imDat[] <- 0
    expect_error(do.call(.canonical_manga_lsf,bad),"no positive")
    bad <- fixture; bad$pre$imDat[1] <- 0
    expect_true(is.na(do.call(.canonical_manga_lsf,bad)$lsf_sigma_angstrom_pre[1]))
  }
})

test_that("native reader excludes derived DISP and verifies HDU identity", {
  file <- tempfile(fileext=".fits");file.create(file);on.exit(unlink(file))
  for (generation in c("v2_7_1","v3_1_1")) {
    f <- lsf_fixture(generation)
    f$flux$hdr <- c("INSTRUME","MaNGA","EXTNAME","FLUX","CTYPE3","WAVE","VERSDRP3",generation)
    f$flux$wavelength <- f$wavelength
    f$pre$hdr <- c(f$pre$hdr,"EXTNAME",f$native_names['pre'])
    f$post$hdr <- c(f$post$hdr,"EXTNAME",f$native_names['post'])
    hdus <- list(f$flux,NULL,NULL,f$post,f$pre,list(imDat=f$wavelength,hdr=c("EXTNAME","WAVE")))
    testthat::local_mocked_bindings(.capivara_read_fits=function(path,hdu)hdus[[hdu]])
    lsf <- read_manga_lsf(file)
    expect_equal(lsf$provenance$native_extension,f$native_names)
    expect_equal(lsf$lsf_sigma_angstrom_post,f$post$imDat)
    # A MEGACUBE's derived DISP has 30 line planes, not WAVE samples.
    hdus[[4]]$imDat <- array(100,c(2,3,30))
    expect_error(read_manga_lsf(file),"dimensions/orientation")
    hdus[[4]] <- f$post;hdus[[5]]$hdr <- c("EXTNAME","SPECRES")
    expect_error(read_manga_lsf(file),"sequence")
  }
})

test_that("air frame changes pass through vacuum before applying redshift", {
  vac <- frame_cube(); air <- vac
  air$wavelength <- convert_wavelength_medium(vac$wavelength,"vacuum","air")
  air$wavelength_medium <- "air"
  bounds <- c(4800,7400)
  air_bounds <- convert_wavelength_medium(bounds,"vacuum","air")
  for (z in c(.02,.1253)) {
    a <- .subset_cubedat_wavelength_range(air,air_bounds,"observed","rest",z)
    v <- .subset_cubedat_wavelength_range(vac,bounds,"observed","rest",z)
    expect_identical(a$wave_idx,v$wave_idx)
    expect_equal(c(a$provenance$selected_rest_wavelength_min,a$provenance$selected_rest_wavelength_max),
      convert_wavelength_medium(range(v$selected_wavelengths)/(1+z),"vacuum","air"),tolerance=1e-12)
  }
})
