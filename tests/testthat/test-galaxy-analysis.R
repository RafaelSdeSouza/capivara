test_that("high-level analysis preserves fixed apertures, measurements and provenance", {
  w <- seq(4800, 6800, length.out=16)
  flux <- array(outer(c(1, 2, 3, 4, 5, 6), seq(1, 2, length.out=16)),
                c(2L,3L,length(w)))
  valid <- array(TRUE, dim(flux))
  valid[1,1,4] <- FALSE
  d <- list(flux=flux,variance=array(.01,dim(flux)),valid=valid,
    host_mask=matrix(TRUE,2L,3L),wavelength=w,wavelength_frame="rest",
    wavelength_medium="vacuum",flux_unit="test flux density",redshift=0,
    object_id="API-FIXTURE",provenance=list(source="synthetic test"))
  rep <- full_flux_representation()
  fit <- analyse_galaxy(d, 3L, rep)
  expect_s3_class(fit,"capivara_galaxy_analysis")
  direct <- segment_spectra(list(imDat=flux,wavelength=w,
    wavelength_frame="rest",wavelength_medium="vacuum"),
    3L,rep,variance=d$variance,sample_validity=valid,redshift=0)
  expect_identical(fit$segmentation$cluster_map,direct$cluster_map)
  expect_identical(fit$segmentation$original_cube$imDat,flux)
  expect_equal(sum(fit$regions$n_spaxels),6)
  tab <- fit$spectra
  bad <- tab$n_valid_spaxels != tab$n_region_spaxels
  expect_true(any(bad))
  expect_true(all(is.na(tab$fixed_aperture_sum[bad])))
  expect_true(all(is.na(tab$diagonal_variance_only[bad])))
  expect_equal(tab$fixed_aperture_sum[!bad],tab$quality_valid_sum[!bad])
  expect_identical(fit$provenance$physical_fits,"NOT_RUN")
  expect_identical(fit$provenance$sky_wcs_export,"NOT_CERTIFIED")
  skip_if_not_installed("jsonlite")
  out <- tempfile("atlas-api-")
  on.exit(unlink(out,recursive=TRUE),add=TRUE)
  path <- explore_galaxy(fit, output_dir=out, open=FALSE)
  expect_true(file.exists(path))
  expect_true(all(file.exists(file.path(out,c("atlas.json","region_features.csv",
    "region_adjacency.csv","segmentation.rds","labels_native.rds",
    "regions.csv","regional_spectra.csv","provenance.rds","segmentation.png",
    "sessionInfo.txt")))))
  j <- jsonlite::read_json(file.path(out,"atlas.json"),simplifyVector=FALSE)
  expect_identical(j$schema,"capivara-semantic-atlas-v1")
  expect_equal(j$actual_k,nrow(fit$regions))
  expect_equal(sum(vapply(j$regions,function(r) as.numeric(r$n_spaxels),numeric(1))),6)
  expect_true(any(vapply(j$regions,function(r) any(vapply(r$spectrum,is.null,logical(1))),logical(1))))
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  expect_match(html,"load(JSON.parse(new TextDecoder()",fixed=TRUE)
  expect_error(export_atlas(out),"Atlas output exists")
  expect_error(explore_galaxy(fit, output_dir=out, open=FALSE),
               "existing results are not overwritten")
})
