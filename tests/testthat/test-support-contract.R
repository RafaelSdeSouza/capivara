quality_fixture <- function() {
  flux <- array(1, c(4, 5, 10))
  variance <- array(1, dim(flux))
  dq <- array(0L, dim(flux))
  flux[1, 1, 3] <- NA_real_
  variance[2, 2, 2:4] <- 0
  dq[3, 3, 8:10] <- 1024L
  list(flux = flux, variance = variance, dq = dq, wavelength = 5001:5010)
}

test_that("quality support retains wavelength-dependent native conditions", {
  x <- quality_fixture()
  strict <- build_quality_support(x$flux, x$wavelength, 1,
    variance=x$variance, dq=x$dq, wavelength_interval=c(5002,5009),
    wavelength_frame="observed", interval_wavelength_frame="observed")
  tolerant <- build_quality_support(x$flux, x$wavelength, .6,
    variance=x$variance, dq=x$dq, wavelength_interval=c(5002,5009),
    wavelength_frame="observed", interval_wavelength_frame="observed")
  expect_false(strict$quality_mask[1,1])
  expect_false(strict$quality_mask[2,2])
  expect_false(strict$quality_mask[3,3])
  expect_true(all(tolerant$quality_mask))
  expect_equal(strict$finite_fraction[1,1], 7/8)
  expect_equal(strict$variance_fraction[2,2], 5/8)
  expect_equal(strict$donotuse_fraction[3,3], 2/8)
  expect_equal(strict$evaluated_wavelength_limits,c(5002,5009))
})

test_that("rest-frame quality selection records exact native channels", {
  x <- quality_fixture()
  rest_interval <- c(5002,5009)/1.02
  q <- build_quality_support(x$flux, x$wavelength, .9, variance=x$variance,
    wavelength_interval=rest_interval, wavelength_frame="observed",
    interval_wavelength_frame="rest", redshift=.02)
  expect_equal(q$channel_index, 2:9)
  expect_equal(q$native_wavelength_limits, c(5002,5009))
  expect_equal(q$interval_wavelength_frame,"rest")
  expect_error(build_quality_support(x$flux,x$wavelength,.9,variance=x$variance,
    wavelength_interval=rest_interval,wavelength_frame="observed",
    interval_wavelength_frame="rest"),"redshift")
})

test_that("voxel-quality summaries do not impose a universal spaxel rejection", {
  x <- quality_fixture()
  q <- build_quality_support(x$flux, x$wavelength, 1,
    variance=x$variance, dq=x$dq, wavelength_interval=c(5002,5009),
    wavelength_frame="observed", interval_wavelength_frame="observed")
  requested <- matrix(TRUE, 4, 5)
  s <- build_capivara_support(q, requested, requested,
    analysis_mask=requested, construction_method="voxel validity fixture",
    source="package regression")
  expect_true(all(s$quality_mask))
  expect_identical(s$quality_diagnostic_mask, q$quality_mask)
  expect_false(s$quality_diagnostic_mask[1,1])
  expect_true(s$analysis_mask[1,1])
})

test_that("support preserves disconnected hosts, contaminants and ambiguity", {
  q <- matrix(TRUE, 9, 9)
  detection <- matrix(FALSE,9,9)
  detection[2:4,2:4] <- TRUE       # primary
  detection[7,7] <- TRUE           # disconnected associated clump
  detection[2,8] <- TRUE           # contaminant
  detection[8,2] <- TRUE           # unresolved candidate
  host <- matrix(FALSE,9,9); host[2:4,2:4] <- TRUE; host[7,7] <- TRUE
  ambiguous <- matrix(FALSE,9,9); ambiguous[8,2] <- TRUE
  s <- build_capivara_support(q,detection,host_mask=host,
    ambiguous_mask=ambiguous,construction_method="external labelled fixture",
    source="external")
  expect_true(s$analysis_mask[7,7])
  expect_false(s$analysis_mask[2,8])
  expect_true(s$ambiguous_mask[8,2])
  expect_equal(nrow(s$component_table),4)
  expect_equal(s$qc$morphology_operations,"none")
  expect_false(any(s$analysis_mask[5:6,5:6]))
})

test_that("support IDs and serialization are stable", {
  q <- matrix(TRUE,3,4); d <- q; h <- q
  s <- build_capivara_support(q,d,h,construction_method="user mask",
    source="external",provenance=list(input_sha256="abc"))
  copy <- unserialize(serialize(s,NULL,version=3))
  expect_identical(copy,s)
  expect_identical(.capivara_support_id(copy),s$support_id)
  s2 <- s; s2$analysis_mask[1] <- FALSE
  expect_false(identical(.capivara_support_id(s2),s$support_id))
})

test_that("analysis support enters both segmentation engines before clustering", {
  set.seed(41)
  cube <- array(10+rnorm(6*6*12),c(6,6,12))
  q <- matrix(TRUE,6,6); d <- matrix(FALSE,6,6); d[2:5,2:5] <- TRUE
  s <- build_capivara_support(q,d,d,construction_method="external fixture",source="external")
  expect_warning(a <- segment(cube,Ncomp=3,support=s), "historical spectral")
  expect_warning(b <- segment_large(cube,Ncomp=3,knn_k=8,support=s,valid_mode="signal"),
                 "historical sparse spectral")
  for (z in list(a,b)) {
    expect_true(all(is.na(z$cluster_map[!s$analysis_mask])))
    expect_true(all(is.finite(z$cluster_map[s$analysis_mask])))
    expect_identical(z$support_id,s$support_id)
    expect_identical(z$support$analysis_mask,s$analysis_mask)
  }
})

test_that("unresolved or conflicting support paths cannot be clustered", {
  cube <- array(1,c(3,3,5)); q <- matrix(TRUE,3,3); d <- q
  unresolved <- build_capivara_support(q,d,construction_method="detector only",source="capivara")
  expect_true(all(unresolved$ambiguous_mask))
  expect_error(segment_large(cube,Ncomp=2,support=unresolved),"empty")
  resolved <- build_capivara_support(q,d,d,construction_method="resolved",source="external")
  expect_error(segment_large(cube,Ncomp=2,support=resolved,mask=q),"explicit support path")
})

test_that("pre-clustering support is not post-hoc deletion", {
  # Invalid edge spectra are indistinguishable from one valid host population.
  # Deleting their historical region therefore deletes valid host pixels too.
  cube <- array(NA_real_,c(5,4,4))
  type_a <- c(1,2,4,8); type_b <- c(8,4,2,1)
  for (i in 1:5) for (j in 1:4) cube[i,j,] <- if (j<=2) type_a else type_b
  q <- matrix(TRUE,5,4); q[c(1,5),] <- FALSE
  d <- matrix(TRUE,5,4)
  s <- build_capivara_support(q,d,host_mask=d,analysis_mask=q,
    construction_method="quality fixture",source="external")
  historical <- segment(cube,Ncomp=2)
  expect_warning(aware <- segment(cube,Ncomp=2,support=s), "historical spectral")
  bad_labels <- unique(historical$cluster_map[!q])
  posthoc <- historical$cluster_map
  posthoc[posthoc %in% bad_labels] <- NA_integer_
  expect_lt(sum(is.finite(posthoc[q])),sum(q))
  expect_equal(sum(is.finite(aware$cluster_map)),sum(q))
})

test_that("MaNGA-like edge-quality ring receives no label", {
  set.seed(4); cube <- array(1+rnorm(11*11*9,sd=.02),c(11,11,9))
  rr <- sqrt((row(matrix(0,11,11))-6)^2+(col(matrix(0,11,11))-6)^2)
  q <- rr<=4; detection <- rr<=5
  support <- build_capivara_support(q,detection,host_mask=detection,
    analysis_mask=q,construction_method="MaNGA edge regression",source="external")
  expect_warning(z <- segment_large(cube,Ncomp=4,knn_k=8,support=support,
                                    valid_mode="signal"), "historical sparse spectral")
  expect_true(all(is.na(z$cluster_map[rr>4])))
  expect_true(all(is.finite(z$cluster_map[rr<=4])))
})
