frame_cube <- function() {
  set.seed(410)
  wave <- seq(4300, 9300, by = 17)
  list(imDat = array(10 + runif(12 * length(wave)), c(3, 4, length(wave))),
       wavelength = wave, wavelength_frame = "observed", wavelength_medium = "vacuum")
}

test_that("both spectral APIs select the requested physical interval on native channels", {
  x <- frame_cube()
  for (fun in list(segment, segment_large)) for (z in c(0, .02, .1253)) {
    extra <- if (identical(fun, segment_large)) list(knn_k = 5) else list()
    call <- function(frame, cube = x, redshift = z) do.call(fun, c(list(
      input = cube, Ncomp = 2, redshift = redshift,
      feature_wavelength_range = c(4800, 7400), feature_wavelength_frame = frame), extra))
    obs <- call("observed")
    rest <- call("rest")
    io <- which(x$wavelength >= 4800 & x$wavelength <= 7400)
    ir <- which(x$wavelength >= 4800 * (1 + z) & x$wavelength <= 7400 * (1 + z))
    expect_identical(obs$feature_wavelength_index, io)
    expect_identical(rest$feature_wavelength_index, ir)
    expect_identical(rest$original_cube, x)
    p <- rest$wavelength_provenance
    expect_named(p, c("input_wavelength_frame", "requested_feature_wavelength_range",
      "requested_feature_wavelength_frame", "systemic_redshift",
      "selected_native_wavelength_min", "selected_native_wavelength_max",
      "selected_rest_wavelength_min", "selected_rest_wavelength_max",
      "selected_channel_indices", "number_of_selected_channels", "native_selection_range",
      "wavelength_medium", "coordinate_kind", "resampled", "channel_index_convention"))
    expect_equal(p$selected_native_wavelength_min, x$wavelength[min(ir)])
    expect_equal(p$selected_native_wavelength_max, x$wavelength[max(ir)])
    expect_equal(p$selected_rest_wavelength_min, x$wavelength[min(ir)] / (1 + z))
    expect_equal(p$selected_rest_wavelength_max, x$wavelength[max(ir)] / (1 + z))
    expect_equal(p$number_of_selected_channels, length(ir))
    expect_identical(p$selected_channel_indices, ir)
    expect_equal(p$systemic_redshift, z)
    expect_identical(p$input_wavelength_frame, "observed")
    expect_identical(p$requested_feature_wavelength_frame, "rest")
    expect_equal(p$requested_feature_wavelength_range, c(4800, 7400))
    expect_false(p$resampled)
    expect_equal(summarize_cluster_spectra(rest)$wavelength_provenance, p)
    # Changing coordinates to rest, without changing sampled flux, must not
    # apply systemic redshift again.
    xr <- x; xr$wavelength <- x$wavelength / (1 + z); xr$wavelength_frame <- "rest"
    rest_input <- call("rest", xr)
    expect_identical(rest_input$feature_wavelength_index, ir)
    expect_equal(rest_input$cluster_map, rest$cluster_map)
    expect_identical(call("observed", xr)$feature_wavelength_index, io)
    if (z == 0) expect_equal(rest$cluster_map, obs$cluster_map)
  }
})

test_that("ambiguous and invalid physical requests fail before clustering", {
  x <- frame_cube()
  for (fun in list(segment, segment_large)) {
    expect_error(fun(x, Ncomp = 2, redshift = .1, feature_wavelength_range = c(4800,7400)),
                 "Ambiguous historical")
    for (bad in list(NA_real_, Inf, -1, c(.02, .12), "0.02")) {
      expect_error(fun(x, Ncomp = 2, redshift = bad, feature_wavelength_range = c(4800,7400),
                       feature_wavelength_frame = "rest"), "redshift")
    }
    expect_error(fun(x, Ncomp = 2, feature_wavelength_range = c(4800,7400),
                     feature_wavelength_frame = "rest"), "redshift")
    y <- x; y$wavelength_frame <- NULL
    expect_error(fun(y, Ncomp = 2, redshift = .1, feature_wavelength_range = c(4800,7400),
                     feature_wavelength_frame = "rest"), "wavelength_frame")
    expect_error(fun(x, Ncomp = 2, wavelength_frame = "rest"), "conflicts")
    expect_error(fun(x, Ncomp = 2, feature_wavelength_frame = "observed",
                     feature_wavelength_range = c(7400,4800)), "increasing")
    expect_error(fun(x, Ncomp = 2, feature_wavelength_frame = "observed",
                     feature_wavelength_range = c(14000,18000)), "No wavelengths")
  }
  expect_error(.subset_cubedat_wavelength_range(array(1,c(2,2,10)), c(4,8),
                 "observed", "observed", 0), "Physical wavelength")
  # Inclusive exact endpoints and unchanged explicit observed selection.
  x$wavelength <- seq_len(dim(x$imDat)[3]) + 4799
  a <- .subset_cubedat_wavelength_range(x, c(4800,4810), NULL, "observed", .02)
  b <- .subset_cubedat_wavelength_range(x, c(4800,4810), NULL, "observed", .1253)
  expect_identical(a$wave_idx, 1:11)
  expect_identical(a$cubedat$imDat, b$cubedat$imDat)
})

test_that("variance and SNR selection use the same native channel indices", {
  x <- frame_cube()
  v <- array(2, dim(x$imDat))
  for (fun in list(segment, segment_large)) {
    extra <- if (identical(fun, segment_large)) list(knn_k=5) else list()
    z <- do.call(fun, c(list(input=x, target_snr=.1, var_cube=v, k_values=2,
      redshift=.1253, feature_wavelength_range=c(4800,7400),
      feature_wavelength_frame="rest"),extra))
    expect_equal(z$Ncomp,2)
    expect_true(all(is.finite(z$cluster_snr)))
    expect_equal(z$wavelength_provenance$selected_channel_indices,
      which(x$wavelength>=4800*1.1253 & x$wavelength<=7400*1.1253))
  }
  p <- segment(array(1:120,c(3,4,10)),Ncomp=2)$wavelength_provenance
  expect_equal(p$coordinate_kind,"channel_or_feature_index")
  expect_true(is.na(p$selected_native_wavelength_min))
  expect_equal(p$number_of_selected_channels,10)
})

test_that("line coordinates and conventional features have the same physical meaning at different z", {
  catalogue <- emission_lines("vacuum")
  expect_equal(catalogue$rest_wavelength[catalogue$name=="halpha"],6564.632)
  expect_equal(catalogue$rest_wavelength[catalogue$name=="oiii5007"],5008.240)
  v <- seq(-2200,2200,20)
  profile <- 1 + 10*exp(-.5*((v-140)/70)^2)
  reference <- NULL
  for (line in c("halpha","oiii5007")) for (z in c(.02,.1253)) {
    lab <- catalogue$rest_wavelength[catalogue$name==line]
    wave <- lab*(1+z)*(1+v/299792.458)
    c1 <- .systemic_line_coordinate(wave,lab,z,"observed",line,"synthetic truth",600,"systemic","vacuum")
    c2 <- .systemic_line_coordinate(wave/(1+z),lab,z,"rest",line,"synthetic truth",600,"systemic","vacuum")
    expect_equal(c1$velocity,v,tolerance=1e-8)
    expect_equal(c2$velocity,v,tolerance=1e-8)
    expect_equal(c1$provenance$observed_line_centre,lab*(1+z))
    cube <- array(rep(profile,each=9),c(3,3,length(v)))
    k <- compute_line_maps(cube,wave,matrix(TRUE,3,3),z,lab,line)
    expect_equal(k$velocity[1,1],140,tolerance=.3)
    expect_equal(k$velocity,k$velocity_systemic)
    expect_equal(k$velocity_median_centered,matrix(0,3,3),tolerance=1e-9)
    expect_gt(k$median_velocity_offset_kms,139)
    kp <- k$frame_provenance
    expect_true(all(c("line_name","line_rest_wavelength","systemic_redshift",
      "systemic_redshift_source","observed_line_centre","velocity_definition",
      "velocity_window_kms","profile_centering_mode") %in% names(kp)))
    expect_equal(kp$velocity_window_kms,c(-600,600))
    conventional <- spectropath::classical_features(cbind(c1$velocity,profile-1))
    if (is.null(reference)) reference <- conventional
    expect_equal(conventional,reference,tolerance=1e-8)
    expect_error(build_path_feature_cube(cube,wave,lab*(1+z)^2,matrix(TRUE,3,3),
       rest_wave=lab,redshift=z,line_name=line),"double application")
  }
})

test_that("local-centroid path representation is explicit and preserves systemic extraction", {
  v <- seq(-2190,2210,20)
  f <- exp(-.5*((v-133)/80)^2) - .05*exp(-.5*((v+150)/50)^2)
  for (z in c(.02,.1253)) {
    lab <- 5008.24; wave <- lab*(1+z)*(1+v/299792.458)
    cube <- array(rep(1+f,each=6),c(2,3,length(v)))
    get_path <- function(mode, x=cube) build_path_feature_cube(x,wave,lab*(1+z),matrix(TRUE,2,3),
      rest_wave=lab,redshift=z,line_name="oiii5007",profile_centering_mode=mode)
    a <- get_path("systemic"); b <- get_path("local_centroid")
    expect_equal(a$selected_channel_indices,b$selected_channel_indices)
    expect_equal(a$line_cube,b$line_cube)
    expect_equal(b$frame_provenance$profile_centering_mode,"local_centroid")
    expect_true(all(abs(b$table$represented_centroid_kms)<1e-8))
    expect_true(all(a$table$represented_centroid_kms>120))
    expect_equal(a$table$centroid_kms,b$table$centroid_kms)
    # Translation-invariant signatures cannot themselves encode the bulk centroid.
    expect_equal(a$feature_cube,b$feature_cube,tolerance=1e-9)
    cube[1,1,a$selected_channel_indices[4]] <- NA_real_
    missing <- get_path("systemic",cube)
    expect_true(all(is.na(missing$feature_cube[1,1,])))
  }
})

test_that("emission-feature API obeys frame and medium metadata in both projection branches", {
  v <- seq(-2500,2500,25); lab <- 5008.24
  for (z in c(.02,.1253)) for (limit in c(Inf,4)) {
    wave <- lab*(1+z)*(1+v/299792.458)
    cube <- array(1,c(2,3,length(v)))
    for(i in 1:2)for(j in 1:3)cube[i,j,] <- 1+(i+j)*exp(-.5*((v-(i-j)*50)/80)^2)
    x <- list(imDat=cube,wavelength=wave,wavelength_frame="observed",wavelength_medium="vacuum")
    a <- segment_emission_lines(x,z,lines="oiii5007",Ncomp=2,knn_k=3,max_pixels=limit,feature_mode="moments")
    x$wavelength <- wave/(1+z); x$wavelength_frame <- "rest"
    b <- segment_emission_lines(x,z,lines="oiii5007",Ncomp=2,knn_k=3,max_pixels=limit,feature_mode="moments")
    expect_equal(a$cluster_map,b$cluster_map)
    expect_equal(a$wavelength_provenance$selected_channel_indices,b$wavelength_provenance$selected_channel_indices)
    expect_equal(a$kinematic_provenance$oiii5007$observed_line_centre,lab*(1+z))
  }
})

test_that("model runners transport resolved redshift provenance and explicit profile controls", {
  control <- .capivara_model_control(list(wavelength_frame="rest",
    wavelength_medium="air",profile_centering_mode="local_centroid"))
  expect_equal(control$wavelength_frame,"rest")
  expect_equal(control$profile_centering_mode,"local_centroid")
  expect_error(.capivara_model_control(list(wavelength_frame="unknown")),"arg")
  # Stop at the native dispatcher to test the real model runner's resolved
  # arguments without fitting a disc or creating synthetic FITS models.
  path<-tempfile(fileext=".fits");file.create(path);on.exit(unlink(path))
  captured<-NULL
  testthat::local_mocked_bindings(
    resolve_manga_redshift=function(...)list(redshift=.1253,source="catalogue fixture",plateifu="test"),
    .capivara_run_native_kinematics=function(...) {captured<<-list(...);stop("captured dispatcher")})
  expect_error(run_kinematic_analysis(path,redshift=NA_real_,show_plots=FALSE,
    output_dir=tempdir(),model_control=list(analysis_mode="preview",profile_centering_mode="local_centroid")),
    "captured dispatcher")
  expect_equal(captured$redshift,.1253)
  expect_equal(captured$systemic_redshift_source,"catalogue fixture")
  expect_equal(captured$profile_centering_mode,"local_centroid")
})
