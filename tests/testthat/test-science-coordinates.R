ellipse_fixture <- function(pa, nr = 61L, nc = 75L) {
  x <- col(matrix(0, nr, nc)) - (nc + 1)/2
  y <- row(matrix(0, nr, nc)) - (nr + 1)/2
  a <- pa*pi/180
  major <- x*sin(a) + y*cos(a)
  minor <- x*cos(a) - y*sin(a)
  mask <- (major/22)^2 + (minor/6)^2 <= 1
  list(mask = mask, weight = exp(-0.5*((major/15)^2+(minor/4)^2)))
}

test_that("horizontal, vertical and oblique structure PAs share one frame", {
  for (pa in c(0, 90, 27, 135)) {
    f <- ellipse_fixture(pa)
    centre <- (dim(f$mask)+1)/2
    a <- .structure_fit_component(f$mask, f$weight, centre)
    b <- .structure_weighted_ellipse(f$mask, f$weight)
    expect_lt(.structure_angle_diff_deg(a$pa_image_deg, pa), 1)
    expect_lt(.structure_angle_diff_deg(b$pa_image_deg, pa), 1)
    expect_true(a$axis_ratio <= 1 && b$axis_ratio <= 1)
    expect_false("pa_deg" %in% names(a))
    rot <- .capivara_display_matrix(f$mask, "rot90_cw")
    wt <- .capivara_display_matrix(f$weight, "rot90_cw")
    ar <- .structure_fit_component(rot, wt, (dim(rot)+1)/2)
    br <- .structure_weighted_ellipse(rot, wt)
    expect_lt(.structure_angle_diff_deg(ar$pa_image_deg, a$pa_image_deg+90), 1e-8)
    expect_lt(.structure_angle_diff_deg(br$pa_image_deg, b$pa_image_deg+90), 1e-8)
    expect_equal(sum(rot), sum(f$mask))
  }
})

test_that("image, mask and point display transforms round-trip on rectangles", {
  m <- matrix(seq_len(15), 3, 5)
  points <- expand.grid(y = 1:3, x = 1:5)
  inverse <- c(identity="identity", transpose="transpose", flip_x="flip_x",
               flip_y="flip_y", rot90_cw="rot90_ccw", rot90_ccw="rot90_cw", rot180="rot180")
  for (mode in names(inverse)) {
    rotated <- .capivara_display_matrix(m, mode)
    p <- .capivara_display_points(points$x, points$y, dim(m), mode)
    expect_equal(rotated[cbind(p$y,p$x)], m[cbind(points$y,points$x)])
    back <- .capivara_display_points(p$x,p$y,dim(rotated),inverse[[mode]])
    expect_equal(back$x, points$x)
    expect_equal(back$y, points$y)
    expect_equal(.capivara_display_matrix(rotated,inverse[[mode]]),m)
  }
})

test_that("bar deprojection and projected vectors are inverse operations", {
  s <- data.frame(x=c(0,1),y=c(0,1),valid=TRUE,velocity=c(0,1))
  for (convention in c("nirvana","capivara_legacy")) {
    g <- estimate_disc_geometry(s,list(x0=0,y0=0,vsys=0,pa_image_deg=23,
                                       inc_deg=50,coordinate_convention=convention))
    for (pa in c(0,90,37,136)) {
      phi <- .capivara_image_to_disc_angle(pa,g)
      result <- list(spaxels=data.frame(valid=TRUE,R=1:10),geometry=g,
                     bar_geometry=list(phi_bar_disc_rad=phi*pi/180))
      axis <- .capivara_bar_axis_line(result)
      angle <- atan2(diff(range(axis$x))*sign(tail(axis$x,1)-axis$x[1]),
                     diff(range(axis$y))*sign(tail(axis$y,1)-axis$y[1]))*180/pi
      expect_lt(.structure_angle_diff_deg(pa,angle),1e-8)
    }
  }
})

test_that("bar-labelled fallback uses the deprojection theta sign", {
  g <- list(x0=0,y0=0,pa_image_rad=0,inc_rad=pi/4,coordinate_convention="nirvana")
  q <- seq(-5,5,length.out=11)
  d <- deproject_coordinates(-q*sin(pi/6)*cos(g$inc_rad),q*cos(pi/6),g)
  d$valid <- TRUE
  d$seg_class <- "bar"
  b <- estimate_bar_geometry(d,geometry=g)
  expect_equal(b$phi_bar_disc_deg,-30,tolerance=1e-8)
})

test_that("sky angle requires and respects WCS handedness", {
  expect_error(.capivara_image_to_sky_angle(30,1,1,NULL),"WCS")
  wcs <- function(x,y)cbind(ra=359.999+x*1e-4,dec=y*1e-4)
  expect_equal(.capivara_image_to_sky_angle(90,1,1,wcs),90,tolerance=1e-4)
  reflected <- function(x,y)cbind(ra=20-x*1e-4,dec=y*1e-4)
  expect_equal(.capivara_image_to_sky_angle(45,1,1,reflected),135,tolerance=1e-3)
})

test_that("automatic white-light bar axis reports minor/major and image PA", {
  f <- ellipse_fixture(27)
  g <- list(x0=(ncol(f$weight)+1)/2,y0=(nrow(f$weight)+1)/2,
            pa_image_rad=0,inc_rad=0,coordinate_convention="nirvana")
  b <- .capivara_white_light_bar_axis(f$weight,g,min_pixels=5)
  expect_true(is.list(b))
  expect_gt(b$photometric_axis_ratio,0)
  expect_lte(b$photometric_axis_ratio,1)
  expect_lt(.structure_angle_diff_deg(b$pa_image_deg,27),2)
  expect_equal(b$phi_bar_disc_deg,.capivara_image_to_disc_angle(b$pa_image_deg,g))
})

test_that("clipped light plateaus cannot move automatic centres, PAs or masks on rotation", {
  light <- ellipse_fixture(27,31,41)$weight
  support <- list(collapsed=light,reconstruction=light,mask=light>0.005)
  rotated_support <- lapply(support,function(m).capivara_display_matrix(m,"rot90_cw"))
  a <- score_structures(NULL,support=support)
  b <- score_structures(NULL,support=rotated_support)
  expect_gt(sum(a$maps$collapsed==1),1L)
  expect_equal(unname(a$center),c(16,21))
  expect_equal(unname(b$center),c(21,16))
  for(n in setdiff(names(a$maps),c("orientation_cos2","orientation_sin2"))) {
    expect_equal(.capivara_display_matrix(a$maps[[n]],"rot90_cw"),b$maps[[n]],tolerance=1e-10)
  }
  da <- detect_bar(scores=a,radius_method="candidate",min_area=8L)
  db <- detect_bar(scores=b,radius_method="candidate",min_area=8L)
  expect_gt(sum(da$bar_mask),0L)
  expect_equal(.capivara_display_matrix(da$bar_mask,"rot90_cw"),db$bar_mask)
  expect_lt(.structure_angle_diff_deg(db$diagnostics$pa_image_deg,
                                     da$diagnostics$pa_image_deg+90),1e-8)
  # Genuine equal maxima are averaged, not selected by row/column scan order.
  support$collapsed[16,22] <- support$collapsed[16,21]
  c <- score_structures(NULL,support=support)
  expect_equal(unname(c$center),c(16,21.5))
})

test_that("vetted geometry is converted once and masks are not display-rotated", {
  g <- list(x0=5,y0=4,pa_image_rad=0,inc_rad=pi/4,coordinate_convention="nirvana")
  mask <- matrix(FALSE,7,9); mask[3:5,2:8] <- TRUE
  b <- list(vetted=TRUE,source="independent imaging",pa_image_deg=45,bar_mask=mask)
  out <- .capivara_vetted_bar_geometry(b,g,dim(mask))
  expect_equal(out$bar_mask,mask)
  expect_equal(out$phi_bar_disc_deg,.capivara_image_to_disc_angle(45,g))
  expect_true(is.na(out$pa_sky_deg))
  b$phi_bar_disc_deg <- out$phi_bar_disc_deg + 10
  expect_error(.capivara_vetted_bar_geometry(b,g,dim(mask)),"disagree")
  b$phi_bar_disc_deg <- NULL; b$vetted <- FALSE
  expect_error(.capivara_vetted_bar_geometry(b,g,dim(mask)),"vetted")
})
