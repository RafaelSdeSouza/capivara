# Audit evidence, not a production test suite. Run with Rscript from any directory.
# Requires installed capivara (including its compiled routine), FITSio and graph
# dependencies. The last probe additionally needs reticulate/ppxf and the two
# local pilot inputs named below. Edit these local paths to reproduce elsewhere.
# R source is loaded from the checkout; compiled code comes from the installed
# namespace. Failures are printed per probe, so inspect output rather than using
# the script exit status as a test-suite pass/fail signal.
repo <- "/Users/rd23aag/Documents/GitHub/capivara"
e <- new.env(parent = asNamespace("capivara"))
for (f in list.files(file.path(repo, "R"), pattern = "[.]R$", full.names = TRUE)) sys.source(f, e)
show_case <- function(name, expr) {
  cat("\n", name, "\n", sep = "")
  tryCatch(print(force(expr)), error = function(err) cat("ERROR: ", conditionMessage(err), "\n", sep = ""))
}
cat("READ-ONLY SOURCE PROBES; generated ", format(Sys.time()), "\n", sep = "")
cat("R: ", R.version.string, "\n", sep = "")
show_case("Installed/source API", list(installed_version=as.character(packageVersion("capivara")), installed=formals(capivara::run_manga_bar_model), source=formals(e$run_manga_bar_model)))
seg <- list(original_cube=list(imDat=array(c(10, NA_real_), c(2,1,1))), cluster_map=matrix(1L,2,1))
show_case("Missing flux with finite variance: expected weighted mean=10, variance=1", e$summarize_cluster_spectra(seg,var_cube=array(1,c(2,1,1))))
missing_seg <- seg
missing_seg$original_cube$imDat[] <- NA_real_
show_case("All-missing channel", e$summarize_cluster_spectra(missing_seg,var_cube=array(NA_real_,c(2,1,1))))
show_case("Counterclockwise point transform on 2x3 array: expected x=1,y=1", e$.capivara_display_points(3,1,c(2,3),"rot90_ccw"))
mat <- matrix(1:6,2,3)
show_case("Matrix rotation target for original [1,3]=5", which(e$.capivara_display_matrix(mat,"rot90_ccw")==5,arr.ind=TRUE))
cap <- list(segmentation_map=matrix(1L,2,3),segment_table=data.frame(label=1L,class="disc",use_for_disc_fit=TRUE,use_for_bar_diagnostics=FALSE))
sp <- e$.capivara_build_spaxels(list(velocity=mat),cap)
flag <- matrix(c(FALSE,TRUE,FALSE,FALSE,FALSE,FALSE),2,3)
sp$wrong_flag <- as.vector(flag)
sp$correct_flag <- flag[cbind(sp$y,sp$x)]
show_case("Imputation flag placement",sp[c("x","y","wrong_flag","correct_flag")])
show_case("Starlet finite pixels, 24x24 J=5", {
  set.seed(31)
  dec <- e$starlet_mask(matrix(runif(24*24),24),J=5)
  vapply(c(dec$w,list(dec$cJ)),function(z)sum(is.finite(z)),integer(1))
})
ellipse_mask <- ((col(matrix(0,41,51))-26)/16)^2+((row(matrix(0,41,51))-21)/5)^2<=1
show_case("Horizontal ellipse PA disagreement", list(catalogue=e$.structure_fit_component(ellipse_mask,ellipse_mask+0,c(21,26))$pa_deg,bar=e$.structure_weighted_ellipse(ellipse_mask,ellipse_mask+0,c(21,26))$pa_deg))
show_case("Bar-labelled fallback and NIRVANA theta", {
  g <- list(x0=0,y0=0,pa_rad=0,inc_rad=pi/4,coordinate_convention="nirvana")
  q <- seq(-5,5,length.out=11)
  X <- q*cos(pi/6); Yd <- q*sin(pi/6)
  d <- e$deproject_coordinates(-Yd*cos(g$inc_rad),X,g)
  b <- data.frame(valid=TRUE,seg_class="bar",X=d$X,Yd=d$Yd)
  list(fallback_phi=e$estimate_bar_geometry(b)$phi_b_deg,theta_at_positive_end=d$theta[11]*180/pi)
})
workflow <- parse(file.path(repo,"inst/extdata/kinematics/native_kinematics_workflow.R"))
w <- new.env(parent=e)
for (ex in workflow) {
  if (is.call(ex)&&identical(ex[[1]],as.name("<-"))&&is.call(ex[[3]])&&identical(ex[[3]][[1]],as.name("function"))) eval(ex,w)
}
show_case("FITS channel names from actual workflow expressions", {
  w$support_starlet <- matrix(1,2,3); w$support <- matrix(2,2,3)
  w$seg <- list(cluster_map=matrix(3,2,3)); w$kin_seg <- list(cluster_map=matrix(10,2,3))
  w$kin <- list(flux=matrix(1e4,2,3),velocity=matrix(5,2,3),sigma=matrix(6,2,3),asymmetry=matrix(7,2,3),h3_proxy=matrix(8,2,3),h4_proxy=matrix(9,2,3))
  w$line <- list(slug="halpha")
  for(ex in workflow) if(is.call(ex)&&identical(ex[[1]],as.name("<-"))&&grepl("^(maps|names\\(maps\\)\\[3:8\\])$",paste(deparse(ex[[2]]),collapse=""))) eval(ex,w)
  data.frame(channel=seq_along(w$maps),name=names(w$maps),sentinel=vapply(w$maps,function(x)x[1],numeric(1)))
})
show_case("Augmented structure result stores feature channels as spectra", {
  x <- list(imDat=array(seq_len(4*5*8),c(4,5,8)))
  scores <- list(maps=list(test=matrix(.5,4,5)),support_mask=matrix(TRUE,4,5),masks=list(structure_mask=matrix(TRUE,4,5)),threshold=list())
  z <- e$segment_structures(x,scores,Ncomp=2,feature_maps="test",feature_repeats=3,knn_k=5)
  list(input_channels=8,returned_channels=dim(z$original_cube$imDat)[3],summary_channels=ncol(e$summarize_cluster_spectra(z)$sum_spectra))
})
show_case("All-missing emission spectrum eligibility", {
  wave <- seq(6500,6620,length.out=201)
  a <- array(1,c(2,3,length(wave)))
  for(i in 1:2)for(j in 1:3)a[i,j,] <- 1+(i+j)*exp(-.5*((wave-6562.8)/2)^2)
  a[1,1,] <- NA_real_
  x <- list(imDat=a,axDat=data.frame(ctype=c("x","y","WAVE"),crpix=c(1,1,1),crval=c(1,1,6500),cdelt=c(1,1,.6),len=c(2,3,201)))
  z <- e$segment_emission_lines(x,0,lines="halpha",Ncomp=2,knn_k=3,max_pixels=Inf)
  list(missing_pixel_label=z$cluster_map[1,1],returned_axDat=z$axDat,original_axDat=z$original_cube$axDat)
})
show_case("FITS header redshift reader", {
  f <- tempfile(fileext=".fits")
  FITSio::writeFITSim(array(1,c(2,3,4)),f)
  full <- FITSio::readFITS(f,hdu=1)
  short <- tryCatch(FITSio::readFITS(f,hdu=1,maxLines=1),error=function(err)conditionMessage(err))
  list(full_dimension=dim(full$imDat),maxLines1_success=is.list(short),note="Short-header fixture only; does not localize earlier real-file connection warnings")
})
show_case("Exact and sparse distance semantics", {
  for(probe_seed in seq_len(100)) {
    set.seed(probe_seed); x <- matrix(runif(60),12,5)
    a <- stats::cutree(fastcluster::hclust(stats::dist(x,"manhattan"),method="ward.D2"),3)
    b <- e$.sparse_ward_cluster_matrix(x,3,knn_k=11)$labels
    c <- stats::cutree(fastcluster::hclust(stats::dist(x,"euclidean"),method="ward.D2"),3)
    if(!identical(outer(a,a,"=="),outer(b,b,"==")))break
  }
  list(seed=probe_seed,same_partition_L1_L2=identical(outer(a,a,"=="),outer(b,b,"==")),same_partition_sparse_full_L2=identical(outer(b,b,"=="),outer(c,c,"==")))
})
show_case("Source/installed function equality", c(detector=identical(body(e$detect_bar),body(getFromNamespace("detect_bar","capivara"))),wrapper=identical(body(e$run_manga_bar_model),body(capivara::run_manga_bar_model))))
show_case("Prototype duplicate function bodies", {
  p <- new.env(parent=e)
  for(f in list.files(file.path(repo,"research/capivaraKinematics_prototype/R"),pattern="[.]R$",full.names=TRUE))sys.source(f,p)
  fs <- ls(p,all.names=TRUE);fs <- fs[vapply(fs,function(n)is.function(p[[n]])&&is.function(e[[n]]),logical(1))]
  data.frame(function_name=fs,identical_body=vapply(fs,function(n)identical(body(p[[n]]),body(e[[n]])),logical(1)))
})
show_case("Dependency provenance", lapply(c("capivara","spectropath","capivaraPPXF"),function(n){d<-packageDescription(n);d[c("Package","Version","RemoteSha")]}))
show_case("pPXF paired rebin regression on saved pilot wavelength arrays", {
  util <- reticulate::import("ppxf.ppxf_util")
  do.call(rbind,lapply(c("11004-12701","11014-3704"),function(id) {
    x <- readRDS(file.path("/Users/rd23aag/Documents/GitHub/iFUN/Capivara_Eat_Manga/results/bars_v1/stellar_populations",id,"ppxf_input.rds"))
    rest <- x$wavelength/(1+x$redshift)
    lam <- as.numeric(rest[rest>=4800 & rest<=7400])
    first <- util$log_rebin(lam,rep(1,length(lam)))
    second <- util$log_rebin(lam,rep(1,length(lam)),velscale=first[[3]])
    data.frame(observation=id,redshift=x$redshift,initial_length=length(first[[1]]),repeat_length=length(second[[1]]),paired_repeat_wave=length(second[[2]]))
  }))
})
