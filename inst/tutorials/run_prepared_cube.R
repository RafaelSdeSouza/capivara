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
analysis <- analyse_galaxy(d, n_regions = as.integer(k),
                           representation = full_flux_representation(),
                           output_dir = out)
message("Completed segmentation and regional extraction: ", nrow(analysis$regions),
        " regions. No physical fitting or calibrated anomaly detection performed.")
