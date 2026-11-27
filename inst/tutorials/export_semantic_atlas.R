#!/usr/bin/env Rscript
# Export a compact, interactive atlas from run_prepared_cube.R results.
# Usage: Rscript export_semantic_atlas.R analysis_directory
# Requires jsonlite; no Python, pPXF, or external fitting engine.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L) stop("Usage: Rscript export_semantic_atlas.R analysis_directory")
if (!requireNamespace("jsonlite",quietly=TRUE))
  stop("Install jsonlite with install.packages('jsonlite').")
dir <- normalizePath(args[1L],mustWork=TRUE)
outputs <- file.path(dir,c("atlas.json","region_features.csv","region_adjacency.csv"))
if (any(file.exists(outputs))) stop("Atlas output exists; choose a fresh analysis directory.")
label <- readRDS(file.path(dir,"labels_native.rds"))
m <- label$cluster_map
p <- readRDS(file.path(dir,"provenance.rds"))
t <- read.csv(file.path(dir,"regional_spectra.csv"),check.names=FALSE)
size <- read.csv(file.path(dir,"regions.csv"))
stopifnot(is.matrix(m),is.numeric(m),
          all(is.na(m) | (is.finite(m)&m>0&m==floor(m))))
ids <- sort(unique(as.integer(m[!is.na(m)])))
stopifnot(length(ids)>0L,!anyDuplicated(size$region),setequal(ids,size$region),
          setequal(ids,t$region))
wave <- sort(unique(t$wavelength_angstrom))
stopifnot(length(wave)>=3L,all(is.finite(wave)),all(diff(wave)>0))
spec <- matrix(NA_real_,length(ids),length(wave))
counts <- spec
sizes <- vapply(ids,function(id) sum(m==id,na.rm=TRUE),integer(1))
stopifnot(all(sizes==size$n_spaxels[match(ids,size$region)]))
for (i in seq_along(ids)) {
  r <- t[t$region==ids[i],,drop=FALSE]
  r <- r[order(r$wavelength_angstrom),,drop=FALSE]
  stopifnot(identical(as.numeric(r$wavelength_angstrom),as.numeric(wave)),
            all(r$n_region_spaxels==sizes[i]),
            all(is.finite(r$n_valid_spaxels)),
            all(r$n_valid_spaxels>=0 & r$n_valid_spaxels<=sizes[i]),
            all(r$n_valid_spaxels==floor(r$n_valid_spaxels)),
            all(is.na(r$fixed_aperture_sum) |
                (is.finite(r$fixed_aperture_sum)&r$n_valid_spaxels==sizes[i])))
  spec[i,] <- r$fixed_aperture_sum
  counts[i,] <- r$n_valid_spaxels
}
# Native four-neighbour contacts; a diagonal touch is not adjacency.
contact <- function(a,b) {
  a <- as.vector(a); b <- as.vector(b)
  ok <- !is.na(a)&!is.na(b)&a!=b
  data.frame(region_a=pmin(a[ok],b[ok]),region_b=pmax(a[ok],b[ok]))
}
nr <- nrow(m); nc <- ncol(m)
edges <- data.frame(region_a=integer(),region_b=integer())
if (nr>1L) edges <- rbind(edges,contact(m[-nr,,drop=FALSE],m[-1L,,drop=FALSE]))
if (nc>1L) edges <- rbind(edges,contact(m[,-nc,drop=FALSE],m[,-1L,drop=FALSE]))
if (nrow(edges)) {
  edges$shared_edges <- 1L
  edges <- aggregate(shared_edges~region_a+region_b,edges,sum)
} else edges$shared_edges <- integer()
edges$common_fraction <- rep(NA_real_,nrow(edges))
edges$cosine_distance <- rep(NA_real_,nrow(edges))
# Equal-native-channel cosine distance on the common fixed-aperture grid.
# This is a descriptive shape comparison, not a hypothesis test or anomaly
# probability. It is NOT the MaNGA/Sandra segmentation amplitude functional.
for (i in seq_len(nrow(edges))) {
  a <- spec[match(edges$region_a[i],ids),]
  b <- spec[match(edges$region_b[i],ids),]
  ok <- is.finite(a)&is.finite(b)
  edges$common_fraction[i] <- mean(ok)
  if (sum(ok)<3L || mean(ok)<.7) next
  a <- a[ok]; b <- b[ok]
  ma <- max(abs(a)); mb <- max(abs(b))
  if (ma<=0 || mb<=0) next
  a <- a/ma; b <- b/mb
  edges$cosine_distance[i] <- max(0,min(2,1-sum(a*b)/sqrt(sum(a*a)*sum(b*b))))
}
finite_median <- function(x) if(any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
features <- data.frame(region=ids,n_spaxels=sizes,
   complete_fraction=rowMeans(is.finite(spec)),
   median_flux_density_per_spaxel=vapply(seq_along(ids),
     function(i) finite_median(spec[i,]/sizes[i]),numeric(1)),
   n_neighbours=integer(length(ids)),max_neighbour_shape_distance=NA_real_)
for (i in seq_along(ids)) {
  use <- edges$region_a==ids[i] | edges$region_b==ids[i]
  features$n_neighbours[i] <- sum(use)
  x <- edges$cosine_distance[use]
  if(any(is.finite(x))) features$max_neighbour_shape_distance[i] <- max(x[is.finite(x)])
}
# Explicit lists keep even singleton arrays as JSON arrays.
rows <- lapply(seq_along(ids),function(i) {
  x <- as.list(features[i,,drop=FALSE])
  x$spectrum <- as.list(as.numeric(spec[i,]))
  x$valid_counts <- as.list(as.integer(counts[i,]))
  x
})
source <- if (identical(p$object_id,"SYNTHETIC-DEMO") &&
   identical(p$input$source,"fixed-seed simulated cube, not observational data"))
  "Synthetic demonstration" else "User-provided cube (not independently verified)"
representation <- p$representation$type
if(is.null(representation)) representation <- "not declared"
payload <- list(schema="capivara-semantic-atlas-v1",object=as.character(p$object_id),
  source=source,nx=nr,ny=nc,representation=as.character(representation),
  wavelength_frame=as.character(p$wavelength_frame),
  wavelength_medium=as.character(p$wavelength_medium),
  flux_unit=as.character(p$flux_unit),redshift=p$redshift,
  requested_k=label$requested_k,actual_k=length(ids),
  wavelength=as.list(as.numeric(wave)),
  labels=lapply(seq_len(nc),function(j) as.list(as.integer(m[,j]))),
  regions=rows,
  edges=lapply(seq_len(nrow(edges)),function(i)as.list(edges[i,,drop=FALSE])),
  provenance=list(backend=as.character(label$backend),
    coordinates="Native pixel grid; no celestial WCS claim",
    spectrum="Fixed-aperture native sum; incomplete channels withheld",
    shape="Cosine distance on jointly complete channels, equal channel weights; at least 3 channels and 70% shared coverage",
    uncertainty="Not propagated; no calibrated errors or anomaly significance",
    interpretation="Descriptive region features only; no physical fits or astrophysical labels",
    input_files=as.list(as.character(tools::md5sum(file.path(dir,
      c("labels_native.rds","regional_spectra.csv","regions.csv","provenance.rds")))))))
jsonlite::write_json(payload,outputs[1],auto_unbox=TRUE,na="null",digits=NA)
write.csv(features,outputs[2],row.names=FALSE,na="")
write.csv(edges,outputs[3],row.names=FALSE,na="")
message("Open the Galaxy Atlas and choose ",outputs[1],".")
