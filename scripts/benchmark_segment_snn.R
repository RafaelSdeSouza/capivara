#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(capivara))

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }
  as.integer(value)
}

env_dims <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  dims <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(dims) != 3L || any(!is.finite(dims)) || any(dims < 2L)) {
    stop(sprintf("`%s` must be three comma-separated integers, e.g. 30,30,80.", name))
  }

  dims
}

env_logical <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  tolower(value) %in% c("1", "true", "yes", "y")
}

env_numeric <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    return(default)
  }

  as.numeric(value)
}

make_synthetic_cube <- function(n_row, n_col, n_wave, seed = 42) {
  set.seed(seed)

  wave <- seq(-1, 1, length.out = n_wave)
  templates <- rbind(
    8 + 5 * exp(-0.5 * ((wave + 0.45) / 0.10)^2),
    9 + 4 * exp(-0.5 * ((wave - 0.05) / 0.13)^2),
    7 + 6 * exp(-0.5 * ((wave - 0.42) / 0.09)^2),
    8 + 2 * sin(2 * pi * seq_along(wave) / n_wave)
  )

  xy <- expand.grid(row = seq_len(n_row), col = seq_len(n_col))
  region <- 1L +
    (xy$row > n_row / 2) +
    2L * (xy$col > n_col / 2)

  mat <- templates[region, , drop = FALSE]
  mat <- mat + matrix(stats::rnorm(nrow(mat) * n_wave, sd = 0.7), nrow = nrow(mat))
  mat <- pmax(mat, 0)

  array(mat, dim = c(n_row, n_col, n_wave))
}

time_call <- function(expr) {
  gc()
  elapsed <- system.time(value <- force(expr))[["elapsed"]]
  list(value = value, elapsed = elapsed)
}

sampled_pair_agreement <- function(labels_a, labels_b, max_pairs = 200000L, seed = 42) {
  ok <- is.finite(labels_a) & is.finite(labels_b)
  labels_a <- labels_a[ok]
  labels_b <- labels_b[ok]
  n <- length(labels_a)

  if (n < 2L) {
    return(NA_real_)
  }

  total_pairs <- n * (n - 1) / 2
  set.seed(seed)

  if (total_pairs <= max_pairs) {
    idx <- utils::combn(n, 2)
    return(mean((labels_a[idx[1, ]] == labels_a[idx[2, ]]) ==
      (labels_b[idx[1, ]] == labels_b[idx[2, ]])))
  }

  i <- sample.int(n, max_pairs, replace = TRUE)
  j <- sample.int(n, max_pairs, replace = TRUE)
  keep <- i != j
  mean((labels_a[i[keep]] == labels_a[j[keep]]) ==
    (labels_b[i[keep]] == labels_b[j[keep]]))
}

dims <- env_dims("CAPIVARA_BENCH_DIMS", c(30L, 30L, 80L))
Ncomp <- env_int("CAPIVARA_BENCH_NCOMP", 12L)
knn_k <- env_int("CAPIVARA_BENCH_KNN", 30L)
seed <- env_int("CAPIVARA_BENCH_SEED", 42L)
run_exact <- env_logical("CAPIVARA_BENCH_RUN_EXACT", TRUE)
max_exact_gb <- env_numeric("CAPIVARA_BENCH_MAX_EXACT_GB", 1)

cube <- make_synthetic_cube(
  n_row = dims[1],
  n_col = dims[2],
  n_wave = dims[3],
  seed = seed
)
input <- list(imDat = cube, hdr = NULL, axDat = NULL)

cat(sprintf(
  "Benchmark cube: %d x %d x %d, Ncomp=%d, knn_k=%d\n",
  dims[1], dims[2], dims[3], Ncomp, knn_k
))
memory <- estimate_segment_memory(input, knn_k = knn_k)
cat(sprintf("Valid pixels: %d\n", memory$valid_pixels))
cat(sprintf("Exact Ward distance vector: %.2f GB\n", memory$exact_distance_gb))
cat(sprintf("Exact Ward conservative estimate: %.2f GB\n", memory$exact_conservative_gb))
cat(sprintf("SNN kNN graph estimate: %.2f GB\n", memory$snn_knn_graph_gb))

snn <- time_call(segment_snn(input, Ncomp = Ncomp, knn_k = knn_k, verbose = FALSE))

required_fields <- c("cluster_map", "header", "axDat", "cluster_snr", "original_cube")
coherence_checks <- list(
  required_fields = all(required_fields %in% names(snn$value)),
  cluster_map_dim = identical(dim(snn$value$cluster_map), dims[1:2]),
  snn_cluster_count = length(unique(stats::na.omit(as.vector(snn$value$cluster_map)))) == snn$value$Ncomp
)

timing <- data.frame(
  backend = "segment_snn",
  elapsed_sec = snn$elapsed,
  clusters = length(unique(stats::na.omit(as.vector(snn$value$cluster_map)))),
  valid_pixels = sum(!is.na(snn$value$cluster_map))
)

agreement <- NA_real_
distance_gb <- memory$exact_distance_gb

if (isTRUE(run_exact) && distance_gb <= max_exact_gb) {
  exact <- time_call(segment(input, Ncomp = Ncomp))
  timing <- rbind(
    data.frame(
      backend = "segment",
      elapsed_sec = exact$elapsed,
      clusters = length(unique(stats::na.omit(as.vector(exact$value$cluster_map)))),
      valid_pixels = sum(!is.na(exact$value$cluster_map))
    ),
    timing
  )
  coherence_checks$valid_pixel_mask <- identical(
    is.na(exact$value$cluster_map),
    is.na(snn$value$cluster_map)
  )
  agreement <- sampled_pair_agreement(
    labels_a = as.vector(exact$value$cluster_map),
    labels_b = as.vector(snn$value$cluster_map),
    seed = seed
  )
} else {
  cat(sprintf(
    "Skipping exact segment(): estimated distance vector is %.2f GB, threshold is %.2f GB.\n",
    distance_gb, max_exact_gb
  ))
}

coherence <- data.frame(
  check = names(coherence_checks),
  pass = unlist(coherence_checks, use.names = FALSE)
)

print(timing, row.names = FALSE)
cat("\nCoherence checks:\n")
print(coherence, row.names = FALSE)
if (is.finite(agreement)) {
  cat(sprintf("\nSampled pair agreement: %.3f\n", agreement))
}
