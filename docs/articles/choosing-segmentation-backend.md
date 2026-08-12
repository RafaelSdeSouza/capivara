# Choosing a Segmentation Backend

Capivara exposes an exact Ward backend and a sparse-graph Ward backend.
They return compatible region maps and regional-product inputs, but they
make different computational approximations.

## Exact Ward

[`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
constructs the full pairwise distance object among valid spaxels. The
distance storage therefore grows quadratically with the number of
eligible spaxels. Estimate the lower-bound memory before a large run:

``` r
estimate_segment_memory(cube)

seg <- segment(
  input = cube,
  Ncomp = 25,
  use_starlet_mask = TRUE
)
```

Use the exact backend when the supported footprint is small enough for
the available memory and the all-pairs Ward construction is part of the
intended analysis.

## Sparse Ward

[`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
avoids storing the full all-pairs distance matrix by building a coherent
nearest-neighbour graph.

``` r
seg <- segment_large(
  input = cube,
  Ncomp = 25,
  use_starlet_mask = TRUE,
  knn_k = 100,
  auto_k = FALSE,
  verbose = TRUE
)

seg$backend_info
```

`knn_k` controls graph density. It is an analysis parameter, not a
hidden performance toggle: record it with the input cube, support
definition, random seed where applicable, and output path. `auto_k` may
increase graph density to obtain a connected graph; inspect and retain
`backend_info` when using it.

## Select spectral channels

By default, all wavelength channels drive the segmentation. Set
`feature_wavelength_range` to learn labels from a chosen interval while
keeping the full cube in the returned object:

``` r
seg_window <- segment_large(
  input = cube,
  Ncomp = 25,
  feature_wavelength_range = c(lambda_min, lambda_max),
  use_starlet_mask = TRUE,
  knn_k = 100
)

full_regional_spectra <- summarize_cluster_spectra(seg_window)$sum_spectra
```

`feature_wavelength_range` selects clustering features. In contrast,
when `target_snr` is used, `wavelength_range` selects the channels used
by the S/N screen. The two intervals answer different questions and
should not be treated as aliases.

## Choose the component count

[`choose_ncomp_by_snr()`](https://rafaelsdesouza.github.io/capivara/reference/choose_ncomp_by_snr.md)
evaluates candidate cuts and returns the largest tested count whose
minimum regional S/N remains above the requested threshold.

``` r
choice <- choose_ncomp_by_snr(
  input = cube,
  target_snr = 30,
  var_cube = variance_cube,
  k_values = c(5, 10, 15, 20),
  wavelength_range = c(lambda_min, lambda_max)
)
```

The target, variance treatment, candidate grid, and wavelength interval
are part of the scientific configuration and should accompany the
result.
