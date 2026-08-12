# Summarize Cluster Spectra with Optional Variance Propagation

This is the recommended post-segmentation summary layer for Capivara. It
keeps the segmentation step separate from downstream spectral products
and can return median spectra, arithmetic mean spectra, summed spectra,
and inverse-variance-weighted mean spectra.

## Usage

``` r
summarize_cluster_spectra(
  cluster_result,
  var_cube = NULL,
  variance_inflation = 1
)
```

## Arguments

- cluster_result:

  A list returned by a segmentation function.

- var_cube:

  Optional variance cube matching the dimensions of the original flux
  cube. It can be a FITS-like list with `imDat` or a raw 3-D array.

- variance_inflation:

  Multiplicative factor applied to propagated variances. Use this to
  account for covariance if needed.

## Value

A list with wavelength coordinates, cluster ids, cluster sizes,
per-wavelength finite counts, median spectra, mean spectra, summed
spectra, and, when `var_cube` is supplied, propagated sum variances and
inverse-variance-weighted means.
