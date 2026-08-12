# Regional Spectra

The segmentation map is categorical. Summing the spaxels in each region
gives the flux-preserving spectrum used for fitting.

## Regional summaries

``` r
products <- summarize_cluster_spectra(
  seg,
  var_cube = variance_cube,
  variance_inflation = 1
)

products$cluster_ids
products$n_spaxels
products$sum_spectra
products$sum_variance
```

For region (k), channel $`\lambda`$, and assigned spaxels $`p`$,
Capivara computes

``` math

F_{k,\lambda}^{\mathrm{sum}} = \sum_{p \in k} F_{p,\lambda}.
```

When a variance cube is supplied, the regional variance is the sum of
valid input variances multiplied by `variance_inflation`. Capivara does
not estimate spatial covariance from the flux cube.

Means and medians describe spectral shape; only the summed spectra
preserve the total regional flux.

## Representative cubes

``` r
representative <- reconstruct_cluster_cube(
  seg,
  template = "median",
  preserve_mask = TRUE,
  return_residual = TRUE
)

flux_preserving <- reconstruct_flux_preserving_cube(
  seg,
  fill_mode = "na",
  return_residual = TRUE
)
```

[`reconstruct_cluster_cube()`](https://rafaelsdesouza.com.br/capivara/reference/reconstruct_cluster_cube.md)
fills a region with a representative spectrum.
[`reconstruct_flux_preserving_cube()`](https://rafaelsdesouza.com.br/capivara/reference/reconstruct_flux_preserving_cube.md)
distributes each summed regional spectrum over its assigned spaxels so
summing the reconstruction over the region recovers the regional total.

## Map–table correspondence

Region identifier `k` in `seg$cluster_map` corresponds to `k` in
`products$cluster_ids` and to the matching row of the returned spectral
matrices. Region identifiers remain categorical even when stored as
integers.
