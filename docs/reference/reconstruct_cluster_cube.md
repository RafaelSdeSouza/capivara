# Reconstruct a Model Cube from Cluster Representative Spectra

This function builds a model cube from a segmentation result by
assigning one representative spectrum to every pixel in each cluster. It
is useful for visualization, denoising, and downstream spectral analysis
where a cluster-level spectral template is desired.

## Usage

``` r
reconstruct_cluster_cube(
  cluster_result,
  template = c("median", "mean", "weighted_mean"),
  var_cube = NULL,
  variance_inflation = 1,
  preserve_mask = TRUE,
  fill_mode = c("na", "zero"),
  return_residual = TRUE
)
```

## Arguments

- cluster_result:

  A list returned by a segmentation function.

- template:

  Representative spectrum to assign to each cluster. `"median"` is
  robust, `"mean"` is the arithmetic cluster mean, and `"weighted_mean"`
  uses inverse-variance weighting.

- var_cube:

  Optional variance cube. Required when `template = "weighted_mean"`.

- variance_inflation:

  Multiplicative factor applied to propagated variances when `var_cube`
  is supplied.

- preserve_mask:

  Logical; if `TRUE`, channels that were non-finite in the original cube
  remain masked in the reconstructed cube.

- fill_mode:

  Either `"na"` or `"zero"` for unassigned pixels and, when
  `preserve_mask = TRUE`, for masked spectral channels.

- return_residual:

  Logical; if `TRUE`, also return the residual cube.

## Value

A list with the reconstructed cube, optional residual cube, the template
spectra, a cluster summary object, and a global flux comparison between
the original and reconstructed cubes.
