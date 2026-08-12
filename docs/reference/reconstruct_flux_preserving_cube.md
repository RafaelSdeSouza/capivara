# Reconstruct a Flux-preserving Model Cube from Cluster Spectra

This is a convenience wrapper around
[`reconstruct_cluster_cube`](https://rafaelsdesouza.github.io/capivara/reference/reconstruct_cluster_cube.md)
that uses the arithmetic cluster mean while preserving the original
spectral mask. In this mode, the reconstructed cube preserves the summed
flux spectrum of the segmented cube on the observed support, which makes
it a sensible default for later spectral fitting.

## Usage

``` r
reconstruct_flux_preserving_cube(
  cluster_result,
  fill_mode = c("na", "zero"),
  return_residual = TRUE
)
```

## Arguments

- cluster_result:

  A list returned by a segmentation function.

- fill_mode:

  Either `"na"` or `"zero"` for unassigned pixels and masked channels.

- return_residual:

  Logical; if `TRUE`, also return the residual cube.

## Value

The same structure returned by
[`reconstruct_cluster_cube`](https://rafaelsdesouza.github.io/capivara/reference/reconstruct_cluster_cube.md).
