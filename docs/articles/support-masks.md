# Support Masks

A support mask identifies the spatial elements eligible for
segmentation. It does not assign region labels. Capivara can use an
explicit missing-data footprint, a white-light starlet support, or an
adaptive multi-band support.

## Starlet support

Both segmentation backends accept the same starlet controls:

``` r
seg <- segment_large(
  input = cube,
  Ncomp = 25,
  use_starlet_mask = TRUE,
  support_method = "starlet",
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  starlet_mode = "soft",
  positive_only = TRUE,
  mask_mode = "na",
  knn_k = 100
)
```

The cube is collapsed to a white-light image, decomposed across starlet
scales, and reconstructed from the selected scales. `include_coarse`
controls the coarse plane; `denoise_k` is expressed in MAD units; and
`positive_only` excludes negative reconstructed values. These settings
determine the support.

Plot the support before segmentation:

``` r
support_result <- build_starlet_mask(
  cube,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  mode = "soft",
  positive_only = TRUE
)

support_result$mask
```

## Adaptive support

The adaptive builder can combine per-band evidence and optional
transforms.

``` r
adaptive <- build_adaptive_support(
  cube,
  transform = "asinh",
  sky_method = "border",
  border_fraction = 0.1,
  z_threshold = 3,
  min_band_persistence = 2
)

seg <- segment(
  cube,
  Ncomp = 15,
  use_starlet_mask = TRUE,
  support_method = "adaptive",
  support_args = list(
    transform = "asinh",
    sky_method = "border",
    z_threshold = 3,
    min_band_persistence = 2
  )
)
```

The transform, sky estimator, persistence, and thresholds determine the
adaptive support. Excluded spaxels remain unassigned.

## Missing and background values

With `mask_mode = "na"`, excluded cube values and map locations are
missing. A FITS export may use zero for the background; the science cube
remains unchanged.
