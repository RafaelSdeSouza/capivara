# Build an adaptive multi-band/spatial support mask

Builds a foreground support from band-wise or wavelength-slice detection
evidence rather than from a single white-light image. This is useful for
IFU cubes where the source is persistent over many wavelength channels
but the background may vary spatially or spectrally.

## Usage

``` r
build_adaptive_support(
  input,
  bands = NULL,
  transform = c("none", "asinh", "signed_log1p", "log1p", "copula_uniform",
    "copula_gaussian"),
  sky_method = c("border", "all"),
  border_fraction = 0.1,
  z_threshold = 3,
  min_band_persistence = NULL,
  single_band_z = Inf,
  smooth_sigma = 0
)
```

## Arguments

- input:

  A FITS-like object with `imDat`, or a 3-D cube array.

- bands:

  Wavelength slices used to build the support. May be `NULL` for all
  slices, a numeric vector of indices, or a character vector matching
  cube dimnames along the third dimension.

- transform:

  Per-slice image transform before robust sky normalization. Options are
  `"none"`, `"asinh"`, `"signed_log1p"`, `"log1p"`, `"copula_uniform"`,
  and `"copula_gaussian"`.

- sky_method:

  Pixels used for robust sky/noise estimation. `"border"` uses image
  borders; `"all"` uses all finite pixels.

- border_fraction:

  Fraction of rows/columns used as image border when
  `sky_method = "border"`.

- z_threshold:

  Per-slice detection threshold in robust sigma units.

- min_band_persistence:

  Minimum number of selected slices above `z_threshold`. If `NULL`,
  defaults to two slices when possible.

- single_band_z:

  Optional single-slice rescue threshold. Use `Inf` to disable.

- smooth_sigma:

  Optional Gaussian smoothing applied to each significance map before
  thresholding. Requires `EBImage` or `imager`.

## Value

A list with `collapsed`, `reconstruction`, `mask`, `evidence`,
`evidence_count`, `z_maps`, `band_stats`, and metadata.
