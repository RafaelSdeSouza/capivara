# Build a Sagui-style starlet mask from a spectral cube

Build a Sagui-style starlet mask from a spectral cube

## Usage

``` r
build_starlet_mask(
  input,
  collapse_fn = collapse_white_light,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  mode = c("soft", "hard"),
  positive_only = TRUE
)
```

## Arguments

- input:

  A FITS-like object with `imDat`, or a 3-D cube array.

- collapse_fn:

  Function used to collapse the cube to a 2-D image.

- starlet_J:

  Number of starlet scales.

- starlet_scales:

  Integer vector of scales kept in the reconstruction.

- include_coarse:

  Logical; include the coarse starlet plane.

- denoise_k:

  Optional denoising threshold in MAD units.

- mode:

  Thresholding mode used while reconstructing selected starlet scales.

- positive_only:

  Logical; keep only positive reconstructed values.

## Value

A list with `collapsed`, `decomposition`, `reconstruction`, and `mask`.
