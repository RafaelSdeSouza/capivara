# Large-cube Capivara segmentation

This backend keeps Ward's merge criterion but restricts candidate merges
to a k-nearest-neighbor graph. It is intended for large cubes where
[`segment`](https://rafaelsdesouza.com.br/capivara/reference/segment.md)
becomes memory-limited because exact Ward needs an all-pairs distance
object.

## Usage

``` r
segment_large(
  input,
  Ncomp = 15,
  redshift = 0,
  scale_fn = median_scale,
  target_snr = NULL,
  var_cube = NULL,
  k_values = NULL,
  wavelength_range = NULL,
  feature_wavelength_range = NULL,
  snr_stat = c("integrated", "median_per_wavelength"),
  variance_inflation = 1,
  use_starlet_mask = FALSE,
  support_method = c("starlet", "adaptive"),
  support_args = list(),
  collapse_fn = collapse_white_light,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  starlet_mode = c("soft", "hard"),
  positive_only = TRUE,
  mask_mode = c("na", "zero"),
  knn_k = 40,
  auto_k = FALSE,
  max_k = NULL,
  feature_scale = c("none", "robust_col"),
  spatial_weight = 0,
  mask = NULL,
  valid_mode = c("sagui", "signal", "finite"),
  return_details = FALSE,
  verbose = FALSE
)
```

## Arguments

- input:

  A FITS-like object with \`imDat\`, or a raw 3D array.

- Ncomp:

  Integer, the number of clusters to form. Defaults to \`15\`.

- redshift:

  Kept for API compatibility with \`segment()\`.

- scale_fn:

  Optional row-wise spectral scaling function. Defaults to
  `median_scale()`, matching
  [`segment`](https://rafaelsdesouza.com.br/capivara/reference/segment.md).

- target_snr:

  Optional minimum accepted SNR per cluster. When supplied, Capivara
  chooses the largest number of clusters whose minimum cluster SNR
  remains above this threshold.

- var_cube:

  Optional variance cube matching the input cube. Used only when
  `target_snr` is supplied.

- k_values:

  Optional candidate cluster counts tested when `target_snr` is
  supplied. For this scalable backend, the default grid is capped at 50
  clusters.

- wavelength_range:

  Optional wavelength interval used to compute SNR when `target_snr` is
  supplied.

- feature_wavelength_range:

  Optional wavelength interval used to select the spectral channels used
  for clustering. The returned `original_cube` remains the full input
  cube so downstream summed spectra are still flux-preserving across the
  full spectral axis.

- snr_stat:

  Either integrated SNR or median per-wavelength SNR when `target_snr`
  is supplied.

- variance_inflation:

  Multiplicative factor applied to propagated variances when
  `target_snr` is supplied.

- use_starlet_mask:

  Logical; if `TRUE`, build a Sagui-style support mask before
  clustering.

- support_method:

  Foreground support builder used when `use_starlet_mask = TRUE`.
  Options are `"starlet"` and `"adaptive"`.

- support_args:

  Optional named list passed to the selected support builder. For
  `"adaptive"`, arguments are passed to
  [`build_adaptive_support`](https://rafaelsdesouza.com.br/capivara/reference/build_adaptive_support.md).

- collapse_fn:

  Function used to collapse the cube to white light when
  `use_starlet_mask = TRUE` and `support_method = "starlet"`.

- starlet_J:

  Number of starlet scales when `use_starlet_mask = TRUE`.

- starlet_scales:

  Integer vector of scales kept in the reconstruction when
  `use_starlet_mask = TRUE`.

- include_coarse:

  Logical; include the coarse starlet plane when
  `use_starlet_mask = TRUE`.

- denoise_k:

  Optional denoising threshold in MAD units when
  `use_starlet_mask = TRUE`.

- starlet_mode:

  Thresholding mode for the starlet reconstruction when
  `use_starlet_mask = TRUE`.

- positive_only:

  Logical; keep only positive reconstructed values when
  `use_starlet_mask = TRUE`.

- mask_mode:

  Either `"na"` or `"zero"` for masked spaxels when
  `use_starlet_mask = TRUE`.

- knn_k:

  Number of nearest neighbours for the sparse Ward graph.

- auto_k:

  If TRUE, increase \`knn_k\` when the graph is too disconnected.

- max_k:

  Maximum k allowed when \`auto_k = TRUE\`.

- feature_scale:

  Optional column-wise feature scaling after row scaling.

- spatial_weight:

  Optional weight for appending normalized x/y coordinates.

- mask:

  Optional logical spatial mask with dimensions n_row x n_col.

- valid_mode:

  Valid-pixel rule. `"sagui"` applies the stricter finite-fraction,
  energy, and row-MAD screen used by Sagui's sparse-Ward backend;
  `"signal"` matches the exact
  [`segment`](https://rafaelsdesouza.com.br/capivara/reference/segment.md)
  support; `"finite"` requires all channels to be finite.

- return_details:

  Return features, labels, and diagnostics.

- verbose:

  Print progress messages.

## Value

A Capivara-style segmentation object with the same core fields as
[`segment`](https://rafaelsdesouza.com.br/capivara/reference/segment.md):
`cluster_map`, `header`, `axDat`, `cluster_snr`, and `original_cube`.
Scalable-backend diagnostics are stored under `backend_info`.
