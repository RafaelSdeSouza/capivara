# Cluster a 2D Representation of a Data Cube

This function processes a data cube (such as an IFU cube) by flattening
it into rows (spatial pixels) and columns (spectral variables), scaling
each row, computing pairwise distances using Capivara's internal
distance helper, and performing hierarchical clustering. The resulting
clusters are rearranged back into a 2D grid consistent with the original
spatial dimensions. The function also retains the original data cube for
reference and post-processing.

## Usage

``` r
segment(
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
  mask_mode = c("na", "zero")
)
```

## Arguments

- input:

  A FITS object representing the input data cube. Typically, this is an
  IFU data cube.

- Ncomp:

  Integer, the number of clusters to form. Defaults to \`15\`.

- redshift:

  Numeric redshift placeholder kept for API compatibility.

- scale_fn:

  A function used to scale each row of the 2D representation of the data
  cube. Defaults to [`scale`](https://rdrr.io/r/base/scale.html). If you
  have a custom scaling function, pass it here.

- target_snr:

  Optional minimum accepted SNR per cluster. When supplied, Capivara
  chooses the largest number of clusters whose minimum cluster SNR
  remains above this threshold.

- var_cube:

  Optional variance cube matching the input cube. Used only when
  `target_snr` is supplied.

- k_values:

  Optional candidate cluster counts tested when `target_snr` is
  supplied.

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

## Value

A list containing:

- cluster_map:

  A `n_rows x n_cols` matrix of cluster assignments (integers), where
  each element corresponds to a spatial pixel in the original cube
  layout.

- header:

  The header metadata from the input FITS file.

- axDat:

  Axis information (e.g., spatial and spectral axes) from the input FITS
  file.

- cluster_snr:

  A numeric vector containing the signal-to-noise ratio (SNR) for each
  cluster.

- original_cube:

  The original FITS data cube as input to the function, for reference
  and post-processing.

- starlet_info:

  When `use_starlet_mask = TRUE`, a list containing the spatial mask,
  white-light image, starlet decomposition, reconstruction, and masked
  cube used before clustering.

This process is often used in IFU data analysis, where clustering is
applied to grouped spectral profiles of spatial pixels to identify
regions with similar characteristics.

## Details

Missing spectral channels are handled automatically: non-finite values
produced during row-wise scaling are replaced with zero before distances
are computed. This keeps the standard exact workflow usable for masked
cubes without a separate public entry point.

Steps performed by the function:

1.  Reads the input FITS data cube.

2.  Converts the cube into a 2D matrix (spatial pixels x spectral
    variables).

3.  Scales the data row-wise using `scale_fn`.

4.  Computes pairwise distances between rows using an internal distance
    helper.

5.  Performs hierarchical clustering using Ward's D2 method via
    [`hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html).

6.  Cuts the dendrogram into `Ncomp` clusters and reshapes the results
    into a 2D cluster map, or, when `target_snr` is supplied, chooses
    the largest cut whose minimum cluster SNR remains above the
    requested threshold.

7.  Calculates the signal-to-noise ratio (SNR) for each cluster.

## See also

[`segment_large`](https://rafaelsdesouza.com.br/capivara/reference/segment_large.md),
[`build_starlet_mask`](https://rafaelsdesouza.com.br/capivara/reference/build_starlet_mask.md),
[`hclust`](https://rdrr.io/pkg/fastcluster/man/hclust.html),
[`cutree`](https://rdrr.io/r/stats/cutree.html)

## Examples

``` r
input_cube <- list(imDat = array(runif(5 * 5 * 12), dim = c(5, 5, 12)))
clustering_result <- segment(input = input_cube, Ncomp = 5)
cluster_map <- clustering_result$cluster_map
original_cube <- clustering_result$original_cube
```
