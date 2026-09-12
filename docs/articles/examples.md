# Examples

These examples use real MaNGA cubes.

## MaNGA segmentation mosaic

![A landscape mosaic of real MaNGA galaxies and their Capivara spectral
segmentation maps, with each categorical region shown in a distinct
colour](../reference/figures/mosaic_segmented_sagui.png)

Real MaNGA inputs and Capivara region maps. Colours identify regions and
do not represent an ordered quantity.

The panels show the IFS cube, spatial support, and spectral regions. Use
[`summarize_cluster_spectra()`](https://rafaelsdesouza.com.br/capivara/reference/summarize_cluster_spectra.md)
to obtain the summed regional spectra.

## Exact and sparse-graph comparison

![MaNGA 8443-6102 comparison of exact and sparse-graph Capivara
segmentation with and without starlet
support](../reference/figures/manga_8443_6102_compare_current.png)

MaNGA 8443-6102 segmented with
[`segment()`](https://rafaelsdesouza.com.br/capivara/reference/segment.md)
and
[`segment_large()`](https://rafaelsdesouza.com.br/capivara/reference/segment_large.md),
with and without starlet support.

For a new observed cube:

``` r
library(capivara)
library(FITSio)

cube <- FITSio::readFITS("/path/to/input_cube.fits")

seg <- segment_large(
  input = cube,
  Ncomp = 25,
  use_starlet_mask = TRUE,
  starlet_J = 5,
  starlet_scales = 2:5,
  include_coarse = FALSE,
  denoise_k = 0,
  positive_only = TRUE,
  mask_mode = "na",
  knn_k = 100
)

regional_spectra <- summarize_cluster_spectra(seg)$sum_spectra
```

Replace the path with a local FITS cube. The [Get
Started](https://rafaelsdesouza.com.br/capivara/articles/getting-started.html)
guide provides a runnable simulation.

## Kinematic and model-specific examples

The [Kinematic
Analysis](https://rafaelsdesouza.com.br/capivara/articles/kinematic-analysis.html)
guide covers line maps and an axisymmetric disc model. The [Bar
Model](https://rafaelsdesouza.com.br/capivara/articles/bisymmetric-bar-model.html)
guide adds bisymmetric terms for a known barred galaxy.
