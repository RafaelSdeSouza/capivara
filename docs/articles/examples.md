# Examples

These panels are existing real-data outputs, not outputs of the
synthetic first-run cube. Their provenance is recorded in
`docs/website_image_provenance.csv`; unresolved generator,
configuration, and machine-readable product mappings are marked there
rather than inferred from filenames.

## MaNGA segmentation mosaic

![A landscape mosaic of real MaNGA galaxies and their Capivara spectral
segmentation maps, with each categorical region shown in a distinct
colour](../reference/figures/mosaic_segmented_sagui.png)

Real MaNGA inputs and Capivara categorical region maps. Colours
distinguish labels and do not encode an ordered quantity. See the
provenance table for the current reproduction status.

The visual sequence is: IFS cube, eligible spatial support, categorical
spectral regions, then summed regional spectra. The image directly shows
the input and region-map stages; regional spectra are obtained with
[`summarize_cluster_spectra()`](https://rafaelsdesouza.github.io/capivara/reference/summarize_cluster_spectra.md)
and are not implicitly encoded by the colours.

## Exact and sparse-graph comparison

![MaNGA 8443-6102 comparison of exact and sparse-graph Capivara
segmentation with and without starlet
support](../reference/figures/manga_8443_6102_compare_current.png)

MaNGA 8443-6102: existing comparison of
[`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
and
[`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
with and without starlet support. This panel is retained with an
explicit incomplete-provenance status until its exact generator and
configuration are verified.

Use the same scientific configuration record for a new observed cube:

``` r
library(capivara)
library(FITSio)

cube <- FITSio::readFITS("/path/to/verified-input-cube.fits")

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

The path is intentionally a user-supplied placeholder. The executable
[Getting
started](https://rafaelsdesouza.github.io/capivara/articles/getting-started.md)
guide contains the self-contained first run with no hidden file
dependency.

## Kinematic and model-specific examples

The [Kinematic
analysis](https://rafaelsdesouza.github.io/capivara/articles/kinematic-analysis.md)
guide separates line-map segmentation from an axisymmetric disc
comparison. The [Bisymmetric bar
models](https://rafaelsdesouza.github.io/capivara/articles/bisymmetric-bar-model.md)
guide applies only when independent evidence supports a bar hypothesis;
neither example changes the interpretation of ordinary spectral-region
labels.
