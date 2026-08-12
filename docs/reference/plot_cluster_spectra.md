# Visualize Scaled Spectra and Median Profiles for Clustered IFU Data

This function creates a faceted plot that displays the individual pixel
spectra (after scaling) and the corresponding median spectrum for each
segmented cluster obtained from an IFU data cube. Each facet corresponds
to a unique cluster, allowing for easy comparison of spectral profiles
across different spatial regions.

## Usage

``` r
plot_cluster_spectra(
  cluster_result,
  scale_fn = median_scale,
  palette = "magma"
)
```

## Arguments

- cluster_result:

  A list produced by a segmentation function (e.g.
  [`segment`](https://rafaelsdesouza.github.io/capivara/reference/segment.md))
  containing at least the following elements:

  - `cluster_map`: A matrix mapping each spatial pixel to a cluster.

  - `original_cube`: The original data cube (typically a FITS object).

  - `axDat`: Axis metadata (including wavelength information).

- scale_fn:

  A function that scales each individual spectrum (a numeric vector).
  The default is `median_scale` which should be defined elsewhere in
  your package.

- palette:

  Either a single viridis palette name or a vector of colors used to
  color the per-cluster median spectra.

## Value

A `ggplot2` object representing a faceted plot of the median profile for
each cluster.

## Details

The function performs the following steps:

1.  Extracts the spatial cluster map, original data cube, and wavelength
    information from `cluster_result`.

2.  Converts the cube to a 2D matrix (assuming `cube_to_matrix` is
    available).

3.  Creates a data frame mapping each valid spatial pixel to its
    cluster.

4.  Extracts each pixel's spectrum and applies the scaling function
    `scale_fn` to it.

5.  Computes the median spectrum for each cluster at every wavelength.

6.  Constructs a faceted ggplot of the median spectrum for each cluster.

## See also

[`segment`](https://rafaelsdesouza.github.io/capivara/reference/segment.md),
`cube_to_matrix()`, [`axVec`](https://rdrr.io/pkg/FITSio/man/axVec.html)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Assuming 'input_cube' is a FITS cube and 'segment' clusters the data:
  cluster_result <- segment(input_cube, Ncomp = 5)

  # Generate the spectra plot
  spectra_plot <- plot_cluster_spectra(cluster_result, scale_fn = median_scale, palette = "magma")
  print(spectra_plot)
} # }
```
