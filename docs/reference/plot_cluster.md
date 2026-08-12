# Plot a Cluster Map with Discrete Cluster Colors

This function visualizes the output of a clustering algorithm applied to
an IFU data cube, where each cluster is assigned a distinct fill color.

## Usage

``` r
plot_cluster(cluster_data, palette = "starry_night")
```

## Arguments

- cluster_data:

  A list containing the clustering results with the following
  components:

  - `cluster_map`: A matrix representing the spatial layout of the
    clusters.

  - `cluster_snr`: A numeric vector with the SNR values for each
    cluster.

- palette:

  Character or vector of colours used for the discrete cluster colors.
  The default `"starry_night"` uses the Van Gogh inspired blue and gold
  palette used in the companion workflows. Viridis palette names such as
  `"magma"`, `"inferno"`, `"plasma"`, `"viridis"`, and `"cividis"` are
  also supported.

## Value

A `ggplot2` object representing the cluster map.

## Details

The function converts the `cluster_map` into a long-format data frame
and colors each cluster with a discrete palette.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example cluster data
cluster_data <- list(
  cluster_map = matrix(sample(1:5, 100, replace = TRUE), nrow = 10, ncol = 10),
  cluster_snr = c(10, 20, 30, 40, 50)
)

# Plot the cluster map
plot <- plot_cluster(cluster_data, palette = "starry_night")
print(plot)
} # }
```
