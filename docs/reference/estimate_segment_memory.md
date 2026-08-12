# Estimate segmentation memory requirements

Estimate the RAM pressure from the exact Ward backend before allocating
the all-pairs distance vector. This is a lower bound: clustering and R
object overhead require additional memory.

## Usage

``` r
estimate_segment_memory(
  input,
  valid_mode = c("signal", "finite"),
  knn_k = 20,
  overhead_factor = 2
)
```

## Arguments

- input:

  A FITS-like object with `imDat`, a raw 3-D cube array, or an integer
  giving the number of valid pixels.

- valid_mode:

  Valid-pixel rule used when `input` is a cube: `"signal"` matches
  [`segment`](https://rafaelsdesouza.github.io/capivara/reference/segment.md),
  while `"finite"` requires every spectral channel to be finite.

- knn_k:

  Number of nearest neighbours used for the sparse graph estimate.

- overhead_factor:

  Multiplicative factor used for a conservative exact Ward RAM estimate
  beyond the condensed distance vector alone.

## Value

A one-row data frame with valid-pixel count and approximate RAM
estimates in GB.
