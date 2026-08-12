# FITS and DS9 Export

[`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
and
[`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
return a two-dimensional categorical map in `seg$cluster_map`. Its
dimensions match the input cube’s spatial footprint.

## DS9 label image

``` r
cluster_map_fits <- seg$cluster_map

# File convention only:
#   0  = outside support or unassigned background
#   >0 = categorical Capivara region identifier
cluster_map_fits[is.na(cluster_map_fits)] <- 0L
storage.mode(cluster_map_fits) <- "integer"

FITSio::writeFITSim(
  cluster_map_fits,
  file = "capivara_segmentation_map.fits",
  c1 = "Capivara map: 0=background; positive integers=region identifiers"
)
```

The zero conversion makes background handling explicit for FITS viewers.
It does not introduce a new Capivara region, and the positive integers
do not form an ordered physical scale.

## Preserve WCS

``` r
spatial_ax <- NA
if (!is.null(seg$axDat) && is.data.frame(seg$axDat) && nrow(seg$axDat) >= 2) {
  spatial_ax <- seg$axDat[1:2, , drop = FALSE]
}

FITSio::writeFITSim(
  cluster_map_fits,
  file = "capivara_segmentation_map_wcs.fits",
  axDat = spatial_ax,
  header = seg$header,
  c1 = "Capivara map: 0=background; positive integers=region identifiers"
)
```

Verify the output orientation and WCS against the original cube in the
target viewer. FITS libraries can expose axis ordering differently from
plotting code, so visual agreement alone is not a substitute for
checking coordinate metadata.

## Regional spectra

``` r
products <- summarize_cluster_spectra(seg)

write.csv(
  data.frame(region = products$cluster_ids, products$sum_spectra),
  "capivara_sum_spectra.csv",
  row.names = FALSE
)
```

The region column links every table row to the same positive identifier
in the FITS map. Preserve the wavelength coordinate, input cube
identifier, support settings, backend parameters, and output paths
alongside these files.
