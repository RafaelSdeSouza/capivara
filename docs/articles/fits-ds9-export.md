# FITS and DS9

[`segment()`](https://rafaelsdesouza.com.br/capivara/reference/segment.md)
and
[`segment_large()`](https://rafaelsdesouza.com.br/capivara/reference/segment_large.md)
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

Zero marks the background in the FITS file. Positive integers identify
Capivara regions.

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

Check the output orientation and WCS against the original cube. FITS
libraries may expose axis ordering differently from plotting code.

## Regional spectra

``` r
products <- summarize_cluster_spectra(seg)

write.csv(
  data.frame(region = products$cluster_ids, products$sum_spectra),
  "capivara_sum_spectra.csv",
  row.names = FALSE
)
```

The `region` column links each table row to the same identifier in the
FITS map. Save the wavelength coordinate and segmentation parameters
with the files.
