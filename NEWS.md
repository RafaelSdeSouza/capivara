# capivara 0.2.0

Released: 2026-03-18

## Added

- `segment_big_cube()` as the user-facing scalable segmentation path for cubes where exact all-pairs distances would exceed available RAM.
- `build_starlet_mask()` for Sagui-style white-light starlet masking before clustering.
- `summarize_cluster_spectra()` for median, summed, and inverse-variance-weighted cluster spectra.
- `choose_ncomp_by_snr()` for variance-aware component selection from an SNR threshold.
- `reconstruct_cluster_cube()` and `reconstruct_flux_preserving_cube()` for representative and flux-preserving model cubes.
- `scripts/run_manga_starlet_comparison.R` for a full-frame MaNGA comparison workflow.

## Changed

- `segment()` now handles missing spectral channels directly in the exact workflow.
- `segment_big_cube()` now uses block medoids rather than block averages, improving compact structures in large cubes.
- `segment()` and `segment_big_cube()` now accept the optional white-light starlet mask directly.
- `torch` is now optional; Capivara falls back to base R distance calculations when `torch` is unavailable.
- The public segmentation API is now limited to `segment()` and `segment_big_cube()` to avoid duplicate entry points.
- GitHub and website documentation now describe the missing-data, starlet-mask, and variance-aware workflows.

## Fixed

- `plot_cluster_spectra()` now uses the intended pixel within each cluster instead of repeatedly reusing the first pixel.
- The Sagui comparison workflow now preserves the full cube dimensions rather than applying the starlet mask on a cropped cube.
- Signal and noise estimation no longer produces warnings when row sums are non-positive.
