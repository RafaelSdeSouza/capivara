# capivara 0.2.0

Released: 2026-03-18

## Added

- `segment_masked()` for missing-data-safe segmentation with the same clustering logic as `segment()`.
- `build_starlet_mask()` and `segment_starlet()` for Sagui-style white-light starlet masking before clustering.
- `summarize_cluster_spectra()` for median, summed, and inverse-variance-weighted cluster spectra.
- `choose_ncomp_by_snr()` for variance-aware component selection from an SNR threshold.
- `reconstruct_cluster_cube()` and `reconstruct_flux_preserving_cube()` for representative and flux-preserving model cubes.
- `scripts/run_manga_starlet_comparison.R` for a full-frame MaNGA comparison workflow.

## Changed

- `segment()` and `segment_masked()` now share the same internal segmentation core.
- `torch` is now optional; Capivara falls back to base R distance calculations when `torch` is unavailable.
- `segment_robust()` and `cube_cluster_with_snr()` are preserved as compatibility wrappers around the new API.
- GitHub and website documentation now describe the missing-data, starlet-mask, and variance-aware workflows.

## Fixed

- `plot_cluster_spectra()` now uses the intended pixel within each cluster instead of repeatedly reusing the first pixel.
- The Sagui comparison workflow now preserves the full cube dimensions rather than applying the starlet mask on a cropped cube.
- Signal and noise estimation no longer produces warnings when row sums are non-positive.
