# capivara 0.3.0

Released: 2026-05-05

## Added

- `segment_large()` as the scalable sparse-Ward backend for large cubes where exact Ward is RAM-limited.
- `estimate_segment_memory()` to estimate the exact Ward all-pairs distance-vector RAM requirement before allocating it.
- `scripts/benchmark_segment_large.R` for RAM-aware sparse-Ward benchmarking. It skips `segment()` automatically when the exact distance vector would exceed the configured threshold.

## Changed

- The public segmentation API is now limited to `segment()` for exact Ward and `segment_large()` for scalable sparse Ward.
- Removed experimental duplicate segmentation entry points from the exported API.

# capivara 0.2.0

Released: 2026-03-18

## Added

- `segment_large()` as the user-facing scalable segmentation path for cubes where exact all-pairs distances would exceed available RAM.
- `build_starlet_mask()` for Sagui-style white-light starlet masking before clustering.
- `summarize_cluster_spectra()` for median, summed, and inverse-variance-weighted cluster spectra.
- `choose_ncomp_by_snr()` for variance-aware component selection from an SNR threshold.
- `reconstruct_cluster_cube()` and `reconstruct_flux_preserving_cube()` for representative and flux-preserving model cubes.
- `scripts/run_manga_starlet_comparison.R` for a full-frame MaNGA comparison workflow.

## Changed

- `segment()` now handles missing spectral channels directly in the exact workflow.
- The scalable segmentation path now uses block medoids rather than block averages, improving compact structures in large cubes.
- `segment()` and `segment_large()` now accept the optional white-light starlet mask directly via `use_starlet_mask`.
- `torch` is now optional; Capivara falls back to base R distance calculations when `torch` is unavailable.
- The public segmentation API is now limited to `segment()` and `segment_large()` to avoid duplicate entry points.
- GitHub and website documentation now describe the missing-data, starlet-mask, and variance-aware workflows.

## Fixed

- `plot_cluster_spectra()` now uses the intended pixel within each cluster instead of repeatedly reusing the first pixel.
- The Sagui comparison workflow now preserves the full cube dimensions rather than applying the starlet mask on a cropped cube.
- Signal and noise estimation no longer produces warnings when row sums are non-positive.
