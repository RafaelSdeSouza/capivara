# Changelog

## capivara 0.3.0

Released: 2026-05-05

### Added

- [`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
  as the scalable sparse-Ward backend for large cubes where exact Ward
  is RAM-limited.
- [`estimate_segment_memory()`](https://rafaelsdesouza.github.io/capivara/reference/estimate_segment_memory.md)
  to estimate the exact Ward all-pairs distance-vector RAM requirement
  before allocating it.
- `research/benchmarks_segment_large.R` for RAM-aware sparse-Ward
  benchmarking. It skips
  [`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
  automatically when the exact distance vector would exceed the
  configured threshold.

### Changed

- The public segmentation API is now limited to
  [`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
  for exact Ward and
  [`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
  for scalable sparse Ward.
- Removed experimental duplicate segmentation entry points from the
  exported API.

## capivara 0.2.0

Released: 2026-03-18

### Added

- [`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
  as the user-facing scalable segmentation path for cubes where exact
  all-pairs distances would exceed available RAM.
- [`build_starlet_mask()`](https://rafaelsdesouza.github.io/capivara/reference/build_starlet_mask.md)
  for Sagui-style white-light starlet masking before clustering.
- [`summarize_cluster_spectra()`](https://rafaelsdesouza.github.io/capivara/reference/summarize_cluster_spectra.md)
  for median, summed, and inverse-variance-weighted cluster spectra.
- [`choose_ncomp_by_snr()`](https://rafaelsdesouza.github.io/capivara/reference/choose_ncomp_by_snr.md)
  for variance-aware component selection from an SNR threshold.
- [`reconstruct_cluster_cube()`](https://rafaelsdesouza.github.io/capivara/reference/reconstruct_cluster_cube.md)
  and
  [`reconstruct_flux_preserving_cube()`](https://rafaelsdesouza.github.io/capivara/reference/reconstruct_flux_preserving_cube.md)
  for representative and flux-preserving model cubes.
- `research/manga/run_manga_starlet_comparison.R` for a full-frame MaNGA
  comparison workflow.

### Changed

- [`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
  now handles missing spectral channels directly in the exact workflow.
- The scalable segmentation path now uses block medoids rather than
  block averages, improving compact structures in large cubes.
- [`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
  and
  [`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
  now accept the optional white-light starlet mask directly via
  `use_starlet_mask`.
- `torch` is now optional; Capivara falls back to base R distance
  calculations when `torch` is unavailable.
- The public segmentation API is now limited to
  [`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
  and
  [`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
  to avoid duplicate entry points.
- GitHub and website documentation now describe the missing-data,
  starlet-mask, and variance-aware workflows.

### Fixed

- [`plot_cluster_spectra()`](https://rafaelsdesouza.github.io/capivara/reference/plot_cluster_spectra.md)
  now uses the intended pixel within each cluster instead of repeatedly
  reusing the first pixel.
- The Sagui comparison workflow now preserves the full cube dimensions
  rather than applying the starlet mask on a cropped cube.
- Signal and noise estimation no longer produces warnings when row sums
  are non-positive.
