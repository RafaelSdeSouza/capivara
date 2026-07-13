# Public API Audit

This note records the intended user-facing surface after the kinematics/bar
split. It is meant as a maintenance checklist before releases.

## Core `capivara`

Core should stay focused on cube segmentation and segment-level spectra.

Keep public:

- `segment()` and `segment_large()`
- `build_starlet_mask()`, `build_adaptive_support()`, and `detect_support()`
- `summarize_cluster_spectra()`
- `reconstruct_cluster_cube()` and `reconstruct_flux_preserving_cube()`
- `plot_cluster()`, `plot_cluster_spectra()`, and
  `capivara_plot_spectra_with_map()`
- Memory/SNR helpers such as `estimate_segment_memory()` and
  `choose_ncomp_by_snr()`

Keep secondary but visible for now:

- `score_structures()`, `segment_structures()`, `catalogue_structures()`,
  and `threshold_structures()`

Hide from the public API:

- `detect_bar()`
- `detect_ring()`

Those routines are exploratory semantic-structure prototypes. Production bar
modelling now belongs to `capivaraKinematics`, where it can use velocity maps,
disc geometry, residuals, and the bisymmetric model consistently.

## `capivaraKinematics`

The companion package should expose complete workflows first, then expert
building blocks.

Primary public workflow:

- `run_manga_bar_model()`

Useful teaching/batch helpers:

- `run_one_galaxy()`
- `run_batch_barred_manga()`
- `read_capivara_output()`

Model and plotting functions worth exposing:

- `fit_disc_model()`
- `fit_axisymmetric_piecewise_model()`
- `fit_bisymmetric_model()`
- `estimate_disc_geometry()`
- `estimate_bar_geometry()`
- `plot_capivara_kinematics()`
- `plot_capivara_component_decomposition()`

MaNGA convenience helpers worth exposing while MaNGA remains a first-class
example:

- `infer_manga_plateifu()`
- `resolve_manga_redshift()`
- `manga_redshift()`
- `manga_metadata()`

## Documentation Policy

Use roxygen comments for function-level documentation and pkgdown for rendered
package websites. Put tutorials in vignettes/articles once the code stabilizes,
but keep source-able R scripts under `scripts/` or `examples/` for hands-on
teaching.

Do not expose experimental functions just because they are convenient during
development. Prefer one polished public workflow plus a few honest expert
helpers over a wide namespace full of half-promises.
