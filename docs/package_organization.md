# Capivara Package Organization

Capivara should be documented as a small family of R packages rather than one
large mixed namespace.

## Core `capivara`

Purpose: segmentation and flux-preserving spectral summaries for IFU cubes.

Public API:

- `segment()` for exact Ward segmentation on small/medium cubes.
- `segment_large()` for sparse-Ward segmentation on realistic IFU cubes.
- `build_starlet_mask()` and `build_adaptive_support()` for support masks.
- `summarize_cluster_spectra()` for segment spectra used by fitting packages.
- `reconstruct_cluster_cube()` and `reconstruct_flux_preserving_cube()` for
  bin-level products.
- Lightweight plotting helpers.

Experimental semantic morphology should not be the headline API. Core can keep
neutral structure-score tools, but production bar modelling belongs in the
kinematics/bar companion package.

## `capivaraKinematics`

Purpose: velocity-map products, disc models, and bisymmetric/bar modelling using
Capivara segmentation products.

Primary entry point:

- `run_manga_bar_model()`

Advanced entry points:

- `estimate_disc_geometry()`
- `fit_axisymmetric_piecewise_model()`
- `fit_bisymmetric_model()`
- `plot_capivara_kinematics()`
- `plot_capivara_component_decomposition()`

## Documentation Recommendation

Use `pkgdown` for package documentation on GitHub Pages. It understands R
packages, roxygen references, function indexes, articles, and examples directly.
Use Quarto for the broader project/paper website if desired, but keep package
reference documentation in pkgdown.

Sphinx is excellent for Python/C++ documentation, but for this R package family
it would add machinery without improving the user experience.
