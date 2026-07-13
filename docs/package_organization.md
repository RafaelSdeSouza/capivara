# Capivara Package Organization

Capivara 2.0 is one user-facing R package with focused modules. This preserves
one installation path and one documentation site while keeping the public API
organised by scientific task.

## Core `capivara`

Purpose: segmentation, regional spectral summaries, and native kinematic
analysis for IFU cubes.

Public API:

- `segment()` for exact Ward segmentation on small/medium cubes.
- `segment_large()` for sparse-Ward segmentation on realistic IFU cubes.
- `build_starlet_mask()` and `build_adaptive_support()` for support masks.
- `summarize_cluster_spectra()` for segment spectra used by fitting packages.
- `reconstruct_cluster_cube()` and `reconstruct_flux_preserving_cube()` for
  bin-level products.
- Lightweight plotting helpers.
- `run_kinematic_analysis()` for the full native velocity-map, segmentation,
  disc, and bisymmetric workflow.
- `fit_disc_model()` and `fit_bisymmetric_model()` for expert modelling.

Experimental semantic morphology should not be the headline API. Core can keep
neutral structure-score tools, but production bar modelling belongs in the
kinematics module, where it uses velocity maps, disc geometry, residuals, and
the bisymmetric model consistently.

`spectropath` is a required dependency of the kinematics module. It remains a
small independent implementation of the path-signature mathematics, but users
do not need to install or configure it separately.

## Documentation Recommendation

Use `pkgdown` for package documentation on GitHub Pages. It understands R
packages, roxygen references, function indexes, articles, and examples directly.
Use Quarto for the broader project/paper website if desired, but keep package
reference documentation in pkgdown.

Sphinx is excellent for Python/C++ documentation, but for this R package family
it would add machinery without improving the user experience.
