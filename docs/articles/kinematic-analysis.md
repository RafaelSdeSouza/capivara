# Kinematic Analysis

Capivara separates three questions that are often mixed together:

1.  Which spaxels belong together in emission-line kinematics?
2.  Is an axisymmetric rotating disc an adequate comparison model?
3.  Is there evidence for a more specific perturbation, such as a bar?

The first two are useful for both barred and unbarred galaxies. A bar
model is therefore not a default setting.

## 1. Kinematic-aware segmentation

[`segment_kinematics()`](https://rafaelsdesouza.github.io/capivara/reference/segment_kinematics.md)
makes the native line maps and clusters flux, velocity, dispersion, and
profile-shape information. It does not run ordinary full-spectrum
segmentation, which answers a different question.

``` r
library(capivara)

segments <- segment_kinematics(
  cube_path = "/data/galaxy.fits",
  redshift = 0.03,
  emission_line = "halpha",
  segmentation_mode = "kinematic",
  knn_k = 50,
  n_segments = 25,
  show_plots = TRUE
)
```

Choose `segmentation_mode = "path_signature"` when the path-signature
feature representation is part of the science question. It additionally
writes that segmentation; it does not replace the native kinematic maps.

`support_mode = "starlet"` is the conservative default and retains the
white-light galaxy footprint. For emission-line work in a field with
appreciable sky, set `support_mode = "line_flux"`; Capivara then uses a
robust border-noise threshold on the measured line-flux map and keeps
its main connected component. `line_flux_sigma = 3` is a useful strict
starting point; lower it if faint extended line emission is
scientifically important.

## 2. Axisymmetric disc comparison

The default model is an axisymmetric disc. It produces an observed
velocity map, disc model, residual map, and circular-speed profile
without making any claim that the galaxy hosts a bar.

``` r
result <- run_kinematic_analysis(
  cube_path = "/data/galaxy.fits",
  redshift = 0.03,
  emission_line = "halpha",
  model = "axisymmetric",
  segmentation_mode = "kinematic",
  knn_k = 50,
  n_segments = 25,
  show_plots = TRUE
)

plot(result, which = "model")
```

For a MaNGA LOGCUBE, leave `redshift = NA_real_` to resolve it from the
local metadata/header. For another IFU cube, provide the redshift
explicitly.

## 3. Model modules

Use
[`kinematic_models()`](https://rafaelsdesouza.github.io/capivara/reference/kinematic_models.md)
to see the installed dynamical modules. Each module uses the same native
maps and kinematic segmentation, so new models such as a spiral
perturbation can be added without changing the upstream workflow.

``` r
kinematic_models()
#>             model             label requires_bar_angle components_plot
#> 1    axisymmetric Axisymmetric disc              FALSE           FALSE
#> 2 bisymmetric_bar   Bisymmetric bar               TRUE            TRUE
```

The bar-specific workflow is documented separately in
[`vignette("bisymmetric-bar-model", package = "capivara")`](https://rafaelsdesouza.github.io/capivara/articles/bisymmetric-bar-model.md).

Return to the [spectral-segmentation first
run](https://rafaelsdesouza.github.io/capivara/articles/getting-started.md),
or inspect the [result-led
examples](https://rafaelsdesouza.github.io/capivara/articles/examples.md).
Kinematic labels remain categorical; the velocity and dispersion maps
are separate measured or modelled quantities.
