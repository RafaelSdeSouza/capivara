# Bisymmetric Bar Modelling

The bisymmetric module is for a galaxy with independent evidence of a
bar. It extends the axisymmetric disc model with the two-sided streaming
terms, following the velocity-field formulation used in the
barred-galaxy literature. It is not a generic substitute for the disc
model.

## Before fitting

Capivara first estimates a photometric bar-angle prior from the central
white-light elongation, then deprojects it relative to the kinematic
disc major axis. It deliberately does not substitute the disc position
angle: doing so silently assumes an aligned bar and gives a misleading
model for most galaxies. Inspect the magenta prior in the component
figure before interpreting the bisymmetric terms. A measured in-plane
angle can always override the automatic estimate with `bar_phi_deg`.

## Fit the bar module

``` r
library(capivara)

bar_result <- run_kinematic_analysis(
  cube_path = "/data/known_barred_galaxy.fits",
  redshift = 0.03,
  emission_line = "halpha",
  segmentation_mode = "kinematic",
  model = "bisymmetric_bar",
  knn_k = 50,
  n_segments = 25,
  model_control = list(
    disc_inc_deg = 38,
    robust_fit = TRUE
  ),
  show_plots = TRUE
)

plot(bar_result, which = "components")
```

The `components` plot is intentionally unavailable for an axisymmetric
run: the circular, tangential second-order, and radial second-order
components only have a physical interpretation within the bisymmetric
model.

## Optional bar-support prior

The model normally uses the full valid kinematic footprint. When a
carefully measured bar footprint is needed as a weak geometric prior,
add:

``` r
model_control = list(
  bar_phi_deg = 41, # optional measured in-plane override
  use_bar_support_mask = TRUE,
  bar_support_width_deg = 25
)
```

This option should be compared against the full-footprint model, not
used to force a bar-like residual pattern.

For the neutral default workflow, start with [Kinematic
analysis](https://rafaelsdesouza.github.io/capivara/articles/kinematic-analysis.md).
A single demonstration is not evidence for a bar or outflow detection;
interpretation requires the model diagnostics and independent
observational context.
