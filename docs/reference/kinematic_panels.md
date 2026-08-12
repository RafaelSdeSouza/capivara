# Extract individual kinematic plot panels

The standard \`plot()\` method creates a compact multi-panel figure.
This helper returns the individual ggplot objects so they can be saved,
restyled, or assembled into a paper-specific layout without rerunning
the analysis.

## Usage

``` r
kinematic_panels(result, view = c("model", "components"))
```

## Arguments

- result:

  A result from \[run_kinematic_analysis()\] or
  \[run_manga_bar_model()\].

- view:

  \`"model"\` extracts the disc or bisymmetric-model panels.
  \`"components"\` extracts the bar-component panels and is available
  only for \`model = "bisymmetric_bar"\`.

## Value

A named list of ggplot objects.

## Details

For an axisymmetric model, \`view = "model"\` returns \`footprint\`,
\`velocity\`, \`disc_model\`, \`disc_residual\`, \`circular_speed\`, and
\`summary\`. For a bisymmetric bar model it returns the corresponding
full-model, axisymmetric, residual, and velocity-profile panels. \`view
= "components"\` returns the white-light/bar-prior and
kinematic-segmentation context panels, the individual circular and
bar-flow component maps, supporting line-property maps, and residual
comparisons.
