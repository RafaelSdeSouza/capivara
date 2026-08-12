# Run the explicit bisymmetric-bar module for a MaNGA cube

This convenience wrapper is intentionally bar-specific. By default it
derives a photometric in-plane bar-angle prior from the white-light
image; supplying \`bar_phi_deg\` overrides that estimate. It never
substitutes the disc position angle for a bar angle. For an ordinary
rotation model, use \[run_kinematic_analysis()\] with its axisymmetric
default.

## Usage

``` r
run_manga_bar_model(
  cube_path,
  redshift = NA_real_,
  emission_line = "halpha",
  segmentation_mode = c("kinematic", "path_signature"),
  bar_phi_deg = NA_real_,
  output_dir = NULL,
  object_id = NULL,
  knn_k = 50,
  n_segments = 25,
  n_path_segments = 45,
  support_mode = c("starlet", "line_flux"),
  line_flux_sigma = 3,
  model_control = list(),
  show_plots = interactive()
)
```

## Arguments

- cube_path:

  Path to an IFU FITS cube.

- redshift:

  Redshift. Leave as \`NA\` for a MaNGA cube with local metadata.

- emission_line:

  Emission-line alias, such as \`"halpha"\` or \`"oiii5007"\`.

- segmentation_mode:

  \`"kinematic"\` for line-map features or \`"path_signature"\` to also
  compute path-signature regions.

- bar_phi_deg:

  Optional manual in-plane bar angle in degrees relative to the disc
  major axis. Leave as \`NA\` to estimate it from white light.

- output_dir:

  Directory for saved products. Defaults beside the cube.

- object_id:

  Optional output identifier.

- knn_k:

  kNN graph size for kinematic clustering.

- n_segments:

  Number of kinematic-aware segments.

- n_path_segments:

  Number of path-signature segments.

- support_mode:

  Foreground support for the kinematic maps. \`"starlet"\` keeps the
  white-light starlet footprint; \`"line_flux"\` trims it with a robust
  emission-line-flux threshold.

- line_flux_sigma:

  Border-noise threshold, in robust sigma units, for \`"line_flux"\`
  support.

- model_control:

  Named list of model controls. For a bar model, \`bar_phi_deg\` is an
  optional manual in-plane prior; when omitted, Capivara estimates it
  from the inner white-light elongation. See \[kinematic_models()\] for
  available modules.

- show_plots:

  Print figures as they are generated.

## Value

A \`capivara_kinematic_result\` with the bisymmetric model and its
component decomposition.
