# Segment an IFU cube using emission-line kinematics

Builds native emission-line maps, then clusters either the kinematic
feature cube or the path-signature feature cube. It deliberately does
not run the ordinary full-spectrum Capivara segmentation: use
\[segment()\] or \[segment_large()\] when spectral segmentation is the
scientific objective.

## Usage

``` r
segment_kinematics(
  cube_path,
  redshift = NA_real_,
  emission_line = "halpha",
  segmentation_mode = c("kinematic", "path_signature"),
  output_dir = NULL,
  object_id = NULL,
  knn_k = 50,
  n_segments = 25,
  n_path_segments = 45,
  starlet_scales = "2:5",
  include_coarse_starlet = FALSE,
  support_mode = c("starlet", "line_flux"),
  line_flux_sigma = 3,
  show_plots = interactive()
)
```

## Arguments

- cube_path:

  Path to a MaNGA LOGCUBE.

- redshift:

  Optional redshift. If \`NA\`, uses FITS header then local MaNGA
  metadata.

- emission_line:

  Emission line alias, such as \`"halpha"\` or \`"oiii5007"\`.

- segmentation_mode:

  \`"kinematic"\` clusters flux, velocity, dispersion, and profile-shape
  maps. \`"path_signature"\` additionally computes the path-signature
  feature segmentation.

- output_dir:

  Directory for saved products. Defaults beside the cube.

- object_id:

  Optional output identifier.

- knn_k:

  kNN graph size for sparse Ward clustering.

- n_segments:

  Number of kinematic-aware segments.

- n_path_segments:

  Number of path-signature segments.

- starlet_scales:

  Starlet support scales.

- include_coarse_starlet:

  Include the coarse starlet plane in support.

- support_mode:

  Foreground support for the kinematic maps. \`"starlet"\` keeps the
  white-light starlet footprint; \`"line_flux"\` trims that footprint
  with a robust emission-line-flux threshold before fitting.

- line_flux_sigma:

  Border-noise threshold, in robust sigma units, when \`support_mode =
  "line_flux"\`.

- show_plots:

  Print the compact kinematic panel.

## Value

A \`capivara_kinematic_segmentation\` object containing native maps,
support, kinematic segmentation, optional path segmentation, and paths.
