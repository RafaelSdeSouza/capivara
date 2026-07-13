# capivaraKinematics

Companion R package for Capivara-native velocity maps and segmentation-aware
disc/bisymmetric modelling of IFU cubes.

This code intentionally lives outside the core Capivara package. It is the
kinematics companion in the same spirit as `capivaraPPXF`.

## Install

Install the core package first, then this companion package:

```r
install.packages("remotes")
remotes::install_github("RafaelSdeSouza/capivara")
remotes::install_local("/path/to/capivara/extensions/capivaraKinematics", dependencies = TRUE)
```

`spectropath` is a required dependency, so it is installed automatically with
`capivaraKinematics`. It provides the path-signature segmentation features.

## Native Workflow

The production path should not depend on DAP products. The native workflow is:

1. read a MaNGA LOGCUBE;
2. build a starlet galaxy mask with Capivara;
3. clean the mask with `better_galaxy_mask()` so the largest galaxy component
   is retained and enclosed holes are filled;
4. extract emission-line flux, velocity, dispersion, and asymmetry maps from
   the cube with Capivara code;
5. nearest-fill small missing emission-line holes inside the starlet footprint
   and save an imputation flag;
6. run Capivara spectral, kinematic-aware, and path-signature segmentations;
7. fit a local NIRVANA-style bisymmetric velocity model to the native velocity
   map;
8. save plots, tables, logs, and RDS files for reproducibility.

DAP MAPS are only useful as calibration targets while validating orientation,
velocity scale, and the bisymmetric equation against published figures.

## Where Things Live

- Local MaNGA redshift table:
  `extensions/capivaraKinematics/inst/extdata/manga_drpall_v3_1_1_redshifts.csv.gz`
- MaNGA metadata helpers:
  `extensions/capivaraKinematics/R/manga_metadata.R`
- One-call native kinematics + bar wrapper:
  `extensions/capivaraKinematics/R/run_manga_bar_model.R`
- One-cube RStudio example:
  `extensions/capivaraKinematics/examples/run_one_cube.R`
- Teaching tutorial (segmentation, regional spectra, quick spectral fit,
  kinematics, and the bisymmetric model):
  `extensions/capivaraKinematics/examples/tutorials/capivara_full_workflow.R`

The bundled redshift table is a compact DRPall v3_1_1-derived table with
`plateifu`, `mangaid`, redshift, sky position, and DRP quality metadata for
the MaNGA science rows. This is DRP metadata, not DAP modelling output.

## Simple R Call

Once installed, the whole workflow is one function call:

```r
library(capivara)
library(capivaraKinematics)

result <- run_manga_bar_model(
  cube_path = "/path/to/manga-8078-12703-LOGCUBE.fits",
  redshift = NA_real_,
  segmentation_mode = "all",
  knn_k = 100,
  n_segments = 25,
  n_path_segments = 45,
  show_plots = TRUE
)

print(result)
plot(result)
plot(result, "components")
```

When `redshift = NA_real_`, the resolver uses:

1. a finite manual value, when supplied;
2. the LOGCUBE FITS header, when present;
3. the local bundled DRPall table.

Quick metadata lookups:

```r
manga_redshift("8078-12703")
manga_metadata("10218-12703")
resolve_manga_redshift("/path/to/manga-8932-3701-LOGCUBE.fits")
```

## Model Ownership

The code here is ours: R readers, starlet/better-mask usage, least-squares
piecewise bisymmetric fitting, plotting, and comparison scripts. The velocity
equation and angle convention follow the published NIRVANA/Spekkens-Sellwood
bisymmetric model, so this should be described as a local NIRVANA-style model,
not a full reimplementation of NIRVANA. It does not yet include NIRVANA's
Bayesian sampler, beam smearing, or full dispersion likelihood.

## RStudio Teaching Script

Open `examples/tutorials/capivara_full_workflow.R`, edit the small input block
at the top, and click **Source**. It shows every stage in RStudio and writes
the resulting figures, tables, and RDS objects beside the input cube.

The 8078 paper-calibration runners are retained under the repository's
`research/manga/` directory; they are not part of the supported user API.
