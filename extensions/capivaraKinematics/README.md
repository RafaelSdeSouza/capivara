# capivaraKinematics prototype

Experimental companion workflow for Capivara-native kinematic modelling of
barred MaNGA galaxies.

This code intentionally lives outside the core Capivara package. The goal is to
test a segmentation-aware method on real Capivara outputs, then split the stable
parts into a companion/plugin package in the style of `capivaraPPXF`.

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
- RStudio visual example:
  `extensions/capivaraKinematics/examples/visual_manga_bisymmetric.R`

The bundled redshift table is a compact DRPall v3_1_1-derived table with
`plateifu`, `mangaid`, redshift, sky position, and DRP quality metadata for
the MaNGA science rows. This is DRP metadata, not DAP modelling output.

## Simple R Call

For development from the Capivara repo:

```r
repo_root <- "/Users/rd23aag/Documents/GitHub/capivara"

source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "capivara_kinematics_utils.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "manga_metadata.R"))
source(file.path(repo_root, "extensions", "capivaraKinematics", "R", "run_manga_bar_model.R"))

result <- run_manga_bar_model(
  cube_path = "/path/to/manga-8078-12703-LOGCUBE.fits",
  redshift = NA_real_,
  segmentation_mode = "kinematic",
  repo_root = repo_root,
  disc_pa_deg = 14.7,
  disc_inc_deg = 36.4,
  bar_phi_deg = 40.9,
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

## 8078 Native Reproduction

The current native example can be reproduced with one clean script:

```sh
Rscript extensions/capivaraKinematics/examples/native_manga_bisymmetric_8078.R
```

The script has one edit block at the top for `cube_path`, `redshift`,
`output_dir`, geometry, and segmentation settings. It then creates all native
Capivara maps, segmentations, the bisymmetric model, and the figures.

For a new cube, copy and edit the generic template:

```sh
Rscript extensions/capivaraKinematics/examples/native_manga_bisymmetric.R
```

This template is designed to be run by changing only the input block at the top.

The lower-level runner is still available for development:

```sh
Rscript extensions/capivaraKinematics/scripts/run_8078_native_starletfull.R
```

By default this writes to `/private/tmp/capivara_native_8078_starletfull`.
The recipe uses:

- `manga-8078-12703-LOGCUBE.fits`
- `z = 0.0281`
- Halpha
- `knn_k = 100`
- `Ncomp = 25`
- path-signature groups `45`
- starlet scales `1:5`
- coarse starlet plane included
- Halpha peak search `350 km/s`
- centroid window `220 km/s`
- nearest-filled line-map holes are tracked with `halpha_imputed`

No VorBin, Voronoi binning, or S/N-driven spatial binning is used.
