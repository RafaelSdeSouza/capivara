# Repository Layout

The repository contains one installable R package and a separate research
area. Keeping those roles distinct prevents experimental outputs and
galaxy-specific scripts from leaking into the user-facing API.

| Location | Role | User-facing? |
| --- | --- | --- |
| `R/`, `man/`, `inst/` | **capivara**: segmentation, regional spectra, kinematic/path-aware maps, and disc/bisymmetric models | Yes |
| `inst/tutorials/` | Minimal RStudio scripts that demonstrate the supported workflow | Yes |
| `research/` | Paper calibration, benchmarks, and one-off galaxy experiments | No |
| `images/` | Documentation images referenced by the README or pkgdown site | Yes |
| `outputs/` | Generated data products; ignored by Git and package builds | No |

The former `capivaraKinematics` prototype is retained under `research/` only
while its validation history remains useful. It is not a second installable
user-facing package. The core kinematics module writes all run products beside
the selected input cube (or to a user-supplied output directory).
