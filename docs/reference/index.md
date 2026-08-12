# Package index

## Core spectral segmentation

Exact and sparse-graph segmentation, component selection, and memory
planning.

- [`segment()`](https://rafaelsdesouza.github.io/capivara/reference/segment.md)
  : Cluster a 2D Representation of a Data Cube
- [`segment_large()`](https://rafaelsdesouza.github.io/capivara/reference/segment_large.md)
  : Large-cube Capivara segmentation
- [`choose_ncomp_by_snr()`](https://rafaelsdesouza.github.io/capivara/reference/choose_ncomp_by_snr.md)
  : Choose the Number of Clusters from a Target SNR Cut
- [`estimate_segment_memory()`](https://rafaelsdesouza.github.io/capivara/reference/estimate_segment_memory.md)
  : Estimate segmentation memory requirements

## Support construction

Foreground support builders used before assigning region labels.

- [`build_starlet_mask()`](https://rafaelsdesouza.github.io/capivara/reference/build_starlet_mask.md)
  : Build a Sagui-style starlet mask from a spectral cube
- [`build_adaptive_support()`](https://rafaelsdesouza.github.io/capivara/reference/build_adaptive_support.md)
  : Build an adaptive multi-band/spatial support mask

## Flux-preserving spectra and reconstruction

Regional spectra and reconstructed cubes for downstream analysis.

- [`summarize_cluster_spectra()`](https://rafaelsdesouza.github.io/capivara/reference/summarize_cluster_spectra.md)
  : Summarize Cluster Spectra with Optional Variance Propagation
- [`reconstruct_cluster_cube()`](https://rafaelsdesouza.github.io/capivara/reference/reconstruct_cluster_cube.md)
  : Reconstruct a Model Cube from Cluster Representative Spectra
- [`reconstruct_flux_preserving_cube()`](https://rafaelsdesouza.github.io/capivara/reference/reconstruct_flux_preserving_cube.md)
  : Reconstruct a Flux-preserving Model Cube from Cluster Spectra

## Segmentation plotting

Discrete region maps and regional spectral displays.

- [`plot_cluster()`](https://rafaelsdesouza.github.io/capivara/reference/plot_cluster.md)
  : Plot a Cluster Map with Discrete Cluster Colors
- [`plot_cluster_spectra()`](https://rafaelsdesouza.github.io/capivara/reference/plot_cluster_spectra.md)
  : Visualize Scaled Spectra and Median Profiles for Clustered IFU Data

## Kinematics

Native line maps, kinematic segmentation, model discovery, and reusable
panels.

- [`segment_kinematics()`](https://rafaelsdesouza.github.io/capivara/reference/segment_kinematics.md)
  : Segment an IFU cube using emission-line kinematics
- [`run_kinematic_analysis()`](https://rafaelsdesouza.github.io/capivara/reference/run_kinematic_analysis.md)
  : Run a Capivara kinematic analysis
- [`kinematic_models()`](https://rafaelsdesouza.github.io/capivara/reference/kinematic_models.md)
  : List the installed kinematic model modules
- [`kinematic_panels()`](https://rafaelsdesouza.github.io/capivara/reference/kinematic_panels.md)
  : Extract individual kinematic plot panels

## Bisymmetric modelling

Explicit modelling for galaxies with independent evidence of a bar.

- [`run_manga_bar_model()`](https://rafaelsdesouza.github.io/capivara/reference/run_manga_bar_model.md)
  : Run the explicit bisymmetric-bar module for a MaNGA cube

## Package

- [`capivara`](https://rafaelsdesouza.github.io/capivara/reference/capivara-package.md)
  [`capivara-package`](https://rafaelsdesouza.github.io/capivara/reference/capivara-package.md)
  : capivara
