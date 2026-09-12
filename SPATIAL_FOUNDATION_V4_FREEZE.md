# CAPIVARA V4 spatial-foundation integration freeze

This release was developed in the isolated `release/spatial-foundation-v4`
worktree from reviewed production commit
`96159fc3d516a6466e9ee76f778a1a15e603c112`. It promotes only contracts already
accepted in the iFUN five-galaxy validation. No bar model or new segmentation
method enters the package.

Explicit representations now derive their default analysis support from their
declared measurement validity and eligibility. `spectral_shape` therefore uses
the frozen MaNGA/Sandra S/N >= 30 eligibility contract, while `full_flux` uses
its declared validity domain. A supplied support object remains an explicit
additional restriction. Starlet is still selectable, but no explicit semantic
analysis intersects it silently. Results record `support_source`,
`validity_contract`, `eligibility_contract`, and `representation`. Calls with
`representation = NULL` retain the reviewed historical algorithms and record
that legacy provenance separately.

Every partition exposes exact labelled-cell unions through `regions_sf()`,
`domain_sf()`, and `region_adjacency()`. The canonical geometry contains
positive-East and positive-North coordinates transformed from native cell
corners through the declared celestial WCS. With redshift, the spatial contract
records the cosmology, angular scale, area-equivalent pixel scale, and adopted
centre and returns areas and centroids in kpc. Disconnected components and
holes remain polygon topology; no smoothing, interpolation, closing, or
simplification is applied. Adjacency uses positive shared-boundary length and
is checked against native four-neighbour labels. Its graph is retained with the
topology table.

`plot_cluster()` now renders exact `sf` regions. Sky and physical panels use a
temporary display copy with one East-left reflection; canonical geometry and
native matrices are unchanged. Explicit semantic results use the approved Van
Gogh 2.0 identity palette by default. Historical results retain the legacy
`starry_night` default, and the explicit legacy palette remains available.
Continuous, signed, ranked, and nominal data use separate public scales.

Validation used the complete package test suite, an asymmetric 5 x 7 fixture
with 15 labels, direct comparison with the reviewed V4 installation, a clean
source build and installation, and `R CMD check --no-manual`. The reviewed and
new installations return identical legacy outputs and identical V4 hierarchy
children, merge costs, tie decisions, and K=8/K=16 cuts in controlled replay.
The clean installed package also replayed the five native MEGACUBEs. The
accepted `spectral_shape` counts are 2082, 748, 2659, 707, and 1843; the
`full_flux` counts are 2695, 761, 2761, 788, and 1918. All ten K=8 and K=16
label matrices match the accepted files exactly. All ten polygon round trips,
native adjacency comparisons, and single-reflection display checks pass. The
external replay table is
`results/production_spatial_foundation_release/five_object_installed_replay.csv`
in the iFUN validation tree.

```ini
REPRESENTATION_DOMAIN_DEFAULT = READY
STARLET_LEGACY_PATH = PRESERVED
SF_REGION_GEOMETRY = READY
SPATIAL_COORDINATE_CONTRACT = READY
VECTOR_RENDERER = READY
VAN_GOGH_2_DEFAULTS = READY
SOURCE_INSTALL_EQUIVALENCE = PASS
READY_FOR_EXTERNAL_BAR_GEOMETRY = YES
```
