# RFC 0006 — Agreed visual direction and first real-data scientific atlas

**Status:** accepted **design direction**, proposed **real-data implementation protocol**, 20 November 2026. This is a documentation-only record; it does not mean a deployed website or validated interactive atlas exists.

Related: [Scientific architecture](0001-platform-architecture.md), [API and mathematical contracts](0002-scientific-api-specification.md), [Visual/editorial standards](0004-scientific-editorial-and-geometric-statistics.md), [External benchmark](0005-external-benchmark-and-editorial-copy.md).

## The direction agreed in review

CAPIVARA will use a **predominantly white scientific editorial site** with vivid, disciplined color from the astronomical data. The central interactive object is the **actual geometrical union of the spaxels assigned to a region**, not a raster image made to look like a smoothed polygon. Each region is selectable and is linked to its **measured** regional spectrum and its provenance.

It should combine:
- documentation quality and onboarding comparable to Astropy and Photutils;
- a vivid, well-captioned real-data gallery as in SunPy, without hiding scientific units or uncertainty;
- the mathematical precision and operator demonstrations of the Astro Informatics/S2FFT/S2WAV work;
- interactive, linked region/spectrum behavior inspired by CARTA and Cubeviz;
- separate and correctly labelled 3D spectral-cube, physical-volume, and *feature-space* views.

**Avoid**: default dark documentation backgrounds, washed-out color, pseudo-depth brightness mountains, fabricated sky geometry, interpolated/smoothed segment boundaries presented as measured, generic “unlock the universe” prose, decorative formulas, arbitrary 3D point clouds presented as physical manifolds, and a homepage showing a simulated galaxy as though observed.

The overall site's written voice is that of an astronomical analysis manual and a good methods paper: short declarative sentences, defined physical quantities, clear reasons for each method, limitations next to results, and an explicit distinction between numerical tests and scientific validation.

## First real-data candidate: MaNGA 7443-12703 (Mrk 848)

**Why this target?** It is a recognisable interacting/merging galaxy in SDSS's own MaNGA/Marvin tutorials. It is also one of the five native MaNGA objects used in CAPIVARA's V4 spatial and semantic-representation validation, avoiding an arbitrary choice made only for appearance. This is a **candidate for the showcase**, not an assertion that the completed linked products exist in this branch.

**Official source:** [SDSS DR17 MaNGA data access](https://www.sdss4.org/dr17/manga/manga-data/data-access/), [MaNGA DR17 data model](https://www.sdss4.org/dr17/manga/manga-data/data-model/), [SDSS Marvin tutorial for 7443-12703](https://www.sdss4.org/dr17/manga/manga-tutorials/marvin-tutorial/marvin-tools).

Official DR17 DRP spectral file location (for reproducibility, **not** an independently downloaded file):

```text
https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/7443/stack/manga-7443-12703-LOGCUBE.fits.gz
```

**Existing source-backed capability:** `regions_sf()` uses exact labelled-cell unions (possibly MULTIPOLYGON and/or holes), `domain_sf()` represents support, and `region_adjacency()` records neighbours. `plot_cluster()` already renders these exact polygons and transforms WCS for display without altering canonical geometry. `summarize_cluster_spectra()` returns region sums and counts; possible variance/covariance limits must remain visible. These are existing R methods, not new planned algorithms.

**Frozen evidence:** `SPATIAL_FOUNDATION_V4_FREEZE.md` reports five-object native validation and exact polygon round trips. Its underlying accepted per-object arrays and results are documented under external validation/worktree paths, and are **not guaranteed to be committed to this public repository**. Before generating the public atlas, either retrieve those accepted products with their provenance/hashes or re-run the frozen method on the official DR17 cube, recording differences. Do not derive a fake segmentation from the SDSS photograph.

**Representation choice:** First viewer should use an explicitly stated, instrument-scoped validated method (e.g., `manga_sandra_spectral_shape()` with its exact rest-frame anchor windows and measured-amplitude S/N ≥30 eligibility, or the frozen full-flux method), then a simple control to show a second existing validated representation. Do not mix incompatible region IDs or claim that two partitions are the same.

## The scientific viewer: minimum genuine prototype

**Layout on a phone:**
1. White CAPIVARA masthead, compact page title: **“Where do the spectra change?”** Short deck describing the actual data/estimand; source and method attached to the figure.
2. Large real MaNGA broadband/white-light image with an optional colored segmentation overlay. Segment boundaries are **stroke-rendered exact polygons** from `regions_sf()`, with selection-friendly hit areas that do not alter the displayed boundaries. No raster approximation at interaction time.
3. Tap a region to outline it prominently and scroll/update its **summed observed regional spectrum** below; optionally switch to mean/variance-aware display if supported. The selected region's color is consistent in map, spectrum and metadata.
4. A compact metadata row gives region ID, contributing spaxel count, spectral coverage, declared frame/units, method and support provenance; show unavailable uncertainty as “not available”, not 0.
5. An optional “What am I looking at?” paragraph explains why a spectral region is not necessarily a physical component.

**Layout on desktop:** The same sky map and regional spectrum sit side by side and both occupy substantial space. Selection is bidirectional where justified: click polygon to select spectrum; select region/spectrum entry to highlight polygon. Keyboard and touch must both work.

**Further views, but not required for the first prototype:**
- **Cube slices**: correctly labelled `(x,y,λ)`; wavelength is *not* line-of-sight physical depth; show rest/observed frame and native-sampling details.
- **Feature-space view**: PCA as a diagnostic with explained preprocessing, explained variance and a measured mapping from spaxel to points; do not call a 3D PCA point cloud a physical manifold. UMAP/t-SNE may be supplementary but not validation evidence.
- **Observation-operator lab**: show transformation actions and resolution mixing with explicit analytic/simulation tests; separate from the stable case study.
- **True volume rendering/ray tracing**: only when data describe a physically justified 3D volume or reconstruction, with rendering settings and transfer function stated.

## Data pipeline and rendering contract

Read the DR17 **DRP** LOGCUBE using a validated MaNGA adapter: FLUX, WAVE, IVAR, MASK, per-spaxel LSF, available WCS and spatial-correlation products; no assumption that the DAP model cube and measured DRP cube are interchangeable.

For `segment_result` with declared method + support:
- Generate **canonical** region geometries via `regions_sf()` in a stated coordinate mode, preserving connected components, holes and true cell footprint edges.
- Export a GeoJSON-like feature collection or SVG-ready polygons with stable `region_id`, topology type and coordinate metadata. **Do not use GeoJSON latitude/longitude conventions for local/native or tangent-plane coordinates without an explicit CRS contract**. For sky plotting, East-left is a *view transform once*, not a coordinate rewrite.
- Associate `region_id` (not ordered color indices or string positions) with matching measured sum, mean, channel count, valid/invalid bits, declared units and uncertainty flags.
- Record `method_id`, `Ncomp`/actual_k, input SHA, source revision, spatial and spectral coordinates, representation, and assumed/included covariance. Preserve holes/IDs through serialization.
- Generate a static vector fallback image with the same science/caption, plus an interactive SVG/canvas map. For 8–25 regions SVG is appropriate; profile before introducing larger WebGL dependencies.
- Fetching data must be reproducible, cached and separate from live site requests. Deliver a small precomputed downsampled view for mobile, with full-resolution downloadable science products.

**No numerical inference in browser event handlers.** Browser interactions only select or display *already computed* science arrays and derived maps.

## First viewer's acceptance tests

- Exactly the same pixel memberships before/after polygon export (including disconnected components and holes); area and adjacency agree with `regions_sf()`.
- Region selection maps to the **correct** summed observed spectrum and region ID, including after changing display orientation and coordinate mode; no wrong indexing due to R/Python order.
- Every region has correct channel counts; missing flux and partially missing variance are never reported as zero certainty.
- Spectral-axis units, frame and medium agree with the input; plotted brightness scaling and data range are disclosed.
- Categorical region colors are stable across map and spectrum and distinguishable; a palette alone is not used for scientific rank.
- Keyboard and tap targets work at roughly 390px mobile width. No horizontal scroll or illegible spectroscopic plot axes.
- Static PDF/SVG scientific figure and interactive view agree on polygon geometry and data labels.
- Product clearly labelled **MaNGA DR17** + **CAPIVARA algorithm/version** + **real observation**, with a direct link to the MNRAS paper and original SDSS source.
- Production website and stable R API are unchanged until a dedicated implementation/validation review.

## Important known limitations

No genuine interactive MaNGA atlas has yet been produced in this RFC branch. The existing white-background prototypes use synthetic data, and the source/validated per-object outputs needed for a real linked viewer must still be acquired or regenerated. This design decision is **not** evidence that a new spectro-spatial operator exists.

The experimental mathematical programme remains separate: representation group actions (with named target), observational factorization conditions, covariance/statistics under PSF/LSF and identifiability of scientific functionals. Classical results are credited; promotion requires a distinct theorem/algorithm with baselines and astronomical validation.

## Next execution decision

**First build and review one real-case viewer for 7443-12703 using a frozen CAPIVARA method, exact polygons and matching measured regional spectra.** Do it on an isolated implementation branch and/or local preview. Only after it passes geometry, WCS, spectral and mobile tests should the website architecture be modified or the viewer become a reusable CAPIVARA visualization capability.
