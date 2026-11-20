# RFC 0001 — CAPIVARA scientific-platform architecture

**Status:** proposed architecture for scientific review; **not** authorization to refactor the stable R code or release the Python prototype.  
**Date:** 2026-11-20  
**Project:** CAPIVARA — Rafael S. de Souza, Creator and Scientific Lead  
**Baseline:** current `main` R package reports `0.4.5.9000`; the existing Python skeleton lives only on `feat/platform-foundations-r-python-math`.  
**Related:** [RFC 0002](0002-scientific-api-specification.md) (scientific API), [RFC 0003](0003-validation-documentation-release.md) (tests, website and release gates).

## Decision

Build **one scientific specification with two native implementations** and an independently evaluated research layer. Preserve the published R methods, numerical history, all useful existing capabilities, legacy names and the live documentation URLs while modernizing the internal architecture. Python is a native implementation, not an R wrapper and not a claim of numerical parity until verified.

**Scientific identity:** geometry and statistics of spatially resolved spectra. The mature scientific use case is segmentation and flux-aware regional analysis of IFU cubes. Kinematics and structured representations are complementary. Observation-aware spectro-spatial statistical operators are research candidates, **not** stable product claims.

**Out of scope for the core:** replacing Astropy/specutils, pPXF, full radiative-transfer solvers or survey reduction pipelines; claiming that a spectral partition identifies physically distinct stellar populations; shipping speculative geometric methods without independent validation.

## Design rules (non-negotiable)

1. **An estimand precedes an algorithm.** Every public scientific method declares its input quantity, intended question, objective/functionals, assumptions, invariances, what it discards, output meaning, and known failure modes.
2. **Measurement and interpretation are different objects.** A wavelength channel, a valid voxel, a detected host, a spatial support mask, a spectral region, a fitted kinematic model and a physical structure are not aliases.
3. **No silent physical operation.** Redshift-frame conversion, air/vacuum conversion, resampling, unit conversion, PSF/LSF homogenization, noise inflation, emission-line selection, spatial mask intersection and model selection must be explicit and logged.
4. **Never turn missing data into a measured zero.** Any numerical fill used inside a historical method is labelled as an algorithmic convention; it is never a new observation.
5. **Historical reproducibility.** `segment()` legacy L1 + Ward.D2, sparse Euclidean centroid Ward on a graph, and the observed-entry spatial hierarchy have distinct mathematical objectives. Each retains a separate method ID and reference fixtures.
6. **A claim has a verification status.** `verified_numeric`, `supported_observational`, `experimental`, `legacy_reproducibility` and `unvalidated` are not interchangeable. Status is specific to a method, instrument, assumptions and test dataset.
7. **Reviewable evolution.** No breaking public-API changes outside an announced migration cycle; changes require scientific regression tests, deprecation guidance and an explicit changelog.
8. **Separation of methods and presentation.** Visualization reads validated result objects; it never re-estimates fluxes, reorders region identities silently or alters scientific data.
9. **No forced single engine.** NumPy/SciPy and R/fastcluster remain independently implementable; optional acceleration must reproduce the declared algorithm to a documented tolerance.
10. **Scientific users first.** A first-time astronomer must be able to understand what a function answers and what it does not establish before supplying 30 options.

## Layered model

```text
   Native instrument products (MaNGA / MUSE / generic IFU)
                       |
                ADAPTERS / INGESTION
                       v
 Observation = flux + coordinates + measured validity + uncertainty
               + PSF/LSF/WCS provenance (if known)
                       |
            QUALITY / DETECTION / SUPPORT
                       |
              EXPLICIT REPRESENTATION
                /                 \
       SEGMENTATION         KINEMATICS (separate)
                \                 /
               REGIONAL PRODUCTS / DIAGNOSTICS
                       |
             SCIENTIFIC FIGURES & EXPORT
                       |
          OPTIONAL EXPERIMENTAL OPERATORS
          (independent, gated, not assumed)
```

### Scientific domains and ownership

| Domain | Purpose | Public/stability treatment | Dependencies |
|---|---|---|---|
| `observation` | Ingest, validate and record measurement coordinates, units, WCS, masks, variance, PSF/LSF | Foundational stable contract, initially adapters may be limited | R FITSio/optional astronomy readers; Python astropy/specutils/spectral-cube adapters |
| `support` | Measurement validity, foreground detection, host selection, representation eligibility and analysis support | Existing masks remain; decomposition of roles mandatory | No required learned source detector |
| `representation` | Construct quantity being compared: full flux, shape, line profile; declare invariance/removed information | Existing representation declarations promoted as the conceptual centre; preserve older names | Instrument/line coordinate specifications |
| `segmentation` | Hierarchy, objective, graph eligibility, cut, region IDs and diagnostics | Legacy, sparse and observed-entry paths remain distinct; do not rebrand one as another | Deterministic reference implementations |
| `regions` | Sum, mean, weighted spectra, uncertainty, flux-conserving products, adjacency | Stable core after scientific regression gates | Observation + SegmentationResult |
| `kinematics` | Line moments/kinematics and separately validated dynamical comparisons | Existing exports retained; validation scope clearly labelled | Spectrum/LSF/covariances |
| `visualization` | Publication-quality maps, spectra, uncertainty and QC displays | Dedicated design system; composable plots, SVG/PDF export | Read-only results |
| `experimental` | Transport-aware statistics, observation-operator/covariance research, identification diagnostics | Python `capivara.experimental`; R research/internal until independent approval | Explicit mathematical specifications and independent baselines |

### Companion project policy

- **SAGUI:** remains a separately identifiable photometric-segmentation research project; share semantics and validation fixtures where practical, do not absorb it into an untested IFU API.
- **SpectroPath:** independent mathematics dependency for ordered profile paths/signatures, with a public mathematical explanation and provenance.
- **capivaraPPXF:** optional downstream fitting integration. No broad fitting claims until LSF, noise, template and uncertainty contracts pass.
- **ModelAtlas:** optional visualization interface, not a hard dependency or substitute for CAPIVARA's own plotting baseline.

## R and Python boundary

**R:** Existing namespace and user workflows remain the current source of truth; keep public `segment()` and `segment_large()` signatures compatible. Add compositional internals and compatibility shims behind existing functions only after fixtures pass. Retain ggplot2 + pkgdown. Do not require Python, reticulate or a GPU for core segmentation.

**Python:** Native `capivara` distribution, with typed `Observation`, `Representation`, `SegmentationResult` and `RegionalSpectra` concepts, standard NumPy array protocol and explicit Astropy/`specutils.Spectrum` adapters. Dependency core is NumPy/SciPy with optional astronomical, visualization, GPU and scientific-fitting extras. The skeleton on the earlier experimental branch is a feasibility spike, not a release or parity certificate.

**Cross-language:** Share a language-neutral contract document and serialized golden fixtures (array shapes, wavelength/medium/frame, validity, labels, objective, uncertainty, provenance). Check equivalence of region co-membership rather than raw label integers; use exact comparisons where meaningful and explicit numerical tolerances otherwise. Differences in tie order, algorithms, missing-data decisions and floating-point libraries must be surfaced, not hidden.

**Repo layout proposal (not implemented):**

```text
R/  src/  tests/testthat/  vignettes/    # existing R package
python/src/capivara/ python/tests/      # independent Python package
contracts/                              # versioned scientific semantics + fixtures schema
validation/                             # reproducible scientific gates and manifests
examples/                               # small public scientific datasets/recipes
rfcs/                                   # design decisions and review
docs/                                   # existing published R site (do not overwrite blindly)
```

Use one repo during development to make parity review tractable; revisit a split only if release cadence or governance requires it. No shared C++ backend is mandated up front: test scientific equivalence first, optimize after profiling.

## Scientific metadata schema

Every scientific result must be able to describe: `method_id`; `method_version`; `software_version` and source SHA; `observation_id` and origin; `coordinate_frame`, `spectral_medium`, `spectral_axis_unit`, `flux_unit`; spectral channel-selection provenance; voxel validity and spatial support lineage; `representation_id` with retained/removed quantities and claimed invariances; uncertainty model; instrument WCS/PSF/LSF if supplied (otherwise explicit unknown); objective and backend parameters; numerical tolerance and tie policy; actual region count; warnings/QC; input and output hashes when available; `validation_status` linked to evidence.

**Unknown is not zero.** Missing PSF, unavailable covariance and unknown flux units must remain explicit unknowns. A warning cannot make an unsupported physical claim valid.

## Versioning and compatibility

- Maintain legacy **public behavior** on the established R paths. Where correctness fixes change a scientifically meaningful result, release with versioned migration notes and frozen before/after fixtures.
- Semantic versioning for eventual tagged R and Python releases separately. Do not declare a “2.0” release merely because design documents exist; `0.4.5.9000` is the repository development version observed when writing this RFC.
- During transition, experimental functions use visibly distinct namespaces and documentation status; no unstable method is advertised as equivalent to published CAPIVARA.
- Maintain previously published website routes and canonical citation URLs; use redirects if the navigation changes.

## Key architecture risks

1. **Scientific drift:** `legacy L1+Ward.D2`, graph Ward and observed-entry SSE have different semantics; enforcing one algorithm would break results.
2. **Implicit observational conventions:** IFU orientation, air/vacuum, observed/rest, flux units and LSF can invalidate analyses while producing plausible plots.
3. **Two implementations, no parity:** shared prose alone does not guarantee numerical equivalence; the reference manifest is a release blocker.
4. **Too much scope:** the core should not grow to replace pPXF or full survey pipelines.
5. **Premature novelty:** commutators, equivariance and information monotonicity are known mathematics; any new spectro-spatial claim needs a prior-art and falsification gate.

## Decisions explicitly deferred to implementation RFCs

- Python scientific minimum versions and acceleration stack after compatibility CI survey (follow Scientific Python SPEC 0 rather than freezing arbitrary minimums).
- A final stable `Observation` class signature, conversion policy, ownership and serialization format after adapter experiments.
- Whether to promote a new observed-entry API to the default (never silently alter historical `segment()`).
- The form of spectro-spatial kernels/operators; novelty, approximation guarantees and validation must be demonstrated before naming a public algorithm.
- Splitting R/Python release repositories, if maintainability justifies it.

## Required authorization gate

This RFC and its companion API/validation RFCs can be reviewed and revised as documents. **No stable R refactor, Python parity claim, public website rebuild, automated package release or experimental-method promotion follows solely from writing them.** Implementation begins with separate scoped PRs after design approval.