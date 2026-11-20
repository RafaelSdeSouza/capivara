# RFC 0003 — Validation, documentation, visualization and release gates

**Status:** proposed quality/release policy; no current release implied.  
**Companion:** [Architecture](0001-platform-architecture.md) / [Scientific API](0002-scientific-api-specification.md).

## 1. Principle: trustworthy science before visual polish, then exceptional visual communication

Use the best patterns of Astropy/Photutils (scientific concepts and worked examples), scikit-learn (consistent estimator API, user guide and testing), specutils (astronomical spectral quantities), Matplotlib/ggplot2 (publication-quality reproducible figures), pkgdown (R package reference) and PyData Sphinx (Python reference). They are **design references**, not claims that CAPIVARA already supports their features or citation counts.

A result should be reproducible without a GUI; a figure should be an unambiguous view of that result. The same method must be understandable through its equation, scientific intuition, practical example, uncertainty and failure mode.

## 2. Documentation information architecture

**One scientific gateway**, with language-specific API and installation instructions. Keep the current canonical R website `https://rafaelsdesouza.com.br/capivara/` and published route structure functioning while iterating; Python docs may initially be hosted at a sibling versioned path. Avoid a disruptive site generator switch before documenting existing R exports.

Proposed navigation (ordered for user goals):

1. **Overview** — what CAPIVARA solves; one honest hero figure from a real documented observation; core input→representation→segmentation→regional products; current validation status; clear R/Python entrypoints.
2. **Quickstart** — one small, downloadable or simulated fixture with actual expected outputs and both languages; installation with supported dependencies and troubleshooting.
3. **Scientific concepts** — *why* segmentation, *what distance means*, observed-entry missingness, spatial support versus physical association, spectral shapes and invariance, PSF/LSF, WCS and covariance, and result interpretation.
4. **Examples / gallery** — reproducible MaNGA IFU workflow, spectral region map + representative spectra, quality/variance demonstration, exact-versus-sparse objective comparison, kinematic fields, simulation failure modes, export to fitting/DS9 and conditional PSF experiment.
5. **User guide** — end-to-end operations organized by task, not alphabetically: ingest, inspect, select support, choose representation, segment, validate, analyze, plot, export.
6. **Reference** — complete public R and Python API, their status, equations, parameters, return types, accepted units and cross-links.
7. **Benchmarks & verification** — runtime/memory alongside algorithmic comparability; numerical parity manifests and observational validations, with negative/null cases.
8. **Research / experimental** — operator hypotheses, references, prerequisites, counterexamples, and separate validation status; no claims of production readiness.
9. **Cite / About / Contributing** — MNRAS original paper, software release DOI when available, citation instructions, changelog, authorship, funding, license, governance.

Keep reference pages accessible from search and deep links. Avoid hiding methods or examples merely to simplify navigation. Distinguish “recommended for this science question” from “exported, therefore certified.”

### Per-function reference template

```text
FUNCTION NAME — one-sentence scientific question
Intuition: what comparison or estimation does it perform?
When to use / when NOT to use
Mathematical definition or precise algorithm, with equation and variables
Inputs: shapes, coordinate units/frames, valid samples, uncertainty and instruments
Output: objects, region/map semantics, provenance and quality diagnostics
Assumptions and invariances; information removed or unavailable
Performance, memory and numerical behavior (including tie policy)
Reproducible minimal example + one real-data scientific example
Validation status: named fixture/instrument/scope; known failure modes
Related methods / original paper / cite this workflow
```

The user guide must not paste 41 long function signature dumps. References carry exhaustive arguments; articles teach decisions and reasoning.

### Example-first philosophy

Each example has a *question*, measured dataset or synthetic generative truth, declared coordinate/uncertainty model, minimum code, expected output figure, possible failure, interpretation and source/citation. Examples should run on modest hardware without remote credentials; larger benchmarks are separate.

### Current-site migration

First audit existing `_pkgdown.yml`, `vignettes/`, `docs/reference/` and saved URLs. Write a redirect map and backlink check. Roxygen is source for R function references; rendered `docs/` output is generated, not a parallel hand-edited truth. Python Sphinx uses the same vocabulary and can cross-link conceptual pages. Generated documentation must be reproducible in CI.

## 3. Scientific visual design system

CAPIVARA should be visually distinctive *because its diagnostic plots encode correct mathematics*, not because they contain decoration.

**Required figure families:**

- **Four-panel scientific audit:** native white-light/line signal; categorical segmentation with numbered linked regions; representative or all-region spectra; spatial residual/diagnostic map, using independent scientific color scales.
- **Representation comparison:** original measured spectra, explicitly transformed spectra, retained amplitudes, eligibility and missing-channel masks; mathematically identified invariances explained in caption.
- **Instrument-resolution demonstration:** intrinsic synthetic truth, convolved observed cube, observed line profiles and the failed commutation of spatial PSF with variable velocity alignment. Display known truth only for simulations.
- **Hierarchy and scale:** dendrogram or graph contacts, merge criterion, actual region count, disconnection/inversion/tie flags with scientifically legible labels.
- **Kinematic scientific panels:** line flux, line-of-sight velocity, observed dispersion, declared instrumental LSF, model residuals, uncertainty and support, where available. Do not equate a velocity moment with intrinsic thermal/turbulent dispersion.

**Consistency rules:**
- Same region has the same categorical color across maps/spectra/legends; no ranked interpretation of arbitrary labels.
- Always retain astronomical pixel/sky conventions, WCS, axis orientation and spatial scale bars where genuinely known.
- Correct spectral units, air/vacuum/frame, flux density versus integrated flux; uncertainty ribbons only when actually propagated.
- Disclose interpolation/smoothing, clipping, log/asinh scaling and whether shown contours represent measured support or segmentation.
- Make curves and figures large enough for the *primary scientific argument*: no tiny essential plot inside an enormous equation/legend; figure-first layout for examples.
- Vector SVG/PDF, publication-ready raster, adjustable DPI/fonts, meaningful colorbars, color-vision-safe palettes, accessible alt text and non-color discriminants for key comparisons.
- Animation *only when it demonstrates a checkable property*: shift covariance, invariance, PSF mixing or information loss; include controls, captions and still-frame equivalent; not autoplay-only evidence.
- Plot helpers must accept existing Matplotlib axes or ggplot-like composition as appropriate, return plotting objects and never change data or silently rerun an estimator.
- The scientific gallery must disclose simulation versus real observation, input release and image/source credit.

**Visual QA** should include sample outputs inspected at phone and desktop widths, snapshot tests for legend/axis/region consistency, and checks that crucial math and figures are readable without speaker notes or tooltips.

## 4. Scientific validation tiers

| Tier | Meaning | Required evidence |
|---|---|---|
| V0 — experimental | Mathematical conjecture, algorithm spike or uncalibrated workflow | Research namespace, explicit failure cases, no validated astrophysical output |
| V1 — numerical | Internal mathematical and software contracts verified | Analytic fixtures, property-based tests, unit/frame checks, tolerance report, CI |
| V2 — R/Python parity (for applicable paths) | Same declared algorithm on same measured samples | Golden fixtures, co-membership, region products, objective values, full provenance, tie/precision policy |
| V3 — instrument-scoped | Tested on a named instrument and data product | FITS round-trips, WCS/LSF/noise, reproducible real cubes, known quality/coverage limits |
| V4 — scientific utility | Demonstrated advantage or calibration for a stated scientific question | Strong baselines, held-out/simulation truth, negative controls and uncertainty/falsification tests |

These are **orthogonal assessments**, not a single misleading “validated” Boolean. Example: a segmentation may pass numerical contract V1 but have no claim of identifying physical components at V4.

## 5. Minimum tests by domain

### Observations and geometry

- Strict shape/axis contracts; R column-major versus Python row-major flattening and round-trip; one-based/zero-based channel indices; native/sky/physical orientation and East-left convention.
- Wavelength air/vacuum, rest/observed (and explicit systemic redshift), frequency/wavelength/f_λ/f_ν; conversion metadata, native-channel selection and no silent resampling.
- FITS masks/IVAR/variance/DQ sentinels with correct data-extension association; input values retained unmodified.
- WCS and geometry for rectangular grids, holes/disconnected components and mixed reference axes.
- LSF pre/post pixel sampling, wavelength-dependent resolution and missing LSF behavior.

### Representations, support and segmentation

- Positive-homogeneous shape invariance and separate amplitude; S/N threshold and missing amplitude window; high-noise negatives and invalid signal.
- Quality support, detection support, host association and representation eligibility distinguishable by exact synthetic labels.
- Historical L1 + Ward.D2, Euclidean sparse Ward and observed-entry hierarchy separate; no deceptive equivalence labels. Fixed ties/inversions and disconnected graphs documented.
- Partition tests compare co-membership, not raw arbitrary region colors/IDs; explicit support and exact tested counting.
- Memory growth baseline and optional accelerated backend agreement, separated from scientific accuracy.

### Region summaries and fitting handoffs

- Sum/mean/median/weighted estimands clearly different; flux conservation only under the appropriate reconstruction.
- Masked contribution with measured zero versus absent sample; partial/no variance must not produce artificially precise outputs.
- Correlated uncertainty test using `Var(sum_i X_i) = 1ᵀΣ1`; diagonal-only algorithms explicitly label approximation.
- FITS/export round-trip, line identity, channel order, pPXF response and wavelength provenance.
- Physical negative controls: velocity gradient without change of intrinsic line shape; changed ratios at fixed velocity; overlapping multicomponent lines; changing PSF/LSF.

### Visualization

- Categorical IDs, WCS orientation, axis labels, units, flux-conserving label consistency, residual semantics and accessibly annotated figures; synthetic truth clearly declared.

## 6. R ↔ Python golden reference manifest (proposed)

The stored fixture must record at least:

```yaml
fixture_id: legacy_manga_small_v1
source:
  repository: RafaelSdeSouza/capivara
  r_commit_sha: "<frozen reference sha>"
  input_release: "<survey/product release or synthetic seed>"
axes:
  order: [row, column, spectral_channel]
  wavelength_unit: Angstrom
  wavelength_frame: observed
  wavelength_medium: vacuum
  channel_index_base: 1
observation:
  flux_unit: "<declared or explicitly unspecified>"
  validity_dataset: "sample_validity"
  variance_kind: diagonal_variance_or_none
method:
  method_id: legacy_l1_wardd2
  representation_id: historic_median_center
  requested_regions: 8
  missingness_policy: historical_zero_fill_after_transform
  tie_policy: "<algorithm-specific>"
verification:
  compare: [co_membership, selected_native_channels, region_sum_flux]
  numerical_tolerances: "<justified per quantity>"
status: proposed_fixture_not_yet_frozen
```

Illustrative schema, not a completed fixture. Every fixture includes its actual datasets, SHA/checksums, command and environment lock information in a future implementation phase.

## 7. Release and CI gates

**Gate A — baseline:** frozen `main` reference SHA, package check, current public export inventory, historical figures, issue audit and known scientific defects. No implementation PR should rewrite history before this.
**Gate B — R:** changes preserve all intended exports or provide versioned deprecation; testthat+R CMD check, source/install equivalence, scientific product round trips and instrument-scoped QC.
**Gate C — Python:** wheel installs cleanly, current minimal-versions test matrix, PyPI name/license/citation confirmed, all supported paths have scientific parity or prominent “Python-specific / not yet paired” labels. No broad parity claim.
**Gate D — docs:** pkgdown/Sphinx sites build, links/redirects pass, examples execute, screenshots reviewed, meaningful scientific figure and cite section ready.
**Gate E — research:** novelty check (2–3 strongest methods), nontrivial math result and independent negative/positive tests; no release of new statistical operator without documented validation.

**Versioned releases:** published software DOI and `CITATION.cff` if the release is archived, original MNRAS `https://doi.org/10.1093/mnras/staf688` always cited with reproducible method ID, separate citation for an independent future methods paper if warranted. Do not claim “peer reviewed software” merely because original method paper exists.

## 8. Proposed reviewable implementation order — no timelines or commitments

1. Freeze scientific baselines and correct/confirm any still-open correctness issues **before** redesigning user-visible workflows.
2. Consolidate scientific metadata and shared golden-fixture schema without breaking old R inputs.
3. Port the smallest native Python feature path with parity, errors and plots.
4. Regenerate R/Python docs and publish an exemplary real-data gallery.
5. Add or improve advanced kinematic analysis modules only when validated.
6. Evaluate the spectro-spatial research proposal independently of the production release.

## 9. Final launch criteria

The scientific platform qualifies for a stable public release when a new user can:
- install either supported language from documented steps,
- run and understand a named tested workflow, including input units/frames and uncertainties,
- reproduce one real-data figure plus an independent failure case,
- discover **every** current public function and why it exists,
- trace a result to algorithmic objective, package version, upstream data and validation evidence,
- see which mathematical operators are proven, which are empirical and which remain research ideas.

Documentation-only changes are allowed now. **Do not merge the earlier Python prototype as a finished port and do not deploy a redesigned website or stable R API until the scientific specifications are approved.**