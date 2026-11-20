# RFC 0002 — Scientific API and mathematical contracts

**Status:** architecture proposal; NOT an implementation or parity certificate.  
**Companion:** [RFC 0001](0001-platform-architecture.md) / [RFC 0003](0003-validation-documentation-release.md).  
**Scope:** all **41 exports** in the current R `NAMESPACE` (20 November 2026, R development version `0.4.5.9000`), the future minimal Python core, and a strictly separate experimental namespace.

## 1. Definition of a scientific API

Every user-facing function must answer:

1. **Question and intuition:** What scientific distinction does the result describe, and why is the operation meaningful?
2. **Input quantity:** What is actually measured; coordinate system, units, spectral frame, instrument response, validity and uncertainty?
3. **Mathematical definition:** Formula, objective, transformation, or precisely specified numerical algorithm and possible nonuniqueness.
4. **Invariances and losses:** What can change without affecting the result; which information is intentionally removed or potentially not recoverable?
5. **Output and limits:** Interpretability, validity domain, failure modes, QC, provenance and whether external physical validation exists.

A mathematically correct algorithm can be scientifically inappropriate for an incompatible estimand. A region label is not a stellar-population class, an emission-line moment is not a decomposed multi-component velocity field, a flux-preserving template is not a high-resolution reconstruction, and a mathematical symmetry is not an observational identifiability guarantee.

## 2. Shared conceptual types (design; not mandated class names)

### `Observation` / `IFUCube`

- Primary signal `F[r,c,k]` with canonical **row, column, spectral-channel** order. R and Python memory order differs: canonical pixel IDs and an explicit flattening/unflattening map are required.
- Coordinate `xi[k]` is native wavelength/frequency/velocity, strictly monotone with declared units and frame; optional spatial WCS, spectral medium (air/vacuum if wavelength), systemic redshift and barycentric convention.
- Explicit sample/voxel validity `V[r,c,k]`, quality/DQ flags, detection support, host membership, representation eligibility and final analysis support **remain different fields**. Do not silently convert NaN to zero.
- Measurement uncertainty: variance, inverse variance or covariance as distinct forms, never simultaneously interpreted without conversion and provenance; declare whether correlations are spectral, spatial or both, and whether the supplied covariance is PSD. Missing uncertainty ≠ zero uncertainty.
- Instrument response: `PSF(λ)`, `LSF(λ,r,c)`, sampling kernels, throughput/filter responses, calibrations and known/unknown flags. Avoid a single unqualified scalar “resolution”.
- Spectral density conventions matter: `f_λ`, `f_ν`, integrated counts and filter-weighted flux densities require Jacobians and normalization contracts.
- Immutability at the scientific level: algorithms do not modify original flux or measurement validity. Supports and transformed features produce new, traceable objects.
- Reproducible summary, input origin and source SHA; original FITS metadata maintained without treating unvalidated header fields as trustworthy.

### `Representation`

Formally, `X = T(F; ξ, V, U, θ_T)`, with declared amplitude, signal subspace, coordinate frame, eligibility and errors. It MUST describe `retains`, `removes`, `invariances`, `conditions` and `validation_domain`. `full_flux` retains measured amplitudes (subject to selected domain); `spectral_shape` may remove a **positive** multiplicative scale `A(F)` and return it separately:

```math
A(aF)=a A(F), \quad X(F)=F/A(F), \quad X(aF)=X(F)\quad (a>0).
```

This identity is conditional on `A` being defined and positive, adequate observed support, and **consistent transformation of noise** when evaluating measured S/N. A genuinely invariant representation cannot simultaneously identify information it discards.

### `SegmentationResult`

- `labels[row,col]`; unassigned label as missing in R, negative sentinel in Python if convenient, **but canonical serialized mask+integer labels independent of language**. Region equality is tested up to relabelling.
- `method_id` and exact dissimilarity/merge criterion; actual vs requested number of regions, graph candidate restrictions, hierarchy and tie diagnostics, coverage and support lineage.
- Cluster membership is a **data-dependent partition**; region IDs are categorical; physical interpretation and spatial connectivity are separate tests.

### `RegionalSpectra`

Per-region native-grid `sum`, `mean`, optionally median/weighted mean, contribution counts, quality flags and covariance contract. Do **not** call an aggregate variance a full covariance matrix. Output tables must be interoperable with pPXF and other downstream fitters without hiding LSF/velocity conventions.

### `ScientificFigure`

Figure is a view of a result and declares native/sky/physical coordinate orientation, region identity semantics, flux units, mask meanings and plotting provenance. A category palette cannot encode ranks unless a separately defined quantity determines ranking. Figures return underlying Matplotlib/ggplot objects and support PDF/SVG; no silent data transformations inside `plot()`.

### `ExperimentalOperator`

Any future transport-aware, kernel or set-valued operator declares domain, observation operators `A_i`, physical constraints, target functional `Q`, validation family, known symmetry or covariance law, and guaranteed versus observed properties. No universal equivariance/identifiability claim from an empirical demo. Experimental work is not part of the stable API.

## 3. Three different segmentation mathematics

### 3.1 Legacy exact `segment()`, method `legacy_l1_wardd2`

For each eligible spaxel, remove the measured spectrum's median; the historical method fills invalid feature entries with zero for distance computation, computes `d(i,j)=sum_k |X_ik-X_jk|`, then applies `Ward.D2` Lance–Williams updates to the resulting dissimilarities. This is **not** equivalent, in general, to optimizing ordinary Euclidean within-cluster SSE. Historical `redshift`, channel-selection and other behavioral conventions require explicit versioned fixtures. “Exact” means all-pairs within this algorithm; it is **not a global optimum** and not necessarily physically accurate.

### 3.2 Sparse `segment_large()`, method `sparse_graph_ward_euclidean`

Create a candidate neighbour graph in a declared feature space. At each allowed merge use the Euclidean centroid Ward increment

```math
\Delta(A,B)=\frac{|A||B|}{|A|+|B|}\,\|\mu_A-\mu_B\|_2^2.
```

Candidate edges and dynamic graph updates mean the chosen hierarchy need not equal unrestricted dense Ward. Record `knn_k`, auto-expansion, disconnected graph status, feature transformations and observed region count. A configured spatial feature/weight is a modeling choice, not evidence that every output region is connected.

### 3.3 Observed-entry hierarchy, method `observed_entry_spatial_ward`

Given validity mask `q_ik∈{0,1}` and features `X_ik`, a within-region criterion over **actually observed entries** is

```math
L(C)=\sum_k\sum_{i\in C:q_{ik}=1}\left(X_{ik}-\bar X_{Ck}^{\rm obs}\right)^2 ,
```

with declared thresholds for shared reliable channels, contributor counts and feature coverage. Candidate merges respect declared spatial contact; missing entries are not fabricated. The algorithm's reported merge increments, admissibility rules and deterministic secondary order must match the R implementation and its validated numerical evidence. **This is not a generic guarantee of missing-not-at-random inference or global optimality.** Document inverse/tie diagnostics and every instrument-specific validated configuration.

### No silent default replacement

An eventual simplified `segment(cube, method=...)` facade may route to distinct algorithms, but cannot change what an existing R call means. Python `segment_legacy` is a parity exercise until gold fixtures pass; `scipy.cluster.hierarchy.linkage` with a precomputed non-Euclidean distance is not a proof of the Ward SSE objective.

## 4. Exhaustive CURRENT R export registry and scientific purpose

The registry below describes current exported **names**, not a claim that every capability has passed scientific validation or that an API rename is already approved. Existing references remain available. Each entry must receive an individual function page with `Question / Intuition / Definition / Assumptions / Returns / Failure modes / Example / Validation / See also`.

| Export | Why a scientist calls it | Required mathematical/scientific warning or contract |
|---|---|---|
| `segment` | Partition a manageable IFU into spectrally coherent regions | Historical median-centred L1 + Ward.D2 distinct from Euclidean SSE; no spatial-connectedness guarantee |
| `segment_large` | Approximate scalable spectral hierarchy via sparse candidate neighbours | Graph-restricted Euclidean Ward; record graph disconnection and transform |
| `choose_ncomp_by_snr` | Choose a partition size subject to region S/N criterion | Requires valid noise model; correlated samples can invalidate naïve S/N |
| `estimate_segment_memory` | Prevent accidental quadratic memory requirements | Memory forecast, not a runtime or allocation guarantee |
| `capivara_representation` | Declare meaning of features fed to an analysis | Must specify retained/removed information, units and validation domain |
| `capivara_amplitude_functional` | Supply a meaningful positive-homogeneous amplitude | Not automatically calibrated or valid on any spectral coordinate |
| `full_flux_representation` | Compare absolute measured spectra | Calibration and distance/uncertainty must be appropriate |
| `spectral_shape_representation` | Compare shape while retaining amplitude separately | Conditional positive-scale invariance; weak/negative amplitude and missing windows fail |
| `manga_sandra_spectral_shape` | Use the specifically evaluated optical MaNGA/Sandra shape definition | S/N≥30, named rest-frame anchor windows; **not** generic all-survey validation |
| `prepare_capivara_representation` | Apply one declared representation to measured data | Preserve voxel validity; publish eligibility causes and dimensional metadata |
| `representation_coverage_qc` | Check adequacy of shared channel coverage | Coverage QC is not a posterior probability of a true physical cluster |
| `spectral_shape_diagnostic` | Diagnose pairwise shape distinguishability under specified noise | First-order/diagonal-Gaussian limits; do not reinterpret as spatial covariance |
| `build_quality_support` | Derive measured validity and flag noncontributing voxels | DQ, variance/ivar and wavelength semantics explicit |
| `build_capivara_support` | Compose eligibility and host/analysis restrictions | Detection evidence, host association and analysis support remain distinct |
| `build_starlet_mask` | Identify foreground support using scale-selected image coefficients | Chosen scales and thresholds, not a certified completeness function |
| `build_adaptive_support` | Alternative data-adaptive support construction | Support choices can induce selection bias; record parameters |
| `capivara_spatial_contract` | Declare pixel geometry, WCS, centre and spatial provenance | Native/sky/physical units, axis order, East-left and PA conventions cannot be assumed |
| `domain_sf` | Construct geometry of the selected spatial analysis domain | Cell geometry, union area and holes from declared coordinates |
| `regions_sf` | Export spatial region polygons | Exact correspondence to pixel IDs; do not interpret geometric area as unbiased physical area |
| `region_adjacency` | Determine which regions share spatial boundaries | Adjacency defined in a named coordinate domain and connectivity convention |
| `convert_wavelength_medium` | Convert air and vacuum wavelength conventions | Explicit physical convention; no silent barycentric/rest-frame correction |
| `emission_lines` | Query available line definitions and wavelength conventions | Rest wavelengths and air/vacuum identity must be traceable |
| `segment_emission_lines` | Segment using selected line measurements/features | Optical-line coverage, sampling, systemic redshift and missing profiles explicit |
| `read_manga_lsf` | Read calibrated MaNGA wavelength-dependent instrumental dispersion | Survey-specific, does not infer intrinsic gas dispersion |
| `select_manga_lsf` | Choose correct pre/post pixel-resolution response for fitting | Fit/sampling mode explicit; invalid combinations fail |
| `segment_kinematics` | Form regional partitions from measured line kinematics or profile shapes | Moments are not automatically intrinsic dispersions; no dynamical assumption |
| `kinematic_models` | Find available physical kinematic models and requirements | A registry is not model evidence |
| `run_kinematic_analysis` | Compare measured velocity field with a specified disc model | Projected geometry/inclination and residual-model assumptions matter |
| `run_manga_bar_model` | Fit an explicitly motivated bisymmetric-bar hypothesis | Requires independent bar geometry; cannot infer a bar from the diagnostic wedge |
| `kinematic_panels` | Extract reproducible model/diagnostic panels | Panels do not constitute independent model validation |
| `prepare_segment_spectra` | Prepare regional measurements and LSF/covariance metadata for fitting | Explicit sparse covariance, response and fit assumptions; optional Python dependency |
| `summarize_cluster_spectra` | Produce region sums, means, medians and propagated quantities | Missing-contributor channels remain missing; complete variance coverage required |
| `reconstruct_cluster_cube` | Produce a representative cube and residual diagnostics | Template choice and masking affect conservation; not new observed data |
| `reconstruct_flux_preserving_cube` | Preserve measured regional sums through supported mean reconstruction | Flux conservation is conditional on defined mask and weighting |
| `plot_cluster` | Inspect categorical segmentation in native/sky/physical coordinates | No ordinal inference from colors; axis directions and WCS explicit |
| `plot_cluster_spectra` | Inspect spectra of assigned regions | Spectrum used, scaling, unit and region link must be visible |
| `theme_capivara` | Use consistent scientific figure typography and spacing | Presentation utility, no scientific estimate |
| `scale_capivara_identity` | Display region identities categorically | Do not turn labels into continuous quantities |
| `scale_capivara_rank` | Display explicitly ranked region property | Rank must be supplied and distinct from ID |
| `scale_capivara_continuous` | Display scalar measured/derived quantity | Units, limits and clipping semantics |
| `scale_capivara_diverging` | Display signed residual/deviation about a reference | Reference zero and symmetric limits when scientifically justified |

This table is an **inventory**: 41 exports. An existing export cannot be removed from this table merely because its documentation becomes less prominent. Stable, expert and experimental are documentation statuses until a deliberate migration plan changes package code.

## 5. Proposed Python scientific facade (not yet implemented)

Python must offer native, transparent functions rather than mimicking R argument order. Tentative design:

```python
from capivara.observation import IFUCube
from capivara.representation import FullFlux, SpectralShape
from capivara.segmentation import segment
from capivara.regions import summarize_regions
from capivara.visualization import plot_diagnostics

cube = IFUCube.from_fits("cube.fits", adapter="manga")
model = SpectralShape(amplitude=..., required_support=...)
prepared = model.transform(cube)
result = segment(prepared, method="observed_entry_spatial_ward",
                 n_regions=20, support=...)
regions = summarize_regions(result, cube=cube, uncertainty="declared")
fig, axes = plot_diagnostics(result, cube=cube, regions=regions)
```

**Illustrative only:** names, defaults, type signatures and adapter availability must be decided and tested in implementation RFCs. No hidden `from_fits` adapter, `FullFlux` class or `segment` stable facade exists because it appears in this proposal.

Avoid implementing the entire API first. Implement one validated path from observation through representation to region products and figures, then grow by verified fixtures.

## 6. Error and epistemic status conventions

Errors must distinguish `InputError` (shape/unit/frame), `UnsupportedAssumption` (e.g., LSF model insufficient), `InsufficientCoverage` (no valid spectral overlap), `UnidentifiableQuantity` (non-recoverable estimand under stated assumptions), and `NumericalFailure` (non-convergence). Exact class taxonomy is an implementation decision, but the **meaning** must be shared across languages.

Return explicit QC: excluded spaxels, zero-contributor wavelengths, missing/invalid variance, graph disconnection, changed frame, unreliable overlap, fitted-model degeneracy, covariance approximation and absence of measured PSF/LSF.

Use `validated numerical implementation` for code-level properties; `observationally validated inference` only with genuine external tests and declared calibration domains. Unit-test success is not a scientific validation claim.

## 7. Research layer: spectro-spatial statistics (candidate, NOT current API)

Physical scene `F(r,u)`, log-wavelength `u=log λ`, local velocity shifts `δ(r)` and instrument mapping `y=A[F]+ε`. At the formal level a transport-corrected derivative

```math
D_\delta = \nabla_r + (\nabla_r\delta)\,\partial_u
```

annihilates a pure shifted template `F(r,u)=f(u-\delta(r))`. This is a known coordinate/gauge-like identity, not novel theory. More importantly, spatial PSF mixing `B` and local shifts `T_\delta` generally do not commute, `[B,T_\delta]\ne0`. Beam smearing and shift-aligned spectral correlations already have substantial prior art.

**Potential new estimand:** an instrumental- and kinematic-response-aware, possibly set-valued functional of cross-line spatial–spectral covariance. Start with precise `Q(F)`, model class `M` and measurement operator `A`; ask if `Q(F)` is determined by observations, or compute admissible range `{Q(F): F∈M, ||AF-y||≤ε}`. Exact observability criteria and data-processing properties are classical; need a new tractable result, approximation guarantee or decisive physical advantage to claim new methodology.

**Promotion gates:** (a) closest prior art incl. spectral correlation functions, beam smearing, forward models and inverse-problem theory; (b) theorem or nontrivial counterexample with explicit assumptions; (c) numerical truth benchmark with two line morphologies, velocity gradient, known PSF/LSF, separate physical change and kinematics, plus independent baseline; (d) real MaNGA / PHANGS-MUSE validation; (e) coverage, stability, scalability. Failure leaves method in research, not package core.

## 8. Specific compatibility and parity fixtures required

1. **Historical exact L1/Ward.D2:** fixed 2×3 and irregular-support spectral cubes; zero-filled missing-channel convention; unambiguous ties and tie-stress case; R reference hierarchy, region co-membership and pixel placement.
2. **Sparse Ward:** known graph, centroid-cost sequence, disconnected candidate graph, auto-expansion and exact small graph agreement where mathematically appropriate; compare Rcpp and Python port on identical explicit neighbor inputs.
3. **Observed-entry:** validity mask, shared-overlap thresholds, semantically named amplitude functional, deliberate partial coverage and deterministic secondary tie ordering.
4. **Physical coordinates:** air/vacuum, observed/rest and systemic redshift, native-channel selection (no interpolation), consistent wavelength-dependent LSF, correct FITS row/col and sky position angle.
5. **Regional products:** observed flux counts, missing contributions, partial variance, spatial covariance toy matrix, all-missing channel, summed flux versus model reconstruction.
6. **Representation invariance:** `X(aF)=X(F)` together with amplitude scaling, valid-domain preservation and transformed noise; failure when amplitude invalid.
7. **Plots:** same region ID rendered consistently across panels, WCS east-left, categorical versus continuous scales, SVG/PDF and accessibility checks.

A comparison is scientific only when source revision, instrument response, feature-selection contract and numerical tolerance are recorded. No default synthetic example proves readiness for all IFU instruments.

## 9. Non-goals and review questions

- No generalized `spin-aware` spectral operator: Doppler shifts, attenuation and filtering have their own transformation structures; imposing the wrong invariance may remove the desired estimand.
- No claim of new physics from the noncommutativity of PSF mixing and spatially variable Doppler alignment.
- Do we promote the observed-entry path as the *recommended* new use pattern while keeping the legacy R default? Recommended: yes, only for already validated domain(s).
- Which physical estimand should guide the first experimental spectro-spatial statistic? Proposed shortlist: line-profile covariance under PSF mixing, line-ratio spatial coherence, or velocity-component multiplicity, with a kill test before code.

**Review criterion:** This RFC defines semantics and inventory; it does not authorize any API break or claim full Python parity.