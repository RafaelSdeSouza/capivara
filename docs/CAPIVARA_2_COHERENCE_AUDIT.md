# CAPIVARA 2.0 coherence audit

Audit date: 10 September 2026. Status: **proposal for review; not approval to refactor or release**.

## 1. Scope, evidence, and decision

CAPIVARA already contains the main components of the proposed scientific package: spectral segmentation, spatial support, regional spectra, structural measurements, line-profile features, native kinematics, and axisymmetric and bisymmetric velocity models. The existing decision in [package_organization.md](package_organization.md) to retain **one user-facing R package** is sound. SpectroPath can remain its mathematical dependency and capivaraPPXF its optional fitting backend. A replacement architecture is not needed.

The immediate obstacle is scientific correctness across component boundaries. Missing flux and variance are combined inconsistently; one workflow mislabels exported FITS channels; another permutes imputation flags before fitting. Structural and dynamical angles do not all use the same convention. The exact and sparse segmentation backends currently use different distance criteria. These are higher priorities than the homepage.

The structural bar detector and the bisymmetric model are **different implementations answering different questions**. The current model workflow does not consume a single, validated result from `detect_bar()`. It can estimate an axis from white-light emission or bar-labelled pixels and can create an angular diagnostic mask. Neither that axis nor a fitted bisymmetric component establishes a photometric bar detection. For iFUN, segment–bar association should use an independently specified, versioned structure mask with mixed-region fractions, not a dynamical model's diagnostic wedge.

### Audited versions and limits

| Component | Version state inspected |
| --- | --- |
| CAPIVARA working tree | `0.3.0`; HEAD `156c60d02e4c014e9895ea7a68f2071ce7400796`, also the remote HEAD at inspection |
| Installed CAPIVARA | `0.3.0`; RemoteSha `467b6ab2f99fddac8ea9cbf1e3daf3eb1f68ada2` |
| SpectroPath checkout | `0.1.0`; HEAD `e1695030c50c09adf5eebdde56336c2342493c25` |
| capivaraPPXF checkout | `0.0.0.9000`; HEAD `3195657b00cf83c2844e74983e7e66444986b533` |
| Installed companion packages | Same displayed versions; no installed RemoteSha recorded |

The installed structural detector body equals the current source body; the installed `run_manga_bar_model()` body does not. Thus the earlier pilot used the current detector implementation, but “both are version 0.3.0” does not identify the model wrapper that ran. Record source SHA and resolved workflow paths before any repeat analysis.

The audit covered all 38 files under `R/`, the sparse-Ward C++ implementation, runtime workflows, package metadata, nine vignettes, both CI workflows, reference coverage, research/prototype organisation, companion interfaces, the original paper, and the rendered website. Research demonstrations were inspected for methods, interfaces, dependencies, and provenance; they were not all executed or independently reproduced. This is a repository-wide coherence audit, not a validation of every scientific estimator.

Evidence labels used below:

- **R — reproduced:** a controlled read-only probe demonstrates the reported behaviour.
- **S — source-confirmed:** the relevant implementation or documentation explicitly establishes the inconsistency.
- **V — validation gap:** a scientific claim needs controlled or observational validation; it is not classified as a demonstrated numerical bug.

The [probe script](CAPIVARA_2_COHERENCE_PROBES.R) and [captured output](CAPIVARA_2_COHERENCE_PROBES.txt) accompany this report. The script evaluates current R source in an environment whose parent is the installed namespace; sparse-Ward uses the installed registered compiled routine. Consequently, these probes are not a clean installation test of the checkout. Their source-level conclusions were checked against the current C++ and R implementations. The pPXF probe uses saved pilot wavelength arrays, not new population fits. No algorithm, public API, homepage, CI, or GitHub metadata was changed for this audit.

Pre-existing uncommitted README, vignette, pkgdown, generated-site, and link-checker edits were preserved. Existing iFUN audits were treated as hypotheses and requirements, not as unquestioned specifications. The local requirement/evidence sources were `results/bars_v1/provenance/capivara_package_audit.md`, `capivara_integration_design.md`, `capivara_api_audit.md`, and `results/bars_v1/bar_detection_validation/bar_detector_validation.md` in the separate iFUN workspace. These are not shipped CAPIVARA APIs or package validation tests.

## 2. Capability and API map

The complete [capability/API map](CAPIVARA_2_CAPABILITY_MAP.md) records implementation, export, object, input, survey and coordinate assumptions, documentation, tests, example, and provenance for every major capability. It also accounts for every current export and all 38 R source files.

“Production public API” in that inventory identifies the established intended package surface; it does **not** certify scientific readiness. The established spectral products remain subject to the defects below. New line-profile and model APIs should be explicitly experimental until their input contracts and validation gates pass.

The most substantial hidden capabilities are `score_structures()`, `threshold_structures()`, `catalogue_structures()`, `detect_bar()`, `detect_ring()`, and `segment_structures()`. They are implemented but not exported. Low-level readers, geometry routines, velocity fits, and residual diagnostics are also internal. SpectroPath is accessible through the native workflow but lacks an adequate scientific entry page. Two actual exports, `emission_lines()` and `segment_emission_lines()`, are absent from the built reference documentation.

## 3. Contradictions, defects, and duplicated implementations

Priorities describe the consequences of using the affected path: **P1** can corrupt scientific products or their interpretation; **P2** obstructs reproducibility, portability, or correct use; **P3** is maintenance or presentation debt. None of the following should be hidden behind a generic “experimental” label instead of being corrected.

### 3.1 Reproduced correctness defects

| ID | Priority / evidence | Finding and consequence | Required response |
| --- | --- | --- | --- |
| C01 | P1 / R | `R/summarize_cluster_spectra.R:74–95` forms weights from variance without jointly masking flux. Flux `(10, NA)` with variance `(1,1)` yields weighted mean 5 and variance 0.5, rather than 10 and 1 for the one usable datum. All-missing channels return summed flux and variance zero. | Use one explicit validity rule for each estimand; return missing for no contribution; report contributing counts and incomplete variance. Test S/N selection against the same rule. |
| C02 | P1 / R | `inst/extdata/kinematics/native_kinematics_workflow.R:885–919` names slots 3:8 as six line maps, but the line maps occupy slots 4:9. The segment map is labelled log flux, flux is labelled velocity, and so on; slot 9 is unnamed. The exported channel catalogue is scientifically wrong. | Construct a named list once and assert names/order before writing; round-trip labelled sentinel arrays through FITS and its catalogue. |
| C03 | P1 / R | `native_bisymmetric_workflow.R:108–111` attaches `as.vector()` masks to a spaxel table ordered with x varying fastest. R's matrix vectorisation instead advances row/y fastest. The `imputed` flag and 0.35 fitting weight reach the wrong pixels. This is not restricted to rectangular images. | Index every map by `[cbind(spaxels$y, spaxels$x)]` or establish a single pixel key; test rectangular and square arrays. |
| C04 | P1 / R | `R/structural_awareness.R:294` measures a component axis with `atan2(y,x)`, whereas the bar ellipse at line 642 uses `atan2(x,y)`. A horizontal fixture returns −180° and 90° respectively. Ring geometry also uses the first convention. These cannot share `pa_deg`. | Define a named image-angle convention and convert all components, rings, bar axes, and orientation features consistently. Preserve legacy conversion explicitly. |
| C05 | P1 / R | `estimate_bar_geometry()` in `R/kinematics_disc_model.R:147–174` uses `atan2(Yd,X)`, but NIRVANA-convention `deproject_coordinates()` uses `atan2(-Yd,X)`. A 30° bar-labelled fallback corresponds to θ=330° in the model. The white-light estimator follows the model convention and is a separate path. | Make fallback estimation convention-aware; test inferred angle → projected vector → recovered angle. Do not reverse every plotted vector indiscriminately. |
| C06 | P1 / R | `R/kinematics_plot_capivara_kinematics.R:100` uses the row count for the y coordinate of a counter-clockwise 90° point rotation. A 2×3 array's point (3,1) maps to (1,0), while the rotated image puts it at (1,1). | Use the correct dimension and one shared transform for image, mask, centre, and axis. The inverse disc-plane projection itself is not implicated by this fixture. |
| C07 | P1 / R | `R/structural_awareness.R:1716–1795` appends structure features as pseudo-wavelength channels and leaves that augmented input in `original_cube`. Eight spectral channels plus three feature maps produce an eleven-channel “spectrum” downstream. | Keep feature matrices separate from observed spectra, or restore the real cube and metadata before returning; assert wavelength–spectrum dimensions. |
| C08 | P1 / R, S | `R/segment_emission_lines.R` replaces missing line-window values by zero before feature validity is assessed. An entirely missing input spaxel receives a segment label. In the non-sampled path, the top-level wavelength metadata are also lost, although `original_cube` is restored. | Propagate original spectral eligibility independently of numerical feature imputation; preserve metadata in both branches. |
| C09 | P1 / R | `R/starlet_layer.R` limits reflected padding to one image extent but later convolution assumes the requested padding exists. On 24×24 input with J=5, each of the first four scales has 576 finite pixels; the fifth scale and coarse image each have only 36. | Support repeated reflection or reject incompatible scale/image combinations. Test small and rectangular cutouts; do not interpret padding failures as faint structure. |
| C10 | P1 / R, S | The exact backend requests Manhattan distance, then `ward.D2`; the sparse C++ backend uses Euclidean Ward variance increments. A fully connected sparse graph matches Euclidean Ward but disagrees with Manhattan Ward on a small fixture. “Exact versus approximate” is therefore incomplete. | Decide and version the intended scientific distance, then compare graph approximation under matched features and distance. Preserve a reproducible legacy mode rather than silently changing old analyses. |
| C11 | P1 / R | capivaraPPXF computes an initial logarithmic galaxy grid and then rebins again using its returned velocity scale, while retaining the first grid. For pilot redshifts 0.01837224 and 0.02453957, lengths change 2647→2646 and 2664→2663. The pattern occurs in `population.R`, `emission.R`, and `diagnostics.R`. | Use each returned rebinned spectrum with its corresponding wavelength grid and enforce one final galaxy/noise/mask grid. Templates retain their own adequate wavelength coverage; do not crop templates arbitrarily to the galaxy length. |

The current exact metric also conflicts with the Euclidean distance stated in the original CAPIVARA method. Resolving this requires a reproducible comparison, not a documentation-only correction. [Original paper](https://doi.org/10.1093/mnras/staf688), [accessible paper text](https://arxiv.org/html/2410.21962v2).

### 3.2 Additional semantic and scientific risks

**C12 — P1/S: sampling changes the emission-feature metric.** Below `max_pixels`, `segment_emission_lines()` passes `scale_fn=NULL` into `segment_large()`, selecting its row-wise MAD scaler. The sampled branch fits the column-scaled features directly. Crossing the default 30,000-pixel threshold therefore changes normalisation as well as using sampled projection; many `...` controls are only honoured in the non-sampled branch. Make preprocessing identical, list supported controls, and identify both approximation and scaling in provenance.

**C13 — P1/S,V: bar-candidate acceptance is not calibrated.** `detect_bar()` accepts `profile_like || body_like`. The body branch can accept an elongated central component without passing the full ellipticity/PA/profile gate. Its confidence is a heuristic score, not a probability. The Ferrer fit is to normalised structure scores, not a calibrated surface-brightness decomposition, and its support criterion does not require successful convergence. Return separate profile/body diagnostics and an explicit candidate/accepted/rejected status; do not advertise a calibrated detector until negative controls are evaluated. A body mask and the broad ellipse need separate geometry fields.

**C14 — P1/S,V: native line measurements are proxies with incomplete units and error handling.** The script sums continuum-subtracted samples without wavelength-bin widths, and estimates moments with positive, locally selected line weights. Its “flux” is not generally a wavelength-integrated line flux; dispersion is an observed moment without an LSF correction. `h3_proxy` and `h4_proxy` are not fitted Gauss–Hermite coefficients. Line rest wavelengths and aliases are duplicated between the public emission code and workflow; air/vacuum conventions are not represented, while the pPXF gas path explicitly uses vacuum templates. Define wavelength medium, units, systemic reference, continuum, integration measure, LSF and uncertainty before presenting physical measurements. MaNGA-specific calibration belongs in its adapter, not in a generic moment routine. [SDSS data model](https://www.sdss4.org/dr17/manga/manga-data/data-model/).

**C15 — P1/S,V: model assumptions are too easy to mistake for measured geometry.** Native fitting permits a default inclination of 60°. A numerical `fit_status="ok"` can coexist with placeholder geometry. The present bisymmetric model is a regularised linear fit at supplied/estimated geometry, not a full NIRVANA inference implementation: it lacks the corresponding joint Bayesian inference and PSF-aware forward modelling. Strong smoothing and weights, including arbitrary imputation downweighting, require injection/recovery calibration before very small non-circular amplitudes can be interpreted physically. [NIRVANA barred-galaxy study](https://arxiv.org/abs/2407.11908).

**C16 — P1/S: bar masks and residuals change meaning between paths.** The optional native bar mask is an angular wedge extending through most of the radial coverage, not `detect_bar()`'s measured footprint and not a likelihood prior. It classifies diagnostic spaxels; it does not itself constrain the model as a photometric prior. `f_bar` is the fraction of residual squared power inside the designated bar region, not its area or luminosity fraction. With no bar mask, zero selected pixels is not evidence for an absent bar. In the bisymmetric result, `v_disc` stores the full bisymmetric prediction and `v_resid` stores residuals after that model. Consequently, `Q_kin` is not automatically a disc-subtracted diagnostic across model choices. Use names tied to the actual prediction and residual reference.

**C17 — P1/S: the selected kinematic segmentation need not reach model labels.** The native bisymmetric script reads `native$kinematic_aware` even when the requested segmentation mode is `path_signature`. Fitting is spaxel-wise, not a fit to the regional spectra or regional velocities. State that distinction and use the selected label map for any region-level diagnostics. No current implementation should be described as fitting a common set of spectral, kinematic and pPXF regions unless those IDs were explicitly linked.

**C18 — P1/S,V: fitting success and scientific QC are conflated.** capivaraPPXF carries input variance but its specialised fitters use a scalar noise estimate instead. Missing samples are filled numerically rather than excluded with a preserved measurement mask. `fit_ok` records an exception-free call; parameter-bound solutions can still be marked successful. Emission QC combines fit/positive-flux/ratio criteria without a full uncertainty-based line-detection assessment. Preserve optimisation status separately from bounds, coverage, residual, noise and science-use flags. Returned gas coefficients remain tied to the spectrum normalisation unless physical units are restored; population `mean_log_age` is a light-weighted logarithmic age, not a linear mean age, and an M/L alone is not a stellar mass.

**C19 — P2/S: feature choice, support, and S/N are not fully separable.** Feature wavelength subsetting occurs before target-S/N evaluation, so a genuinely independent S/N window outside those features is unavailable. Variance arrays without wavelength metadata can be subset on channel numbers while flux is subset in physical wavelength. The no-variance S/N proxy treats flux like counts; its scale is not generally meaningful for calibrated IFU flux densities. Missing features are centred then zero-filled, not compared using a validated overlap-aware distance. Exact and sparse paths also differ in valid-spaxel rules. Document the estimand and make measurement validity independent of clustering preprocessing.

**C20 — P2/S: generic algorithms exist, but a generic input contract does not.** `.wavelength_axis()` uses FITSio `axDat` or channel indices; it does not honour an independent explicit wavelength vector. Wavelength subsetting can reconstruct a linear axis using a median spacing even if the original sampling was non-linear. `plot_cluster_spectra()` assumes `axDat` and labels Å. Core segmentation accepts arrays, but native kinematics is a filename/HDU workflow with MaNGA metadata assumptions and hard-coded Hα provenance strings even for other selected lines. MUSE/MAGNUM demonstrations establish prepared-array use, not a validated MUSE DATA/STAT reader. [ESO MUSE product overview](https://www.eso.org/observing/dfo/quality/MUSE/pipeline/pipe_gen.html).

**C21 — P2/S: masks and IDs need explicit meanings.** A starlet support mask is a selected analysis footprint, not a probability or bar mask. Native line-flux support is restricted by the starlet footprint, so it cannot recover line-only emission outside it. Cleaning to a main connected galaxy can intentionally discard companions or tidal structures; that choice is consequential for environment studies. `original_cube` may contain the supported/masked copy, not the full original data. Segment IDs can be non-contiguous after rejection and must be joined by ID, not vector position. The older RDS reader expects `segmentation_map` rather than core `cluster_map`; it rounds labels rather than requiring valid integer IDs. Observation ID, galaxy ID and segment ID must not be interchangeable.

**C22 — P2/S: configuration and code provenance can leak across runs.** `.capivara_workflow_file()` checks the installed package before an explicit checkout. Wrappers set selected environment variables and restore those keys, but settings they do not set can be inherited from an existing session. Several geometry controls exist only through environment variables. A returned NA control can therefore differ from the effective setting. Prefer explicit ordinary function arguments; immediately record all effective controls and resolved source paths while compatibility wrappers remain. Silent `try()` around FITS writing and reusable output prefixes can conceal failed or overwritten exports; return a verified output manifest.

**C23 — P2/S: object and display semantics still diverge.** `plot_cluster()` maps Row to x and Col to y; other plots map column to x and row to y. Display reversal, transpose and rotation are mixed with inferred axes in different paths. An axial PA is modulo 180°; an eigenvector's sign is not an observed receding direction. `axis_ratio` should consistently mean b/a; ellipticity should mean 1−b/a. An image axis becomes a sky PA only through the local celestial WCS, not by relabelling an image angle.

### 3.3 Documentation and duplication

| ID | Evidence | Contradiction / duplication | Disposition |
| --- | --- | --- | --- |
| D01 | S | `docs/public_api_audit.md` says `bar_phi_deg` is required; current source defaults it to NA and can estimate an axis. `kinematic_models()` calls itself a registry of installed modules but returns two fixed rows with a hard-coded dispatch branch. | Document actual angle-source policy and fixed model catalogue. Keep future module discovery planned rather than implied. |
| D02 | R, S | The historical kinematics package repeats 59 same-named current functions: 52 bodies match and seven diverge. Its README still teaches a second installation. | Freeze as explicitly non-production history; do not repair both copies. Retain provenance and the migration comparison. |
| D03 | S | `run_one_galaxy()`/batch/DAP orchestration coexists with native script orchestration; reports still mention an arctan baseline although the runner now uses piecewise fits. The arctan primitive also remains implemented. | Preserve useful model primitives and historical reproducibility; retire redundant orchestration after parity tests, not before. |
| D04 | S | Public `capivaraPPXF::fit_ppxf()` always stops as an unimplemented scaffold, while specialised population and emission fitters exist. Preparation, rebinning and noise code are duplicated across three fitters. | One common preparation function and a thin CAPIVARA fitting entry point, after C11/C18. Deprecate or implement the scaffold honestly. |
| D05 | S | `classify_structures()` aliases a neutral catalogue rather than classifying physical components. `detect_ring()` can return a candidate mask when `ring_like=FALSE`, unlike the bar default. | Prefer `catalogue_structures()`; standardise candidate and accepted fields before export. |
| D06 | S | README identity is narrower than the implemented package, while some vignette claims are stronger than the code: “exact vs approximation”, bar footprint/prior, and portable cube input. | Rewrite claims after scientific decisions; keep observed quantities, proxies and model assumptions distinct. |
| D07 | S | `DESCRIPTION` omits much of the actual line-profile/kinematic scope. NEWS has no post-0.3 history despite structural, kinematic and emission additions; “added segment_large” appears under both 0.2 and 0.3. | Reconstruct changes from commits; distinguish the old block-medoid backend from sparse Ward without inventing release history. |
| D08 | S | `segment()` documentation describes a default base `scale` transformation, but the default is `median_scale`, which subtracts a median. “Median-normalised” can therefore be ambiguous. | Name centring, division and robust rescaling separately; list the actual ordered transformations in results. |

### 3.4 What the existing pilots do and do not show

The five positive-galaxy bar pilot and its perturbations support stability on those selected examples. Reported maximum PA changes were 0.29–3.38°, fractional radius changes 0.101–0.319, and median mask IoUs 0.958–1. Those are not completeness or purity estimates. The 90° checks used precomputed structure-score maps; they did not validate the entire image-to-score pipeline or every display transform.

The previous smooth-ellipse negative example supplied synthetic score maps directly. It demonstrates that the body branch can bypass the profile gate conditional on those maps; it does **not** measure the full detector's false-positive rate on a smooth raw image or cube. The next benchmark must include independently labelled real unbarred galaxies and raw-image/cube controls.

The pPXF pilot's 36 failed regions are consistent with the now-reproduced rebin regression. Other saved fits reached parameter bounds; none was refitted in this audit. A low non-circular amplitude under fixed 60° inclination and regularisation is not a demonstration that the non-circular model is validated or that a galaxy lacks streaming motions.

An earlier suspicion that `FITSio::readFITS(..., maxLines=1)` necessarily fails is not supported: a short-header fixture succeeded. The real-file connection warnings remain unresolved. Do not report a proven reader failure or leak without a reproducer.

## 4. Proposed scientific architecture

### Retain the existing components

Keep `segment()` and `segment_large()`, both support builders, regional summaries, both reconstruction products, and lightweight plotting. Retain the distinction between measurement-only `segment_kinematics()` and model-selecting `run_kinematic_analysis()`. Preserve `run_manga_bar_model()` as a survey convenience wrapper, with its experimental and geometry requirements explicit. Keep SpectroPath as a separately tested mathematical package. Keep Python pPXF and template installation optional and user-controlled.

Do not implement the earlier proposed seven mandatory object classes, a workflow DSL, or a large new dispatcher. The existing list results can acquire validated fields and a small number of optional S3 classes without forcing ordinary users to manage a new software system.

### The minimum common data contract

Validate input at the boundary, while retaining a compatibility adapter for current `imDat`/`axDat` inputs:

- Observed flux has declared dimensions `[row, column, wavelength]`; explicit wavelength samples have a unit, air/vacuum convention, and observed/rest-frame state.
- Variance, inverse variance and quality flags are distinct. Convert inverse variance only where finite and positive; retain the original validity reasons. A missing variance cannot silently become a zero error.
- WCS, pixel scale, flux unit, LSF and PSF metadata are carried when available. Record unknown quantities rather than inventing survey defaults.
- `observation_id` identifies a cube; optional `galaxy_id` identifies the physical target. Survey identifiers are metadata. No ID is inferred by silently rounding a numeric value.
- Preserve measured arrays and measurement masks. Feature imputation, analysis support and display transforms are separate fields.
- Provenance records package/dependency versions and SHAs where available, input identity/checksum, effective controls, wavelength selection, scaler, distance, graph settings, seed, and adapter/version.

A common result need only retain its input reference, labelled map, ordered `segment_ids`, relevant products, diagnostics and provenance. Regional spectra must always index genuine wavelength samples, never augmented features. Optional columns or nested lists should carry uncertainty and QC without breaking existing consumers.

### Geometry contract

Use `x_pixel=column` and `y_pixel=row` at native array pixel centres. State the direction of each image axis; plotting alone decides whether y increases upward or downward on screen. Define `pa_image_deg` as an axial angle from +y toward +x modulo 180°. Define `pa_sky_deg` separately as a celestial east-of-north angle obtained from WCS. In the existing NIRVANA convention, `phi_bar_disc_deg` is measured using the same signed θ returned by deprojection; it is not an image or sky PA.

Retain `axis_ratio=b/a` and `ellipticity=1-b/a`. Do not infer inclination from b/a without recording the assumed intrinsic thickness and geometry model. Distinguish an unoriented structural axis from a receding kinematic direction. Every conversion must have inverse and rotation tests, including parity-changing WCS and non-square fields.

### Structures and segment membership

Promote neutral score/catalogue functions only after their coordinate and result contracts are tested. Export bar/ring detection as explicitly experimental until negative-control calibration supports a stronger status. Keep the candidate score, candidate mask, accepted mask, selected geometry, reasons and calibration identifier separate. External imaging or expert-vetted masks remain valid inputs; CAPIVARA need not manufacture its own bar detection when suitable independent measurements exist.

Add one small segment–structure association function. For each segment it should report the fraction of supported area inside each accepted structure, optional flux fractions in a declared band, overlap/coverage counts, and an unresolved or mixed state. Area fraction and flux fraction answer different questions and need different names. A binary bar label is a downstream rule with an explicit threshold, not a property guaranteed by segmentation.

For iFUN, use the original MaNGA spectra, variance/IVAR, wavelength and masks stored inside a MEGACUBE only after verifying their extension identities. Run a fresh CAPIVARA segmentation. Existing MEGACUBE bin IDs, population fits and structural labels are not CAPIVARA outputs. Environment catalogues and the scientific comparison between barred/unbarred regions remain in iFUN; the reusable overlap calculation belongs in CAPIVARA.

### Line profiles, models, and fitting

Separate native line extraction from filesystem orchestration. Return conventional measured maps and SpectroPath features as distinct products sharing the same spaxel key and quality masks. SpectroPath represents ordered line-profile information; its features are not interchangeable with fitted velocity, dispersion, or Gauss–Hermite coefficients.

Let the existing model selector consume those measured maps with explicit geometry and uncertainties. Keep the current two-model catalogue simple; extend it with stability, required measurements, geometry-source policy and output definitions. A model result should identify `velocity_model`, `velocity_residual`, and their model reference, not overload `v_disc`. Report numerical status, geometry adequacy, measurement coverage, regularisation sensitivity and science-use flags separately. Model comparison must use a common valid sample and error model; lower residual RMS alone does not select a bar.

Wrap the specialised capivaraPPXF implementation behind one optional regional-fitting entry point only after common preparation and QC are corrected. Retain backend-native details, templates, LSF assumptions, wavelength ranges, normalisation and fit configuration. Join fits back to regions by validated IDs. Do not install Python, download SPS models, or accept external software terms implicitly when CAPIVARA loads or runs ordinary segmentation.

## 5. Proposed pkgdown navigation and visual policy

Use a small top-level navigation: **Get started**, **Science**, **IFU data**, **Validation & papers**, **Reference**, and the repository link. Science is a dropdown, not a grid of promotional cards.

| Navigation / article group | Scientific question | Existing material to retain or revise |
| --- | --- | --- |
| Get started | How do I obtain a labelled map and interpretable regional spectra? | `getting-started`; one executable wavelength-aware synthetic cube; status and citation links |
| IFU data | What must my cube contain, and which measurements are retained? | New common contract; MaNGA adapter; prepared MUSE arrays, then validated MUSE adapter; custom FITS/arrays; existing FITS/DS9 guide |
| Science → Spectral structure | Which spectra are grouped, and what does the grouping preserve? | Backend choice, support masks, regional spectra, reconstruction; explicit distance and missing-data semantics |
| Science → Galaxy structure | Which measured structures overlap a region? | New neutral structure, experimental bar/ring, and segment-association pages, created only with implemented APIs |
| Science → Line profiles and kinematics | What information is in a resolved emission-line profile? | Revise kinematics; separate native measurements, SpectroPath, and kinematic segmentation |
| Science → Dynamical models | Which velocity component is measured relative to an axisymmetric comparison? | Axisymmetric and bisymmetric articles, explicit geometry, matched model comparison and limitations |
| Science → Physical inference | How are regional spectra fitted and results returned to the cube? | Optional pPXF page with uncertainty/QC and template assumptions; no fictitious implemented engines |
| Validation & papers | Which claims have been tested, on what data, and how are papers reproduced? | Tests/controlled benchmarks/limitations; existing examples and paper-reproduction article |
| Reference | Which function answers each question? | Scientific groups covering all exports; visible experimental badges and minimal runnable examples |

Retain existing URLs where possible and add redirects for moved articles. Reference groups should separate support, spectral segmentation, regional products, line measurements, kinematic segmentation, models, and plotting. Do not publish an empty MUSE or fitting page as evidence that an adapter or facade already exists.

Preserve the restrained navy/light/dark appearance. Within the first screen, place a short scientific opening, current version/status, start/data/analysis links and citation. Move the large five-galaxy mosaic below those decisions or reduce it to one legible real example. Use a few scientific figures: a labelled IFU with regional spectra; observed velocity, baseline prediction and residual with units; and a line-profile example distinguishing conventional and path features. Synthetic or smoothed illustrations must be labelled as such. A bar figure must show which geometry was supplied and which was inferred.

### Baseline rendered-site QA

The deployed homepage, all nine principal articles and the function index were inspected structurally in the browser, including mobile width 390 px and desktop layouts (initial homepage 1280 px; final article pass 1422 px as reported by the browser). Representative pages were visually inspected in light/dark mode. The article/index DOM checks showed no page-wide horizontal overflow at these widths. Some images were initially still loading; the example mosaic and backend-comparison image both loaded on reinspection. Mathematical notation is rendered as native MathML; absence of MathJax containers is not a rendering failure.

Confirmed presentation problems are the oversized opening mosaic, small scientific labels on mobile, low-contrast inline code/API links in dark mode, and missing scientific figures on the kinematics/model articles. The current homepage does not expose the implemented analysis paths and stability distinctions within its first screen. The function index lacks the two new emission exports.

The local link checker passed 37 HTML and 48 text files. It checks internal targets/anchors, page structure and alt text, but is not an image-file checker, external-link validator, executable-example test, or scientific-content validator. No site changes were made here; a **fresh local render and full desktop/mobile visual QA remain required after approved changes**, not certified by this baseline.

## 6. Proposed homepage opening

The following opening is proposed text, not a claim that the correction and release gates have already passed:

> CAPIVARA is an R package for spatially resolved analysis of integral-field spectroscopy. It groups spaxels by spectral similarity, measures regional spectra, and connects those regions to galaxy structure and emission-line kinematics. The scientific task is to relate variation across a galaxy to the information in its spectra, while keeping the measured flux, spatial support and assumptions of each analysis explicit.
>
> The core segmentation and regional-spectrum routines operate on prepared IFU arrays. MaNGA workflows are included; MUSE and other instruments require input with correctly specified wavelength sampling, masks, uncertainties and spatial coordinates. Spectral segmentation can be followed by structural measurements, line-profile analysis or velocity modelling. SpectroPath supplies a distinct description of ordered line-profile structure alongside conventional kinematic measurements.
>
> The established spectral APIs are accompanied by experimental structure and dynamical analyses. Bar geometry inferred from an image is separate from a bisymmetric model fitted to a velocity field. Stellar-population and emission-line fitting are downstream analyses supported by the optional capivaraPPXF backend, with their own template and uncertainty assumptions.
>
> Start with the synthetic-cube example, then choose the data and analysis guide relevant to your observation. Consult the validation pages before interpreting experimental results. Cite the original CAPIVARA paper for the published segmentation method and record the software version used for subsequent developments.

Keep README.qmd shorter than the website: this scientific scope, a minimal executable example, installation/dependency boundary, status, citation, and links. Move extensive synthetic-data construction into the article/helper. Derive README.md from README.qmd, but let pkgdown carry the fuller explanation rather than reproducing the entire README literally.

## 7. Version, release, and metadata plan

The checkout and installed package expose materially different implementations as `0.3.0`. The first approved maintenance change should mark the current branch **`0.3.0.9000`**, with an unreleased NEWS section accounting for post-0.3 support, structural, line-profile and kinematic work. No new release date should be invented. The inspected repository has no tags; preserve the known SHA mapping and establish annotated release tags prospectively.

Recommend a **0.4.0 correctness/development release**, not an immediate 2.0 relabelling. It should include the reproduced regression tests and fixes, explicit experimental status, complete reference coverage, matched version metadata and migration notes. Distance/coordinate changes need named compatibility controls or an explicit breaking-change migration, plus saved legacy fixtures.

Reserve **2.0.0** for a tested common data/geometry contract, coherent installed workflows, validated MaNGA and MUSE fixtures, propagated uncertainty/QC, bar controls and model-recovery evidence, and the optional fitting interface. Scientific calibration may leave individual modules experimental in 2.0, but the boundaries must be honest. A change of major version does not itself validate bar detection.

For each release: generate roxygen → build/install in a clean library → run package and integration tests → build vignettes/pkgdown → check links/images and execute examples → visually inspect → tag and archive the tested software. Record dependency versions; use reproducible environments for published analyses without unnecessarily pinning every ordinary user installation.

Proposed DESCRIPTION title: **Spatially Resolved Analysis of Integral-Field Spectroscopy**. Its description should state spectral segmentation, regional flux/variance products, spatial support, ordered line-profile features and kinematic comparisons, marking structural/model methods experimental where appropriate. Do not describe current zero-filling as a validated missing-data-aware metric.

The inspected GitHub About panel supplied no description, website or topics. Proposed values, **not applied**:

- Description: “R tools for IFU spectral segmentation, regional spectra, galaxy structure, line profiles and kinematic modelling.”
- Website: [CAPIVARA documentation](https://rafaelsdesouza.com.br/capivara/).
- Topics: `r`, `astronomy`, `integral-field-spectroscopy`, `spectral-segmentation`, `galaxy-kinematics`, `spectral-fitting`.
- Citation metadata: preserve `inst/CITATION` and the original paper; add a consistent `CITATION.cff` with software version/repository and preferred paper citation. Add a software DOI only if a real archived release is created.
- Badges: license, paper/citation, documentation and R CMD check only when their targets exist and report real status. Avoid a “production-ready” badge.

The original 2025 paper supports the published spectral-segmentation method and its MaNGA demonstration, including downstream spectral synthesis; it does not validate the later structural detector, SpectroPath integration or bisymmetric implementation. Keep the published method and new development evidence visibly separate. [Paper](https://doi.org/10.1093/mnras/staf688), [repository](https://github.com/RafaelSdeSouza/capivara).

## 8. Test and scientific-validation matrix

CAPIVARA has no `tests/` directory. Its CI publishes documentation and checks README rendering; neither workflow is an R CMD check gate. The existing local check log reports **0 errors and 3 warnings** under R 4.5.2/macOS ARM with `--no-manual --no-build-vignettes`: two undocumented exports, stale built vignettes, and missing built output for seven vignettes. That is not a full vignette rebuild or a scientific validation result.

SpectroPath has a small wrapper/parity and feature test set, useful but not a CAPIVARA integration benchmark. capivaraPPXF's existing tests cover preparation/setup/I/O, not sufficient end-to-end fit recovery. The probe script supplied here is audit evidence, not a substitute for a maintained testthat suite.

| Layer / proposed test group | Cases and assertions | Scientific acceptance / scope |
| --- | --- | --- |
| Core object validation | Dimension mismatch, explicit/non-linear wavelengths, units, air/vacuum, variance vs IVAR, masks, unknown metadata, non-contiguous/duplicate IDs | Invalid combinations fail before analysis; no silent axis or unit invention |
| Spectral segmentation | Known small partitions; exact distance and legacy mode; full-graph sparse equality under matched Euclidean features; disconnected graph and infeasible Ncomp; seed, valid-mask and scaler provenance | Separate algorithm correctness from sparse approximation; benchmark stability against graph size |
| Emission segmentation | Sampled/non-sampled common preprocessing; threshold boundary; every supported control; all-missing spaxels; restored wavelength metadata | Crossing max_pixels changes only the documented approximation |
| Support/starlets | Constant/noise/point/extended images; tiny and rectangular arrays; large J; rotations; edge support; disconnected companions | No manufactured missing pixels; report support completeness on faint/extended controlled sources |
| Regional flux/variance/SNR | C01; all-missing/partial coverage; negative flux; zero/non-finite variance; weighted mean analytic cases; covariance inflation; independent S/N window | Flux and its uncertainty refer to the same contributors; covariance assumptions are explicit |
| Reconstruction/export | Sum-preserving and representative products; real wavelength dimensions after feature augmentation; ID mapping; WCS and named FITS channel round-trip | No feature channel in a spectrum; exported arrays retain their physical identities |
| Coordinates and plots | Horizontal/vertical/oblique axes; centre offsets; 90°/180° rotations; every display transform; non-square arrays; WCS handedness and sky/image conversion | Same physical line overlays the same pixels; axial and directed angles remain distinct |
| Structure/bar/ring geometry | Candidate/accepted schema; ellipse/ring limits; Ferrer convergence; both profile/body branches; raw image→scores→mask; rotation equivariance | No branch bypass hidden by a single confidence value |
| Negative bar controls | Unbarred discs, elliptical bulges, spiral arms, rings, lopsidedness, mergers, companions, noise and PSF/inclination sweeps; independently classified real negatives | Predeclare detection definition and thresholds; report completeness/purity by resolution and inclination |
| Segment–structure overlap | Inside/outside/mixed regions, partial coverage, multiple masks, declared flux band, absent/unaccepted bar | Area and flux fractions verified analytically; unclassified does not become “non-bar” |
| SpectroPath mathematics/integration | Ordered/reversed profiles; analytic paths; irregular grids, units, translation/scaling, noise, masked samples; R/Python or trusted analytic parity where available | Mathematical invariants in SpectroPath; observational robustness and feature selection in CAPIVARA |
| Native line maps | Injected Gaussian/asymmetric/blended profiles, absorption continuum, no line, missing channels, non-uniform wavelength bins, rest/observed frames, air/vacuum, variable LSF | Quantify flux/centroid/observed-width bias; fitted GH parameters cannot be asserted from moment proxies |
| Model catalogue/configuration | Every advertised model runs; unknown model rejection; explicit required geometry; clean vs contaminated environment; installed vs checkout provenance | Effective controls equal reported controls; no silent inclination or axis fallback in science mode |
| Imputation and weights | C03 index fixtures; measured-only comparison; known-error cases; no-data areas; covariance and imputed pixels | Primary scientific likelihood uses measured data or a justified imputation model; arbitrary weights labelled exploratory |
| Axisymmetric/bisymmetric recovery | Zero-bar controls and injected Vt/V2r/V2t; varying geometry, radial sampling, noise, PSF and regularisation; matched fit footprint | Report bias, degeneracy and false-positive rates; assess interval coverage only if intervals are actually implemented |
| pPXF grids/templates | C11 redshifts plus synthetic boundaries; paired galaxy/noise/mask grids; template coverage and LSF; masked gaps; independent template wavelength grid | No off-by-one or arbitrary cropping; controlled recovery before pilot refits |
| pPXF QC and units | Bound hits, optimiser failure, poor coverage, wrong/no variance, low line S/N, normalisation restoration, known populations and gas lines | Fit status, QC, uncertainty and science eligibility remain separate; no physical-flux or mass claims from unconverted coefficients |
| MaNGA fixture | Small named-extension FLUX/IVAR/MASK/WAVE fixture and metadata; DAP distinct from DRP; same original arrays in a verified MEGACUBE wrapper | Fresh segmentation and consistent WCS/units; no import of old binning as new output |
| MUSE fixture | Small DATA/STAT fixture; header spectral conventions; non-square footprint; variance, units and LSF; compare with equivalent prepared array | Demonstrated adapter equivalence, not merely ability to accept an R array |
| Synthetic IFU integration | Support → spectral/line segmentation → regional spectra → structure overlap → measured kinematics → both models → optional fitting/export | Shared IDs, coordinates and provenance; known latent components used only for benchmark assessment |
| Documentation/release | All exports have Rd/index entries; executable synthetic examples; clean vignette/site rebuild; internal/external links; image existence/alt; desktop/mobile/light/dark | No stale version or claims stronger than tested implementation |

Add testthat edition 3 and a standard `.github/workflows/R-CMD-check.yaml` after approval. Run PR checks on Linux, Windows and macOS release R, with an additional Linux R-devel lane as resources permit. Test base-R and optional torch distance paths separately. Keep deterministic tiny fixtures in the repository; keep survey-scale data and expensive recovery grids out of ordinary checks. Build at least one CI lane with vignettes enabled.

Place mathematical unit tests in SpectroPath and wrapper/measurement tests in CAPIVARA. Put Python rebin and fitting regressions in capivaraPPXF with an explicitly configured optional environment, no implicit template downloads or license acceptance. CAPIVARA core CI must still pass when that optional backend is unavailable; a separate integration job must actually exercise it. Publishing should depend on successful code/documentation checks rather than treating a rendered homepage as release approval.

## 9. Proposed deletion, archive and deprecation list

These are proposals, not deletions performed during this audit.

| Exact target / group | Proposed action | Condition / preservation |
| --- | --- | --- |
| `Rplots.pdf` at repository root | Delete accidental generated plot | Tracked, excluded by .Rbuildignore; no use found in inspected source/docs. Final reference check and recoverable deletion after approval |
| Former Manim trial and its five planning documents | No further action | Already removed under the earlier explicit request; recoverable in `/Users/rd23aag/.Trash/capivara-manim-video-trial-20260910/`. Do not repeat deletion |
| `research/capivaraKinematics_prototype/` | Freeze and prominently mark historical/non-installable; optionally archive outside the active checkout | Preserve origin SHA and divergence map; retire second-package installation instructions |
| `research/capivara2/` demonstrations | Keep outside package build; add a concise provenance/status index | They contain useful science experiments, not a second public package |
| `research/capivara2/demo_capivara2_clean_pngs.R` | Repair its moved source path or archive the convenience launcher | It still sources `scripts/make_capivara2_paper_figures.R`, while the script is now in `research/capivara2/` |
| `research/capivara2/export_capivara2_workflow_assets.R` | Retain as figure-generation research code | This is not the Manim trial; synthetic/idealised assets must remain labelled |
| `inst/tutorials/run_centa_magnum_*` | Move observation-specific, path-dependent variants to `research/` after review; retain one portable MUSE/prepared-array example | Do not delete the evidence of generic input use; consolidate only after a validated adapter example exists |
| `inst/tutorials/run_full_manga_science_workflow.R` | Revise or move research orchestration out of the installed tutorial set | It assumes local backend/template checkouts and is not the ordinary optional-backend user path |
| `inst/tutorials/capivara_full_workflow.R` | Keep only with explicit teaching/simplified-fit labelling | Its demonstration fit must not represent validated stellar-population inference |
| `inst/extdata/kinematics/native_kinematics_workflow.R` and `native_bisymmetric_workflow.R` | Eventually move runtime bodies into ordinary internal R functions | They are currently executed by public wrappers: **do not delete or archive them now** |
| `R/kinematics_run_one_galaxy.R` and `R/kinematics_run_batch.R` | Deprecate duplicate internal orchestration after fixture parity | Retain shared readers and fitting primitives; preserve historical model comparisons |
| Internal `classify_structures()` alias | Deprecate in favour of neutral `catalogue_structures()` | Search research callers and provide a migration note |
| `R/SpecMean.R`, `R/cube_cluster_with_snr.R` and older plotting helpers | Candidate internal consolidation, not immediate deletion | Resolve callers and prove numerical/plot parity first |
| capivaraPPXF `R/fit.R` scaffold and duplicated preparation blocks | Replace/deprecate the stopping entry point; consolidate preparation | Implement and test the common fitting path first |
| `inst/doc/` and generated `docs/` pages | Regenerate from authoritative sources | Preserve current user edits; do not hand-patch generated HTML as the durable fix |
| Repeated figure assets / ignored outputs and build products | Inventory hashes and references before optional cleanup | No blanket deletion of research data, cached pilot results, tarballs or unrelated user files |

Research and docs are already excluded from the R source build. Repository organisation should improve discoverability, not merely reduce package size. Retaining a small frozen prototype is acceptable if users cannot confuse it with the maintained installation.

## 10. Implementation sequence ranked by scientific risk

| Order | Work package | Exit criterion before the next dependent stage |
| --- | --- | --- |
| 1 | Preserve source/pilot provenance; mark unreleased development version; add failing regression tests for C01–C11 | Every confirmed defect has a small reproducer; current and installed implementations are distinguishable |
| 2 | Correct flux/variance validity, wavelength/feature separation, eligibility and FITS channel/flag indexing | Summaries, maps and fits refer to the same measured samples and physical quantities |
| 3 | Resolve the distance contract and coordinate conventions; implement versioned conversions and compatibility | Exact/sparse comparisons isolate approximation; image/sky/disc angles and overlays pass analytic tests |
| 4 | Repair native line units/masks/uncertainties and remove hidden configuration; validate MaNGA and MUSE adapters | Same physical prepared input yields adapter-independent products; all effective settings are reported |
| 5 | Correct pPXF rebin/preparation, variance use, template coverage, normalisation and QC | Synthetic fits and both pilot rebin regressions pass; bound-hit solutions are not silently science-approved |
| 6 | Standardise structure candidates/accepted geometry and add segment overlap; calibrate bar negatives | Bar association has explicit coverage and mixed states; detector performance is assessed beyond selected positives |
| 7 | Calibrate axisymmetric/bisymmetric fitting, geometry assumptions, imputation treatment and regularisation | Controlled zero-bar and injected-bar tests quantify recovery and false positives; diagnostic names match their reference model |
| 8 | Consolidate script orchestration behind existing APIs and optional fitting facade; archive duplicate paths | Installed-package integration reproduces the approved reference fixtures without external global state |
| 9 | Complete metadata, NEWS, reference and scientific articles; render and visually inspect the local site | Public claims and every example match the tested implementation; desktop/mobile and dark/light QA documented |
| 10 | Release 0.4 or 2.0 according to completed gates, not naming preference | Tagged, archived, reproducible software with scoped validation and explicit remaining experimental modules |

Independent work packages can overlap once their contracts are agreed, but do not regenerate interpreted science products before their upstream coordinate, measurement and uncertainty fixes. In particular, rerun the five-galaxy pilot after the mapping/geometry and pPXF corrections; do not use its current outputs to tune bar thresholds.

**Review stop:** approve or amend this architecture, distance/legacy policy, scientific status boundaries and risk sequence before implementation. The audit recommends extending the package already present, not starting a new one.
