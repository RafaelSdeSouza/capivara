# CAPIVARA capability and API map

Companion to [CAPIVARA_2_COHERENCE_AUDIT.md](CAPIVARA_2_COHERENCE_AUDIT.md), 10 September 2026. This maps the inspected implementation, not a hypothetical final API.

## Classification and provenance keys

The requested classes are abbreviated only to keep the tables readable:

- **PUB:** production public API — established intended surface, not a scientific-readiness certification.
- **EXP:** experimental public API — exported, but contracts/calibration remain incomplete.
- **INT:** internal implementation.
- **RES:** research prototype.
- **DUP:** duplicated/obsolete implementation.
- **EXT:** external backend/dependency.
- **PLAN:** planned/not yet implemented.

Version keys: **C** = CAPIVARA source `0.3.0`, `156c60d02e4c014e9895ea7a68f2071ce7400796`; **I** = installed CAPIVARA `0.3.0`, `467b6ab2f99fddac8ea9cbf1e3daf3eb1f68ada2`; **S** = SpectroPath `0.1.0`, `e1695030c50c09adf5eebdde56336c2342493c25`; **P** = capivaraPPXF `0.0.0.9000`, `3195657b00cf83c2844e74983e7e66444986b533`. CAPIVARA result objects do not consistently store these SHAs. C and I are not interchangeable; the main audit explains workflow-file resolution.

Three linked tables use the same capability IDs. Read their rows together to obtain all requested fields without an eleven-column table.

## A. Implementation, public entry point, and result

| ID / capability | Class | Implementation files | Exported entry point | Current result/object |
| --- | --- | --- | --- | --- |
| A01 Exact spectral segmentation | PUB | `R/segment.R`, `internal_segment_core.R`, `torch_dist.R`, `cube_to_matrix.R`, `median_scale.R` | `segment()` | List: labelled `cluster_map`, representatives, S/N, input/header metadata and `original_cube` |
| A02 Sparse/large-cube segmentation | PUB | `R/segment_large.R`, `src/sparse_ward_cut.cpp`, registered Rcpp bridge | `segment_large()` | Segmentation list with backend diagnostics, graph settings and optional feature-window metadata |
| A03 S/N component selection and memory planning | PUB | `R/choose_ncomp_by_snr.R`, `estimate_segment_memory.R` | `choose_ncomp_by_snr()`, `estimate_segment_memory()` | Selection result/grid; one-row memory-estimate data frame |
| A04 Spatial support/starlets | PUB; helper INT | `R/starlet_layer.R`, `adaptive_support.R`; `detect_support()` in `structural_awareness.R` | `build_starlet_mask()`, `build_adaptive_support()` | Spatial mask/support products and method diagnostics; internal starlet decomposition |
| A05 Regional spectra | PUB | `R/summarize_cluster_spectra.R` | `summarize_cluster_spectra()` | List: wavelength, cluster IDs, spaxel/finite counts, median/mean/sum spectra and optional variance/weighted products |
| A06 Representative and flux-preserving reconstruction | PUB | `R/reconstruct_cluster_cube.R` | `reconstruct_cluster_cube()`, `reconstruct_flux_preserving_cube()` | List: model cube, templates, summary, flux check and optional residual cube |
| A07 Spectral/spatial visualisation | PUB; alternate helper INT | `R/plot_cluster.R`, `plot_cluster_spectra.R`, `capivara_plot_spectra.R` | `plot_cluster()`, `plot_cluster_spectra()` | ggplot/composed plot products; no independent physical inference |
| A08 Neutral structural awareness | INT | `R/structural_awareness.R`: `score_structures`, `threshold_structures`, `catalogue_structures` | None | Score/mask lists and component catalogue; `classify_structures` is a catalogue alias |
| A09 Structure-aware segmentation | INT | `R/structural_awareness.R`: augmentation and `segment_structures`, sparse backend | None | Segmentation list, currently with augmented feature channels in `original_cube` |
| A10 Bar and ring detection/geometry | INT | `R/structural_awareness.R`: `detect_bar`, `detect_ring`, ellipse/profile/Ferrer helpers | None | Geometry, masks, scores, profile/body diagnostics and heuristic acceptance fields |
| A11 Emission-line catalogue and spectral-feature segmentation | EXP | `R/segment_emission_lines.R` | `emission_lines()`, `segment_emission_lines()` | Line catalogue; segmentation plus line/window/feature and sampling metadata |
| A12 Native line extraction and conventional maps | INT, accessed through EXP workflow | `inst/extdata/kinematics/native_kinematics_workflow.R`, `R/kinematics_clean_galaxy_support.R` | Through `segment_kinematics()` / `run_kinematic_analysis()`, not an independent extraction API | Native list/maps/table: line-sum proxy, velocity, observed sigma, asymmetry, moment proxies, support and imputation flags |
| A13 SpectroPath mathematical representation | EXT; integration INT/EXP | SpectroPath `R/path.R`, `chen.R`, `logsig.R`, `streamsig.R`, `path-features.R`, `classical-features.R`, `algebra-helpers.R`; native workflow integration | CAPIVARA exposes `segmentation_mode="path_signature"` through workflow APIs; dependency has its own mathematical API | Ordered-profile path/log-signature and feature products; native feature cube and region map |
| A14 Kinematic segmentation | EXP | `R/kinematics_run_manga_bar_model.R`, native workflow | `segment_kinematics()` | `capivara_kinematic_segmentation` list: selected segmentation, native products, paths, plot, redshift metadata |
| A15 Axisymmetric modelling | EXP; primitives INT | `R/kinematics_disc_model.R`, `kinematics_geometry.R`, native bisymmetric workflow; older arctan primitive in `kinematics_fit_disc_model.R` | `run_kinematic_analysis(model="axisymmetric")` | `capivara_kinematic_result`, also assigned legacy `capivara_manga_bar_result` class; model/geometry/spaxel diagnostics and paths |
| A16 Bisymmetric/bar modelling | EXP; primitives INT | `R/kinematics_disc_model.R`, `kinematics_geometry.R`, `kinematics_run_manga_bar_model.R`, `inst/extdata/kinematics/native_bisymmetric_workflow.R` | `run_kinematic_analysis(model="bisymmetric_bar")`, `run_manga_bar_model()` | Same wrapper result; radial Vt/V2r/V2t profiles, projected components, fitted/residual velocity and geometry |
| A17 Model catalogue, residual diagnostics and panels | EXP; diagnostics INT | `R/kinematics_modules.R`, `kinematics_residual_diagnostics.R`, `kinematics_panels.R`, `kinematics_plot_capivara_kinematics.R` | `kinematic_models()`, `kinematic_panels()`; registered print/plot methods | Fixed two-row model table; diagnostic list; named ggplot panels / composed figures |
| A18 Regional pPXF preparation and mapping | EXT | capivaraPPXF `R/spectra.R`, `io.R`, `fit.R` | Backend `extract_capivara_spectra()`, `as_ppxf_input()`, `write_ppxf_input()`, `map_ppxf_results()`; no CAPIVARA facade | `capivara_ppxf_input`; exported tables and map lists |
| A19 Stellar-population/gas fitting | EXT, experimental implementation | capivaraPPXF `R/population.R`, `emission.R`, `diagnostics.R`, `quality.R` and environment/template helpers | Backend `fit_ppxf_population()`, `fit_ppxf_emission_lines()`, `fit_ppxf_population_diagnostics()`, `emission_quality_flags()` | Per-bin fit/diagnostic tables with population or gas quantities; numerical status and partial QC |
| A20 Generic fitting entry point | PLAN | capivaraPPXF `R/fit.R` contains a stopping scaffold; CAPIVARA facade absent | Backend `fit_ppxf()` is exported but always stops | No fitted result from the generic scaffold |
| A21 MaNGA ingestion and metadata | INT, used by EXP workflows | `R/kinematics_manga_metadata.R`, `kinematics_read_manga_maps.R`, FITS helpers; native cube reader; tutorial/research readers | No stand-alone exported cube adapter; workflow `cube_path` inputs | FITS-like cube/native data, DRP redshift/IDs, or DAP map list depending on reader |
| A22 MUSE and generic IFU ingestion | Generic arrays implemented; prepared MUSE RES; adapter PLAN | `R/internal_segment_core.R` input conversion; `inst/tutorials/run_centa_magnum_*` prepared-array demonstrations | Arrays accepted by spectral/emission APIs; no dedicated exported MUSE reader | FITS-like `imDat` list or raw array; no validated common observation object |
| A23 Segment–structure association | PLAN | Research masks/catalogues are available, but no general association API | None | Proposed ID-keyed area/flux overlap table with mixed/unclassified states |
| A24 Legacy orchestration / duplicate package | INT/DUP, historical copy RES | `R/kinematics_run_one_galaxy.R`, `kinematics_run_batch.R`, `kinematics_read_capivara_output.R`; `research/capivaraKinematics_prototype/` | None in current CAPIVARA for the legacy runners | File/config-driven reports, maps, model results; second package surface survives in prototype |
| A25 CAPIVARA 2 science/figure demonstrations | RES | `research/capivara2/` seven R scripts and assets; related research MaNGA/M83 experiments | None | Demonstration segmentations, structure/fit maps and figure assets, not a separate public analysis package |

## B. Input, survey, and coordinate assumptions

| ID | Input assumptions | Survey assumptions | Coordinates / units / semantic caveat |
| --- | --- | --- | --- |
| A01 | 3-D array or FITSio-like list; valid positive summed signal; median-centred, zero-filled features | Algorithm not intrinsically survey-specific | Array `[row,col,wave]`; flattening in R order; distance is currently L1→Ward.D2; physical wavelength often relies on `axDat` |
| A02 | Same basic cube; additional valid modes, finite/MAD thresholds, graph and optional spatial features | Generic prepared arrays | Euclidean Ward increments; native pixel coordinates if spatial weight used; graph clusters need not be spatially connected |
| A03 | Spectral matrices and variance/target S/N, or valid-pixel count/cube for memory | Generic; flux-as-counts noise fallback is not generically physical | Wavelength-window semantics depend on available metadata; memory is an estimate, not allocation guarantee |
| A04 | Image/cube; chosen multiscale/band settings; finite data and sufficient image size | Generic, but tuning is in pixel/scales rather than a uniform physical resolution | Mask aligned to native rows/columns; selected footprint, not probability or bar evidence |
| A05 | Label map plus retained real cube; optional variance with matching samples | Generic in principle | Sum vs mean vs inverse-variance mean differ; no-contribution and incomplete variance defects C01; wavelength may fall back to channel index |
| A06 | Regional result and template choice; mask/fill policy | Generic | Mean with preserved measurement mask can conserve sums on observed support; median reconstruction is not the same guarantee; no new independent data |
| A07 | Segmentation result; spectra plot expects `axDat` | Examples mainly MaNGA | Row/column plotting differs between helpers; spectral axis labelled Å without a complete unit contract |
| A08 | Whitelight/image/cube and support; ridge/concentration/symmetry choices | Generic image mathematics, tuned mainly on MaNGA demonstrations | Native pixel sizes; component PA from +x differs from bar PA from +y; score is not physical morphology probability |
| A09 | Scores and observed cube; repeated structure-feature channels appended | Generic mathematics | Feature axis is currently mistaken for wavelength axis downstream; spectral flux products unsafe until C07 |
| A10 | Scores/image, central support and threshold/profile options | Tested on selected barred MaNGA examples; not calibrated across IFUs | b/a and ellipticity recorded; PA conventions differ; candidate vs accepted mask differs for ring/bar; disc-plane angle not supplied by image PA alone |
| A11 | Explicit redshift, line selection/windows; cube and physical wavelength; optional sampling | Catalogue and windows assume optical-line coverage | Air/vacuum not explicit; missing profiles zero-filled; sampled/non-sampled scaling differs |
| A12 | Cube filename readable by native script, redshift, selected optical line and support controls | MaNGA-oriented cube/HDU/header conventions | Native row/column maps; centred velocities; channel-sum flux proxy, uncorrected width, not GH h3/h4; imputation must remain distinct |
| A13 | Ordered profile coordinates and chosen transformations/feature names | Mathematical dependency generic; native adapter MaNGA-oriented | Profile order and normalisation matter; features are not velocities or a calibrated bar classifier |
| A14 | Same as A12 plus segmentation mode and graph settings | Public name broader than current file adapter | Selected label map is returned; result not guaranteed to share spectral-region IDs; maps inherit C02/C08/C12 risks |
| A15 | Native velocity map, support, geometry and regularisation controls | Fitting mathematics generic; public runner MaNGA-oriented | NIRVANA/legacy deprojection; fallback inclination possible; image PA not sky PA; spaxel-wise regularised fit |
| A16 | Same plus bar-angle source or estimator, second-order controls, optional wedge | Generic equations; wrapper and provenance still MaNGA-heavy | Disc-plane bar angle separate from image bar PA; sign mismatch in label fallback; diagnostic wedge not independently measured bar footprint |
| A17 | Known result structure or selected model string | Generic in principle; panel metadata can still say MaNGA/Hα | Residual reference depends on fitted model; `f_bar` is squared-residual-power fraction; point/image display mismatch C06 |
| A18 | Correct regional spectra and wavelength/IDs; optional variance; fitting metadata | Backend not intrinsically MaNGA-specific | Sum/mean/ivar mean differ; variance is carried but not necessarily used by fitters; mapping requires unique validated IDs |
| A19 | Python pPXF environment, SPS templates, redshift, optical fit range, LSF/noise assumptions | Nominal fixed FWHM/default ranges are not valid for every IFU | Separate template/galaxy grids; coefficients normalised; bound-hit QC missing; lightweighted log age and M/L must not be relabelled linear age/mass |
| A20 | Intended pPXF input/templates | Planned survey-neutral facade | No present scientific output to interpret |
| A21 | MaNGA filename/headers/catalogue; DRP cube vs DAP MAPS depends on path | Explicitly MaNGA; MEGACUBE requires verified original-array adapter | Preserve WAVE, IVAR and mask definitions; DAP results are external measurements, not CAPIVARA's fresh segmentation |
| A22 | Correctly prepared arrays currently; future MUSE DATA/STAT and wavelength/header adapter | Demonstrations use Cen A/MAGNUM local paths; no general MUSE reader certification | STAT is variance, not IVAR; wavelength medium/sampling/WCS/LSF must be resolved in adapter |
| A23 | Shared grid or explicit WCS resampling, accepted structure mask and segment IDs | Survey-neutral reusable operation; environment science belongs outside package | Area fraction vs band-specific flux fraction; mixed regions and missing coverage explicit |
| A24 | Config files, label/segment tables, MaNGA DAP or saved native products | Survey-specific historical path | Alternative result schemas, outdated arctan report text, duplicated geometry implementations |
| A25 | Local cube/cache/backend paths and research assumptions | Mostly MaNGA examples, some other targets | Smoothed physical-result displays and idealised illustrations are not new independent measurements or validation truth |

## C. Documentation, tests, examples, and provenance

CAPIVARA has no maintained test suite at C. “Probe” below refers only to the controlled audit checks, not an existing testthat regression. Existing package-check success is not repeated as a test for each scientific capability.

| ID | Documentation page/state | Existing example | Tests / validation at inspection | Provenance / version state |
| --- | --- | --- | --- | --- |
| A01 | `reference/segment`; Getting started; Backend choice; original paper | Synthetic vignette and research MaNGA comparisons | No package tests; audit metric probe; paper describes Euclidean method, not every current path | C/I both display 0.3.0; distance/scaler need explicit history |
| A02 | `reference/segment_large`; Backend choice | Backend comparison figure and research benchmark | No package tests; full-graph audit probe; empirical demos not approximation guarantees | C; sparse backend described in 0.3 NEWS; previous 0.2 implementation differs |
| A03 | Both function references; Backend choice | Memory estimate and target-S/N examples | No maintained tests; validity/variance propagation needs analytic cases | C; existing output does not identify a complete uncertainty model |
| A04 | Both references; Support masks | Synthetic vignette, research starlet/MAGNUM comparisons | No tests; small-image padding failure reproduced | C; masks have settings, not calibrated completeness/purity |
| A05 | Function reference; Regional spectra | Sum/weighted summaries in vignette; backend preparation | No tests; missing-data failure reproduced | C; no full source/input provenance embedded |
| A06 | Both references; Regional spectra; FITS and DS9 | Reconstructed cube and flux-check examples | Internal flux check, no maintained round-trip or uncertainty tests | C; input summary and masks determine interpretation |
| A07 | Both public references; Getting started | Synthetic map/spectra and real mosaic | No coordinate snapshot/axis tests | C; older alternative plotting helper remains internal |
| A08 | R roxygen/internal comments; organisation document; no public article | `demo_capivara2_structural_awareness.R` | No maintained tests; PA discrepancy reproduced | C, deliberately hidden; research demonstrates rather than validates |
| A09 | Internal comments; no reference page | `demo_capivara2_structure_support_segmentation.R` | Augmented-spectrum failure reproduced | C; not safe as fitting input without correction |
| A10 | Internal detector documentation; no public reference | `demo_capivara2_bar_detection.R`; iFUN five-positive pilot | Positive/perturbation evidence only; no complete raw-input negative-control calibration | C detector body matches I; score-map pilot has narrower scope than full pipeline |
| A11 | R roxygen exists; both exported Rd/reference pages missing | `inst/tutorials/run_emission_line_segmentation.R`, MAGNUM line-mode comparisons | No tests; missing-spaxel and metadata probe; branch semantics source-audited | C, August post-0.3 addition without new package version |
| A12 | Kinematic Analysis article describes workflow, not full estimator semantics | `inst/tutorials/run_kinematic_analysis.R` | No tests; exported-label defect reproduced; line-profile recovery absent | C runtime script; hard-coded Hα/MaNGA provenance needs correction |
| A13 | Brief Kinematic Analysis mention; dependency README/API | Explicit path_signature workflow mode | S has a small feature/wrapper test set; no CAPIVARA observational integration/calibration suite | C+S; no pinned installed S SHA available |
| A14 | Function reference; Kinematic Analysis | Native workflow tutorial | No tests; no generic adapter fixture | C source signature differs from older installed workflow surface |
| A15 | `run_kinematic_analysis` reference; Kinematic Analysis | Axisymmetric tutorial | No injection/recovery suite; placeholder geometry possible | C; provisional model inference, not a published validation of this implementation |
| A16 | `run_manga_bar_model` reference; Bar Model article | `inst/tutorials/run_bisymmetric_bar_model.R`; selected pilot | Index/sign probes fail; no calibrated model-recovery/false-positive study | C/I wrapper differs; prototype also diverges; record effective estimator/angle source |
| A17 | Both public references, Kinematic Analysis and Bar Model | Named panel extraction in vignettes | No catalogue/dispatch/diagnostic invariance tests; rectangular overlay probe fails | C; fixed catalogue, no genuine plug-in discovery |
| A18 | Backend roxygen/README and core Regional spectra | `inst/tutorials/run_full_manga_science_workflow.R` | Backend setup/I/O/preparation tests; core missing-variance defect propagates | C+P; external Python not required for preparation alone |
| A19 | Backend references; CAPIVARA research examples | Saved pilot population/emission products | Rebin regression reproduced on two wavelength arrays; no fresh fitting or recovery validation in audit | P; Python/template provenance and normalisation need recording |
| A20 | Backend scaffold documentation correctly says future | Export/import route suggested by scaffold | Always stops; not an implemented fit | P; CAPIVARA facade remains PLAN |
| A21 | Kinematic Analysis; tutorials; internal metadata documentation | MaNGA LOGCUBE and historical DAP workflows | No maintained MaNGA fixture; short-header read succeeds, real warning unresolved | C; catalogue release and extension sources must be explicit |
| A22 | Broad generic wording; prepared-array tutorials, no dedicated MUSE guide | `run_centa_magnum_whole_vs_emission_windows.R` and related scripts | Demonstrations only; no DATA/STAT/WCS adapter parity test | C research use, not a validated survey adapter release |
| A23 | Existing iFUN integration design proposes associations | Ad hoc masks and scientific intent only | No general implementation/test | PLAN; do not add to current export list or advertise as shipped |
| A24 | Prototype README/man/configs; repository-layout and organisation docs disagree with second install narrative | Prototype examples and batch configs | 59 shared functions compared: 52 identical bodies, seven divergent | Historical copy and C; freeze rather than maintain parallel science |
| A25 | Research script comments/organisation docs | Seven CAPIVARA 2 scripts, cached pilot/figure products | Demonstrations, not maintained package tests; one moved source path stale | RES; provenance varies by local file/cache; excluded from package build |

## D. Complete R-source coverage

Each row gives the scientific or technical role of inspected source files. This is a coverage index, not a recommendation to preserve the present file layout indefinitely.

| Files under `R/` | Role / linked capabilities |
| --- | --- |
| `capivara-package.R` | Package-level scope/documentation |
| `globals.R` | Static-check declarations; not a scientific estimator |
| `RcppExports.R` | Registered compiled sparse-Ward bridge (A02); paired with `src/RcppExports.cpp` |
| `segment.R` | Exact user entry point (A01) |
| `internal_segment_core.R` | Input conversion, wavelength and common segmentation helpers (A01–A03, A22) |
| `cube_to_matrix.R` | Array-to-spectrum matrix ordering (A01–A06) |
| `median_scale.R` | Median centring, not generic unit normalisation (A01–A02) |
| `torch_dist.R` | Optional torch/base distance implementation (A01) |
| `segment_large.R` | Sparse graph, row/column scaling, validity and output assembly (A02) |
| `choose_ncomp_by_snr.R` | Target-S/N selection and variance aggregation (A03) |
| `estimate_segment_memory.R` | Memory estimates (A03) |
| `starlet_layer.R` | Multiscale filtering/support (A04) |
| `adaptive_support.R` | Alternative data-driven footprint (A04) |
| `summarize_cluster_spectra.R` | Regional flux/variance summaries (A05) |
| `reconstruct_cluster_cube.R` | Representative/mean-preserving cube products (A06) |
| `plot_cluster.R`, `plot_cluster_spectra.R`, `capivara_plot_spectra.R` | Public and alternative spectral/map displays (A07) |
| `SpecMean.R`, `cube_cluster_with_snr.R` | Older internal representative/SNR helpers (A01, A03, A07); consolidate only after caller checks |
| `structural_awareness.R` | Support wrapper, neutral structure scores/catalogue, bar/ring geometry and feature-augmented segmentation (A08–A10) |
| `segment_emission_lines.R` | Optical line catalogue, window/profile features and sampled/full segmentation (A11) |
| `kinematics_capivara_kinematics_utils.R` | Dependency/file/FITS/map conversion utilities (A12–A17, A21) |
| `kinematics_dependencies.R` | Native workflow dependency declaration |
| `kinematics_clean_galaxy_support.R` | Connected support cleaning and mask morphology (A04, A12) |
| `kinematics_manga_metadata.R` | MaNGA IDs/redshift/table/header resolution (A21) |
| `kinematics_read_manga_maps.R` | MaNGA DAP map reader (A21, A24), distinct from raw cube analysis |
| `kinematics_read_capivara_output.R` | Legacy label/table result reader (A24) |
| `kinematics_geometry.R` | Disc geometry and NIRVANA/legacy deprojection (A15–A16) |
| `kinematics_disc_model.R` | Arctan prediction helper, piecewise linear solves, bar-axis estimation and axisymmetric/bisymmetric models (A15–A16) |
| `kinematics_fit_disc_model.R` | Older nonlinear arctan fitting primitive (A15, A24) |
| `kinematics_residual_diagnostics.R` | Model-residual diagnostics and designated bar-region comparison (A17) |
| `kinematics_modules.R` | Fixed two-model catalogue/selector (A17) |
| `kinematics_run_manga_bar_model.R` | Public native segmentation/model wrappers, environment/file orchestration and S3 methods (A14–A17) |
| `kinematics_run_one_galaxy.R` | Legacy config runner and shared map-to-spaxel assembly (A24, A15–A16) |
| `kinematics_run_batch.R` | Legacy batch runner (A24) |
| `kinematics_plot_capivara_kinematics.R` | Image/point transforms, white-light/model/component plots and bar-axis overlay (A17) |
| `kinematics_panels.R` | Named panel extraction (A17) |

## E. What should become visible, and when

1. Immediately after the documentation/version correction: expose existing line-segmentation exports, the actual model assumptions, and the separate role of SpectroPath. Do not imply that documentation creates new functionality.
2. After geometry and result-schema regressions pass: expose neutral structure scores/catalogues, followed by explicitly experimental bar/ring geometry.
3. After implementation and analytic overlap tests: expose segment–structure association for bar/non-bar/mixed regional comparisons.
4. After adapter fixtures and ordinary-function extraction: document native kinematics as general IFU analysis, with MaNGA and MUSE as adapters rather than scientific defaults.
5. After pPXF preparation/QC validation: provide a single optional fitting entry point in CAPIVARA while retaining the specialised external backend.

The full implementation and release order is in section 10 of the [main audit](CAPIVARA_2_COHERENCE_AUDIT.md).
