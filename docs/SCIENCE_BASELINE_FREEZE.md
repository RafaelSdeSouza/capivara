# CAPIVARA scientific baseline freeze

Validated on 2026-09-11 (UTC). Recommendation: `NOT_READY_FOR_SANDRA_PILOT`.

The corrected software passes its regression tests and the five-object
reproducibility checks. All five spectral partitions are unchanged, and their
automatically detected bar masks now rotate exactly with the input cubes.
The joint bar, kinematic and stellar-population science analysis is not yet
validated: the population fits retain substantial residual mismatch, the bar
classifier has no calibrated false-positive test here, and the deprojected
kinematics still use an explicitly labelled inclination placeholder.

Automatic bar-mask construction and projected major-axis measurement require
no manually supplied inclination. Inclination enters the conversion to disc-plane
angles and velocities, not the image-plane mask.

## Frozen sources and installation

| Package | Version | Validated code commit |
| --- | --- | --- |
| capivara | 0.4.0.9000 | `2cdce2f54526ecc9b04eb3aa7f2a54b51f1e7100` |
| capivaraPPXF | 0.0.1.9000 | `9bfdea212e1704e80ca74dd715801df893e8ca07` |

Both use the local branch `fix/capivara-science-baseline`. The authoritative
sources are the `worktrees/capivara-science-baseline` and
`worktrees/capivaraPPXF-science-baseline` directories under the iFUN
`Capivara_Eat_Manga` project. A subsequent CAPIVARA commit records this report
only; the table pins the tested implementation, not a self-referential report SHA.

The original CAPIVARA checkout remains at `156c60d02e4c014e9895ea7a68f2071ce7400796`,
with its pre-existing website and documentation edits preserved. The old
installed CAPIVARA reported 0.3.0 and remote SHA
`467b6ab2f99fddac8ea9cbf1e3daf3eb1f68ada2`; important wrapper bodies differed
from the source that also reported 0.3.0. The previous capivaraPPXF source was
`3195657b00cf83c2844e74983e7e66444986b533`, version 0.0.0.9000.

Both corrected packages were installed from their authoritative sources into
the isolated validation library and the normal R library,
`/Users/rd23aag/Library/R/arm64/4.5/library`. Old normal-library installations
were copied to `results/science_baseline_freeze/previous_library`.
Restart any R session that already loaded an older namespace.

`tools/check_source_install.R` verifies version, function bodies and formals,
namespace exports, and bundled R/Python scripts. It passes for CAPIVARA's 187
top-level functions, 18 exports and 12 bundled scripts, and capivaraPPXF's 17
functions, 13 exports and one Python script. It reports whether the source
worktree is modified. Compiled code is rebuilt and exercised by the tests;
the checker does not claim binary reproducibility across compilers.

capivaraPPXF pins its CAPIVARA dependency to the code commit above. Neither
branch has been pushed; installing that remote SHA or running the remote CI
requires a later, explicitly authorised publication of the commits.

## Defects fixed

1. **Inconsistent coordinate meanings.** Component PCA angles used a different
   zero axis from the weighted bar ellipse. Component and ring outputs now use
   the canonical image convention; disc-plane bar angles have separate names.
   The white-light estimator now reports minor/major, rather than its inverse.
   A rectangular counterclockwise display transform used the wrong dimension
   for a transformed coordinate; point, image and mask transforms now agree.

2. **Orientation-dependent bar centre.** The detector clipped the brightest
   two percent of the collapsed light before choosing its first maximum.
   The resulting plateau made the centre depend on array traversal order.
   In 11004-12701, that chose row 40, column 32 instead of the raw light peak
   at row 37, column 37. The centre is now determined before clipping; genuinely
   tied peaks are averaged. Before this fix, 90-degree cube rotations gave
   additional PA errors of 15.73, 65.39, 55.13, 28.50 and 55.94 degrees in
   the five pilots. The acceptance thresholds and profile model were not tuned.

3. **Unpaired pPXF log grids.** A second `log_rebin()` call with the velocity
   scale from the first call could return one fewer pixel. All region spectra
   are now rebinned together once and retained with that call's returned
   log-wavelength grid. Galaxy, noise and wavelength lengths are checked.
   Templates retain their own wavelength grid and sufficient coverage; they
   are not arbitrarily cropped to hide the mismatch. This follows the
   [pPXF wavelength and template contract](https://pypi.org/project/ppxf/).

4. **Discarded available variance and misleading fit success.** Spectral sums,
   inverse-variance means and their errors use compatible finite contributors.
   Missing contributions no longer become measured zeros or artificially
   precise sums. With rebinning operator W, diagonal input variances propagate
   as the diagonal of W diag(V) W-transpose. Missing spectral support is excluded
   from `goodpixels`. Empirical noise is allowed only as a labelled preview
   when variance is absent. Numerical completion, optimizer convergence,
   kinematic bounds and necessary scientific-eligibility checks are separate.
   Undefined age/metallicity summaries cannot pass population QC.

5. **Single-region interface failure.** Reticulate can turn a length-one
   region-ID vector into a scalar. The backend now retains a region dimension
   and validates IDs/counts in all three active fitters. Diagnostic spectra
   expose `pixel_used_in_fit`; missing observed flux, residuals and noise are
   not plotted as numerical zero-fill measurements. Emission quality flags no
   longer infer BPT science validity from completion and positive template
   coefficients alone; explicit line-flux calibration is also required.

6. **Unsafe kinematic defaults and geometry transport.** Scientific mode
   requires an explicit inclination before cube I/O. Placeholder inclination,
   imputed velocities and geometric wedge masks are confined to explicit
   preview use where applicable. Supplied vetted geometry is converted once
   and its native-grid mask is preserved. Conflicting bar angles are rejected.
   Native-map imputation flags now reach the correct row/column spaxels.
   Missing bar support yields undefined bar diagnostics with a QC reason,
   rather than a measured-looking zero. Completion/convergence flags recognise
   the actual piecewise fitter status strings.

7. **Contaminated workflow state and output labelling.** Explicit source
   workflow paths take precedence when requested; inherited numerical
   environment controls are cleared and restored even after errors.
   The requested kinematic segmentation mode and centre/ring controls are
   honoured. FITS output plane names no longer have an off-by-one mapping.

8. **Other demonstrated spectral/data defects.** Structural feature channels
   are no longer retained as physical wavelength samples. Emission-feature
   eligibility excludes absent spectra and uses consistent preprocessing in
   sampled and unsampled branches. Explicit wavelength vectors are retained
   and validated. Starlet reflection supports fields smaller than its padding
   width while preserving the old padding for ordinary fields.

9. **Reproducibility hygiene.** Runtime and installed tutorial developer-path
   fallbacks were removed. Owned FITS connections close on success and error;
   header-only reads avoid loading an entire cube for metadata. Gzip expansion
   is streamed instead of relying on an estimated uncompressed size.
   Missing exported documentation was generated. Both repositories now have
   ordinary R CMD check CI workflows; the backend CI installs the tested Python
   versions and does not download SPS templates.

## Coordinate and geometry contract

Native arrays are [row, column, wavelength], with one-based x = column and
y = row. `pa_image_deg` runs from +row toward +column, modulo 180 degrees:
vertical is 0 and horizontal is 90. It is clockwise in an x-right, y-up
drawing. A display rotation changes plotted endpoints, not the stored
measurement. `pa_sky_deg` is east of north and requires an explicit WCS
conversion, including the FITS-axis to array-axis mapping. It remains missing
for these five outputs; no north-up frame was assumed.

`phi_bar_disc_deg` is an in-plane angle relative to the disc major axis and
uses the same handedness as the selected deprojection. It is not an image PA.
Axis ratio means b/a, and ellipticity means 1-b/a. See
[the coordinate specification](SCIENCE_COORDINATES.md) and
[the pre-rename PA inventory](SCIENCE_PA_INVENTORY.txt).
Ambiguous historical input aliases are warned or rejected; historical saved
catalogues must not be reinterpreted by blind field renaming.

The pilots call the current internal `detect_bar()`, then explicitly supply
its automatic mask and image PA to `run_kinematic_analysis()`. They do not use
the separate legacy white-light angle fallback. Detector-generated geometry
is labelled `vetted=FALSE`; synthetic rotation consistency does not make it
independent morphological truth. The model preserves the supplied mask and
PA in all five runs and does not estimate a second unrelated bar angle.

## Five-pilot validation

Only 11004-12701, 11014-3704, 8602-12705, 8932-3701 and 9869-9102 were run.
Input SHA-256 values were verified against the manifest. The analysis reads
the original MaNGA FLUX, IVAR, MASK and WAVE HDUs inside each MEGACUBE container;
it uses no MEGACUBE-derived segmentation, stellar populations or fitted spectra.

Full-spectrum segmentation retains the prior settings: 18 sparse-Ward
regions, k = 40 with automatic expansion up to 100, robust column scaling,
spatial weight 0.15, the 4800–7400 Angstrom rest-frame feature window, and the
detector's galaxy-support mask. All five runs return 18 regions with exactly
the same assigned footprint and label-invariant partition as the pre-fix
baseline. No scientific change to the segmentation was needed for these data.

### Automatic bar geometry

| Object | Image PA (deg) | b/a | Profile radius (pixels) | Accepted-mask pixels | Old/new mask IoU |
| --- | ---: | ---: | ---: | ---: | ---: |
| 11004-12701 | 137.6610 | 0.5505 | 17.6711 | 607 | 0.9436 |
| 11014-3704 | 80.6599 | 0.7136 | 12.7544 | 95 | 0.4502 |
| 8602-12705 | 91.1943 | 0.6359 | 11.7187 | 378 | 0.7500 |
| 8932-3701 | 108.3717 | 0.5700 | 10.4398 | 230 | 0.6610 |
| 9869-9102 | 146.3353 | 0.8044 | 16.3812 | 147 | 0.3568 |

All five satisfy the detector's existing bar-like acceptance rule. The table's
PA, b/a and radius describe its candidate ellipse/profile, not a calibrated
physical bar length. For 11014-3704 and 9869-9102, acceptance uses the existing
bright-body fallback rather than the profile criterion; their body PAs are
80.2763 and 146.1232 degrees. Those body diagnostics remain separately named.

Every full-cube 90-degree rotation gives mask IoU = 1 and zero mask-area
change. The maximum additional axial PA error is 2.84e-14 degrees.
The corrected PAs differ from the original outputs by 12.06, 7.05, 12.65,
9.61 and 39.31 degrees, respectively. These are material bar-geometry changes
caused by the established centre defect, despite unchanged spectral regions.

Each `segment_bar_overlap.rds` records region size, overlapping bar-mask pixels
and fractional bar coverage. No arbitrary binary membership threshold is
imposed on partially overlapping regions.

### Population fitting

All fits use summed region spectra, native positive IVAR and MASK == 0, the
existing E-MILES template archive, a 4800–7400 Angstrom rest-frame fit window,
the unchanged approximate 2.76 Angstrom FWHM setting and multiplicative
polynomial degree 8. The conservative fitting mask does not alter segmentation.
The paired log grids contain 2647, 2664, 2683, 2667 and 2673 pixels, respectively.

| Object | Numerical / converged | Bound regions | Necessary fit-level eligible | Failed regions | Bound region IDs |
| --- | ---: | ---: | ---: | --- | --- |
| 11004-12701 | 18 / 18 | 7 | 11 | none | 1, 4, 7, 8, 12, 17, 18 |
| 11014-3704 | 17 / 17 | 6 | 11 | 1 | 3, 6, 10, 11, 13, 17 |
| 8602-12705 | 17 / 17 | 3 | 14 | 16 | 1, 2, 5 |
| 8932-3701 | 16 / 16 | 0 | 16 | 2, 4 | none |
| 9869-9102 | 18 / 18 | 0 | 18 | none | none |

There are 86 numerical successes among 90 regions, 16 bound-limited fits and
70 fits passing the necessary fit-level gates. Each failed region has zero
usable rebinned samples after the native data-quality mask; none fails from
the repaired one-pixel grid mismatch. Previously, the first two objects had
zero completed population fits; the other three each had 18, using empirical
noise without the current native-variance/mask treatment.

`parameter_on_bound` tests stellar and gas velocity/dispersion limits, not
whether individual nonnegative SSP weights equal zero. The
`scientifically_usable` column applies necessary fit-level gates; it does not
certify template or LSF adequacy, population recovery, or uncertainty coverage.
Median reduced chi-square values are 28.66, 29.43, 14.21, 31.12 and 21.75.
They use the new diagonal-variance noise and are not directly comparable to
the previous empirically scaled-noise fits.

After the single-region interface repair, all 90 population fits were rerun
from the same saved inputs. Every fitted value and QC column agrees with the
preceding corrected run to tolerance 1e-12. Additional real-data smoke tests
pass for emission bin 2 and diagnostic bin 8 of 11004-12701. The latter retains
its complete 2647-pixel grid while identifying its 333 used samples and
leaving invalid observed samples explicitly missing.

### Kinematics

All five bisymmetric preview fits complete numerically and converge in the
robust solver; none triggers the implemented amplitude-instability flag.
The automatic velocity-gradient disc PAs are 165.7237, 131.4692, 3.6685,
102.5297 and 114.1874 degrees in the image frame. They are separate from the
bar PAs in the table.

All use the explicitly labelled 60-degree inclination placeholder, and all
have `scientifically_usable=FALSE`. Their deprojected bar angles and velocities
are therefore preview outputs. Bar masks and image angles are transported
unchanged into the model; bar-support diagnostics use the supplied mask.

The pilot warning logs contain benign numerical coercion warnings from the
optional centre/systemic-velocity environment values represented as the string
`"NA"`. These become missing values and invoke the documented estimators;
they are not fit failures. They remain a minor interface-cleanup item.

## Tests and package checks

CAPIVARA: 21 test cases, 155 passing expectations, no failures or skips.
Coverage includes horizontal/vertical/oblique geometry, all display-transform
round trips on rectangular arrays, actual clipped-plateau rotation failure,
mask rotation, WCS parity, image-to-disc-to-image axis recovery, b/a semantics,
vetted-geometry transport, missing-mask QC, variance contributors, missing
emission spectra, feature/wavelength separation, small-field padding, FITS
connection ownership, plane labels, reversible environment controls and
unchanged exact/sparse partition fixtures.

capivaraPPXF: nine test cases, 42 passing expectations, no failures; the one
skip is the deliberately missing-pPXF installation case because Python pPXF
is installed. Tests reproduce both pilot redshift/length grid failures,
check overlap-operator variance propagation and missing pixels, enforce
template coverage, retain scalar region IDs, and distinguish convergence,
bounds, undefined populations and emission-flux calibration.

Both final source archives pass `R CMD build` and `R CMD check --no-manual`:
**0 errors, 0 warnings, 0 notes**. CAPIVARA's vignettes were built and rechecked.
The local environment is R 4.5.2 on aarch64 macOS 15.6.1; Python pPXF 9.4.5,
NumPy 1.26.4 and SciPy 1.13.1. The installed SpectroPath is 0.1.0
(source SHA `e1695030c50c09adf5eebdde56336c2342493c25`).

Repository index access was unavailable during checks, but required local
dependencies were present and the final checks completed normally. The new
Linux/macOS/Windows CI workflows have not been run remotely because nothing
has been pushed.

## Remaining science limitations and readiness decision

- The automatic masks are heuristic candidates. These five known-bar pilots
  and rotation tests do not measure false positives, completeness, or mask
  fidelity against independent morphology. The profile and bright-body
  acceptance branches retain their existing distinct geometry diagnostics.
- No vetted inclination is available. The current native velocity estimator,
  velocity-error propagation, PSF effects and bisymmetric parameter recovery
  have not been validated for scientific inference in this freeze.
- Spatial covariance between MaNGA spaxels and rebinning-induced off-diagonal
  covariance are not included in fitting. Summed IVAR therefore supplies only
  a diagonal approximation. A spatial correction cannot be justified from
  region size alone for disconnected spectral regions. See the
  [SDSS covariance and LSF guidance](https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/).
- The wavelength-dependent native LSF has not replaced the approximate scalar
  fitting resolution. Large residual chi-square and boundary hits require
  an LSF/noise/residual assessment and population-recovery tests before
  interpreting ages and metallicities. No posterior intervals or calibrated
  population uncertainties are produced here.
- Gas template coefficients, tied Balmer templates and fixed-ratio doublets
  are not certified individual physical line fluxes. Scientific BPT use is
  gated off unless an external calibration is explicitly declared.
- Exact L1/Ward.D2 segmentation and sparse Euclidean/Ward segmentation remain
  different scientific configurations. Both previous behaviours are locked
  by regression fixtures; this task does not make them equivalent.
- Real-data validation here is MaNGA-only. The array-based core remains
  survey-independent, but this does not validate every MUSE/FITS adapter.
  The generic `fit_ppxf()` facade remains a scaffold; the tested active
  population, emission and diagnostic entry points are listed above.

The corrected code can be pinned for further controlled validation. The joint
environment/bar/population science pilot remains
`NOT_READY_FOR_SANDRA_PILOT` until bar-mask acceptance, spectral-model/noise
adequacy and any required deprojection are constrained. No Sandra-46 run,
hierarchical segmentation, website/README narrative revision or push was made.

## Reproduction and retained evidence

All paths below are relative to the iFUN `Capivara_Eat_Manga` project.
Final products are under `results/science_baseline_freeze/corrected`; the
first validation pass is retained separately in its parent directory.

- `five_pilot_results.rds`: complete per-object validation rows.
- `pre_fix_reference.rds`: frozen old masks, partitions and fit tables.
- `run_provenance.rds`: input manifest, settings, package paths and SPS hash.
- `backend_refit_validation.rds`: unchanged results after final backend repair.
- `additional_backend_smoke.rds`: single-region emission/diagnostic outputs.
- `pilots/<id>/`: masks, geometry, rotation QC, region/bar overlaps, fitted
  input spectra, population tables, preview kinematics and diagnostic plots.
- Parent `build/capivara.Rcheck/00check.log` and
  `build/capivaraPPXF.Rcheck/00check.log`: final package checks.
- Parent `source_install_checks.txt` and `freeze_manifest.rds`: installed/source
  comparisons, code pins and artifact fingerprints.

Use an existing project/SPS archive and a fresh R process. The runner has a
fixed five-object allowlist and accepts an output directory separate from the
old pilot products:

```sh
Rscript worktrees/capivara-science-baseline/tools/check_source_install.R \
  worktrees/capivara-science-baseline

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
Rscript worktrees/capivara-science-baseline/tools/validate_five_pilots.R \
  /path/to/Capivara_Eat_Manga /path/to/new-validation-output \
  /path/to/validated-R-library /path/to/spectra_emiles_9.0.npz \
  /path/to/science_baseline_freeze/pre_fix_reference.rds
```

The final backend-only refit and extra single-region checks are retained as
`results/science_baseline_freeze/refit_final_backend.R` and
`check_backend_entrypoints.R`. No new templates or survey cubes were downloaded.

## Changed files

CAPIVARA implementation, tests and technical documentation, relative to
`156c60d`:

```text
.github/workflows/R-CMD-check.yaml
DESCRIPTION
NAMESPACE
R/choose_ncomp_by_snr.R
R/internal_segment_core.R
R/kinematics_capivara_kinematics_utils.R
R/kinematics_disc_model.R
R/kinematics_geometry.R
R/kinematics_manga_metadata.R
R/kinematics_plot_capivara_kinematics.R
R/kinematics_read_capivara_output.R
R/kinematics_residual_diagnostics.R
R/kinematics_run_batch.R
R/kinematics_run_manga_bar_model.R
R/kinematics_run_one_galaxy.R
R/science_coordinates.R
R/science_fits.R
R/science_kinematic_flags.R
R/segment_emission_lines.R
R/starlet_layer.R
R/structural_awareness.R
R/summarize_cluster_spectra.R
docs/SCIENCE_COORDINATES.md
docs/SCIENCE_PA_INVENTORY.txt
inst/extdata/kinematics/native_bisymmetric_workflow.R
inst/extdata/kinematics/native_kinematics_workflow.R
inst/tutorials/capivara_full_workflow.R
inst/tutorials/run_bisymmetric_bar_model.R
inst/tutorials/run_centa_magnum_emission_line_modes.R
inst/tutorials/run_centa_magnum_no_mask.R
inst/tutorials/run_centa_magnum_no_mask_preview.R
inst/tutorials/run_centa_magnum_starlet_comparison.R
inst/tutorials/run_centa_magnum_whole_vs_emission_windows.R
inst/tutorials/run_full_manga_science_workflow.R
inst/tutorials/run_kinematic_analysis.R
man/emission_lines.Rd
man/run_kinematic_analysis.Rd
man/run_manga_bar_model.Rd
man/segment_emission_lines.Rd
tests/testthat.R
tests/testthat/test-science-coordinates.R
tests/testthat/test-science-products.R
tests/testthat/test-segmentation-baseline.R
tests/testthat/test-workflow-contracts.R
tools/check_source_install.R
tools/validate_five_pilots.R
```

capivaraPPXF, relative to `3195657`:

```text
.Rbuildignore
.github/workflows/R-CMD-check.yaml
DESCRIPTION
R/diagnostics.R
R/emission.R
R/population.R
R/quality.R
R/science_baseline.R
R/spectra.R
inst/python/capivara_baseline.py
man/emission_quality_flags.Rd
man/fit_ppxf_population.Rd
man/fit_ppxf_population_diagnostics.Rd
tests/testthat/test-science-baseline.R
tools/check_source_install.R
```

This freeze report is an additional documentation-only file.

