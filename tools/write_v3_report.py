"""Write the V3 scientific decisions from retained evidence and exact source pins."""
from pathlib import Path
import csv,json,subprocess
import numpy as np
P=Path(__file__).resolve().parents[3];O=P/'results/science_baseline_freeze_v3';D=Path(__file__).resolve().parents[1]/'docs/SCIENCE_BASELINE_FREEZE_V3.md'
def read(p):return list(csv.DictReader(open(O/p)))
def table(rows,cols):
 out=['| '+' | '.join(cols.values())+' |','| '+' | '.join('---' for _ in cols)+' |']
 for row in rows:out.append('| '+' | '.join(str(row[k]) for k in cols)+' |')
 return '\n'.join(out)
pins=json.loads((O/'manifests/source_pins.json').read_text())
pilot=read('sandra_pilot/pilot_summary.csv');recovery=read('synthetic/stratified_recovery.csv')
strata=[x for x in recovery if x['stratum'] in ['snr_in','sigma']]
for x in strata:
 for key in ['log_age_bias','log_age_scatter','metal_bias','metal_scatter','velocity_bias','velocity_scatter','sigma_bias','sigma_scatter']:x[key]=f'{float(x[key]):.3f}'
for x in pilot:
 for k in ['median_snr','median_chi2']:x[k]=f'{float(x[k]):.2f}'
strategy=[];targets=read('resolution/target_strategies.csv')
for name in ['local_max','upper_envelope','galaxy_common','five_galaxy_common','robust95_unconstrained','local_max_template_floor']:
 a=[x for x in targets if x['strategy']==name];strategy.append(dict(strategy=name,ratio=f"{np.nanmedian([float(x['median_resolution_ratio']) for x in a]):.6f}",added=f"{np.nanmedian([float(x['median_added_variance']) for x in a]):.6f}",violations=sum(int(x['violating_samples']) for x in a)))
mixtures=read('synthetic/regional_mixtures.csv');mix=[]
for name in ['homogenized','median_scalar','rms_variable']:
 a=[x for x in mixtures if x['approximation']==name];mix.append(dict(method=name,n=len(a),bound=sum(x['bound_hit']=='True' for x in a),sigma_bias=f"{np.median([float(x['delta_sigma']) for x in a]):.3f}",sigma_p90=f"{np.quantile([abs(float(x['delta_sigma'])) for x in a],.9):.3f}",age_bias=f"{np.median([float(x['delta_log_age']) for x in a]):.4f}",metal_bias=f"{np.median([float(x['delta_metal']) for x in a]):.4f}"))
final=dict(REGIONAL_RESOLUTION_CONTRACT='CONDITIONAL',VARIANCE_CONTRACT='FAILED',
 STELLAR_POPULATION_RECOVERY='CONDITIONAL_TESTED_EMILES_DOMAIN',SANDRA_STELLAR_POPULATION_PILOT='NOT_READY',
 controlled_variance_subgate='DIAGONAL_APPROXIMATION_VALIDATED_UNDER_INDEPENDENT_NATIVE_NOISE',
 source_pins=pins,scope='Five frozen pilots, 90 regions. No population catalogue entries certified; no Sandra-46 production run.')
(O/'manifests/freeze_decisions.json').write_text(json.dumps(final,indent=2))
text=f'''# SCIENCE BASELINE FREEZE V3

*Resolution-aware regional spectroscopy and stellar-population inference*

V3 constructs a separate regional fitting spectrum with a tested instrumental response under a Gaussian LSF assumption. Its stellar-only backend recovers the specified E-MILES population summaries in controlled experiments. **The Sandra stellar-population pilot is NOT_READY.** Of 90 fixed regions, 62 yield numerical fits, 16 of those reach a parameter bound, and 28 lack sufficient valid wavelength support. No region passes the population QC gate. The native flux sums remain conserved; neither numerical convergence nor flux conservation establishes a calibrated population measurement.

| Decision | Verdict | Scientific scope |
| --- | --- | --- |
| REGIONAL_RESOLUTION_CONTRACT | CONDITIONAL | Gaussian POST response; uniform native grid; native sigma at least 0.9 pixel; local LSF variation at most 5% over the retained kernel neighbourhood. |
| VARIANCE_CONTRACT | FAILED | Sandra native spatial covariance is not propagated into the fitting variance. The independent-native-noise control separately passes DIAGONAL_APPROXIMATION_VALIDATED. |
| STELLAR_POPULATION_RECOVERY | CONDITIONAL_TESTED_EMILES_DOMAIN | The tested same-family, masked stellar-continuum calculation; details below. |
| SANDRA_STELLAR_POPULATION_PILOT | NOT_READY | All stellar catalogue quantities are withheld. |

## Preserved physical-coordinate baseline

The immutable runtime baseline is CAPIVARA 0.4.2.9000 at `4e6263af3fd0bfc556d4a2273c76e4ce535f7829`, capivaraPPXF 0.0.2.9000 at `cad2e5259f98f57570d39d0d515a0d9e35f1823c`, and SpectroPath 0.1.0 at `e1695030c50c09adf5eebdde56336c2342493c25`. The CAPIVARA V2 documentation head is `b8e5de79a858a7521fbfe2b2515409204c53a6dd`; its difference from the runtime pin contains only the report and two evidence tools.

Before source changes, the V1/V2 reports, validation scripts and retained evidence were read. The audit parsed 1,118 serialized, table or text objects, decoded 368 PNGs, verified all 347 final V2 artifact hashes and the source/report hashes, and recorded a before/after ledger for all 2,091 files in both evidence trees. V2 tests reproduced 471 CAPIVARA and 52 backend expectations, with one intentional backend skip. Exact V2 source archives and a separate baseline installation are retained under `baseline_archives/` and `baseline_library/`. Neither historical evidence directory was regenerated. Development uses the isolated `fix/science-baseline-v3` branches.

The five native inputs and the saved V2 rest-frame partition maps define 18 regions per object. No segmentation, bar detector, morphology, SpectroPath, hierarchical analysis or environment calculation was rerun. The V2 observed/rest, vacuum/air and systemic-relative coordinate contracts remain unchanged. Native gas moments and SpectroPath widths remain observed widths with no instrumental correction. The existing joint stellar-plus-gas native-LSF rejection is unchanged.

## The two regional spectra

The native product is the finite native FLUX sum, $F_{{R,\mathrm{{native}}}}(\lambda)=\sum_{{p\in R}}F_p(\lambda)$, on the complete 6,732-channel observed-vacuum WAVE vector. It retains native counts and the variance sum conditional on independent contributors; missing contributor variance remains missing. A Gaussian effective FWHM is not assigned to this sum. The maximum relative L1 discrepancy from the frozen V2 sums is $1.951\times10^{{-15}}$ over all 90 regions. Flux units are the native $10^{{-17}}$ erg s$^{{-1}}$ cm$^{{-2}}$ Angstrom$^{{-1}}$ units summed over spaxels.

For the fitting product, each contributing spaxel is smoothed before summation. The adopted target is the per-region, per-wavelength maximum POST sigma, raised only where necessary to reach the empirical template resolution after the coordinate transformation. Its additional kernel obeys

$$\sigma_{{\mathrm{{conv}},p}}^2(\lambda)=\sigma_{{\mathrm{{target}}}}^2(\lambda)-\sigma_p^2(\lambda).$$

The implementation uses a column-normalized variable Gaussian operator, truncated at six sigma. Its discrete second moment is calibrated to the requested variance, including subpixel broadening; sigma zero gives exact identity. Column normalization conserves the unmasked integrated flux on a uniform grid. A finite sampled kernel is not itself an exact Gaussian. The independent response experiment, rather than the formula alone, determines the permitted approximation.

The fitting aperture remains the entire fixed region. A wavelength is accepted only when every contributor and its kernel support have finite flux, positive variance, eligible MASK values and valid native LSF. Invalid LSF values are not interpolated. Edge halos and neighbourhoods with more than 5% LSF variation are excluded; finite native sums are retained even when the fitting spectrum is wholly unusable. This conservative policy avoids silently changing the population mixture with wavelength. It can reject regions containing permanently unusable spaxels; any future fit to a smaller aperture must receive a distinct aperture definition.

The pilot uses `MASK == 0`. The observed rejected value 1026 includes LOWCOV and DONOTUSE. It must not be reinstated merely to obtain a fit. [SDSS DRP3PIXMASK definitions](https://www.sdss4.org/dr14/algorithms/bitmasks/).

## Sampling, medium and template resolution

The stellar model treats the sampled empirical spectrum as point sampled and introduces no separate detector-pixel top-hat. The E-MILES empirical FWHM and the MaNGA POST estimates therefore describe the responses being matched in this calculation. pPXF 9.4.5 applies its analytic Fourier LOSVD to the sampled templates; rebinning to the fitting grid is applied to both data and templates. This is a method-specific Gaussian-response approximation, tested by the same-family recovery experiment. A gas template with explicit pixel integration would require a separate PRE contract and is not enabled here. [MaNGA method-dependent LSF guidance](https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/), [pPXF documentation](https://pypi.org/project/ppxf/9.4.5/).

The native data remain observed vacuum wavelengths. The fitter computes rest vacuum coordinates using the supplied systemic redshift, then converts every sample to air with the explicit SDSS Morton relation used in V2. The flux-density Jacobian and its square for variance accompany this coordinate change. The inverse conversion is tested. No constant air/vacuum scale approximation is used. The installed pPXF population example explicitly identifies the E-MILES/MILES wavelength medium as air; the local template archive supplies its wavelength-dependent FWHM array.

The exact archive is `spectra_emiles_9.0.npz`, SHA256 `6a1b1a70db4d95f67afe0adbe254f111eb2b50412b9f311714c54ffba18c6f0e`. The backend refuses another archive. Its optical sampling is 0.9 Angstrom and FWHM is 2.51 Angstrom over the retained 3550--7800 Angstrom air template interval. Data narrower than the templates are rejected; negative broadening variance is not silently clipped into compatibility.

The fit uses 3700--7400 Angstrom rest air, one paired 90 km/s logarithmic grid, moments 2, no additive polynomial, multiplicative degree 6, no regularization, and no iterative residual clipping. Stellar velocity is the pPXF logarithmic velocity relative to the supplied systemic redshift; the optical equivalent is $c[\exp(v/c)-1]$. Velocity bounds are ±500 km/s and sigma bounds are 0.01--400 km/s. The lower bound only keeps the Fourier calculation nonsingular; it is not a physical resolution limit. A zero bound had produced seven failed least-squares evaluations in the original low-S/N control; the final run retains these as explicit bound hits.

Fourteen V2 registry lines receive ±800 km/s stellar-fit masks after explicit conversion to air. Observed sky masks are transformed separately. There are no gas templates, extinction parameter, mass estimate or M/L calculation. The 36-template basis contains ages 0.1, 0.2512, 0.631, 1, 1.9953, 3.9811, 6.3096, 10 and 12.5893 Gyr, and template metallicities −0.71, −0.4, 0 and +0.22. Weights are normalized in 5070--5950 Angstrom air. Reported summaries are the luminosity-weighted mean log10(age/year) and the luminosity-weighted arithmetic mean template [M/H]. They do not identify a detailed star-formation history or mass-weighted age.

## Target-resolution cost and independent response tests

The table reports the median, across the 90 regions, of the wavelength-median resolving-power ratio relative to the native regional median sigma and the corresponding added variance in Angstrom squared. The five-galaxy/galaxy comparisons use the existing observed grid; equal observed-grid targets across galaxies do not assert equal rest-frame responses.

{table(strategy,{'strategy':'Target strategy','ratio':'Median R ratio','added':'Added variance [Angstrom²]','violations':'Samples below a contributor'})}

The unconstrained 95th percentile fails the physical inequality. Constraining a robust envelope by the actual maximum removes those violations but adds broadening. The regional maximum plus template floor retains a median resolving-power ratio of 0.996889; a galaxy-wide target gives 0.980522 and the five-object observed-grid target 0.939162. The latter is unnecessary for this experiment and is not adopted. Missing LSF coverage remains recorded rather than used to infer an envelope.

The independent test integrates intrinsic Gaussian lines through a fine-grid variable response, then applies the production operator only to the sampled output. It covers five LSF shapes (constant, smooth, sharp, and one actual profile from each DRP generation), two spatial scale factors, six line centres and intrinsic widths 0.15, 1 and 5 Angstrom: 180 cases. Widths are checked in Angstrom and km/s. The largest fitted-width error is 0.5641%, below the 2% tolerance, with flux conservation to floating-point precision. Missing LSF, edge rejection, exact no-op, target inequalities and both native extension conventions have separate regressions. The line experiment validates this approximation over its tested sampling and response domain; it is not a calibration of the DRP LSF estimates themselves.

![Independent response recovery]({O}/figures/resolution_recovery.png)

*Fractional fitted-width error relative to the independently forward-modeled target. Different symbols at the same wavelength represent spatial response factors and intrinsic widths. The dotted ±2% limits are the declared numerical tolerance.*

## Variance and covariance

For independent native errors, the retained spectral covariance is the sum of $A_p\,\mathrm{{diag}}(V_p)\,A_p^T$ across contributing spaxels. The rest-air Jacobian, common normalization and log-rebin operator are composed with that covariance before taking its diagonal. Propagating the rebinned diagonal alone would discard correlations already produced by homogenization. No empirical scalar noise vector replaces the input uncertainty.

Twenty thousand noise realizations at each of five convolution widths verify the propagated diagonal: the largest median relative variance discrepancy is 0.1944%. At smoothing widths 0.2, 0.5, 1 and 2 native pixels, lag-one correlations are 0.0416, 0.3100, 0.7786 and 0.9394. Ignoring them understates the error of the summed spectral measurement by factors 1.041, 1.291, 1.872 and 2.632 in these controls. These are measurement-specific factors, not a prescription to rescale all fit errors.

An initial log grid oversampled native samples and produced singular covariance matrices. The final 90 km/s grid avoids that rank defect without diagonal jitter on science pixels. Eighteen paired full-induced-covariance and diagonal-objective fits cover both DRP profiles and sigma 20, 80 and 200 km/s at nominal native S/N 30. Median absolute changes are 0.00128 dex in log age, 0.00227 dex in [M/H], 0.409 km/s in velocity and 1.002 km/s in sigma; the largest sigma change is 7.864 km/s. This comparison supports the diagonal objective within the controlled domain, conditional on the independent-native-noise premise.

A 20-draw correlated-noise bootstrap of the 1 Gyr, solar-metallicity, 80 km/s, S/N 30 control gives standard deviations 0.0232 dex in log age, 0.0253 dex in metallicity, 2.052 km/s in velocity and 3.466 km/s in sigma. These are conditional on the fitted spectral model and input covariance. They are not model-systematic errors or a coverage validation. QC-rejected pilot fits receive missing uncertainties, with the method explicitly recorded as not estimated.

The Sandra containers include GCORREL, RCORREL, ICORREL and ZCORREL. A separate 360-region/band diagnostic uses their sparse spatial correlations and the native variances for MASK==0 spaxel subsets. The median correlated-to-independent noise ratio is 2.688, spanning 1.000--3.887. These tables are therefore materially relevant to the regional uncertainty. Their BBINDEX coordinates do not match the native linear WAVE vectors: the maximum discrepancy from BBWAVE is 1396.71 Angstrom. The diagnostic selects the nearest native wavelength to BBWAVE and retains both coordinates, index conventions and subset counts. It does not assert that the four tables reconstruct the full joint spatial/spectral covariance. [SDSS covariance and table-coordinate definitions](https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/).

The fitting products still assume independent native spaxels and spectral samples before the newly modeled convolution. Consequently the **Sandra variance contract fails**, even though the operator calculation and its controlled diagonal approximation pass. Reduced chi-square cannot distinguish template inadequacy from this missing noise calibration. No fit is certified by dividing its chi-square by a generic spatial factor.

## Controlled stellar recovery

The final grid contains **810 same-family fits**: both DRP profiles; ages 0.1, 1 and 10 Gyr; metallicities −0.71, 0 and +0.22; sigma 20, 80 and 200 km/s; nominal native regional S/N 10, 30, 60, 150 and 300; three noise realizations per cell. The three realizations also use velocities −120, 0 and +120 km/s. Thus the noise and velocity effects are not independently crossed. Each synthetic region contains three spaxels with different instrumental responses. The high-S/N extension covers the formal S/N encountered in the pilots; the controlled grid uses the real native wavelength spacing.

Truth is generated from the same empirical SSP family on an independent fine logarithmic grid, with an injected Gaussian LOSVD, followed by a fine observed-wavelength instrumental calculation and native sampling. Cubic interpolation prevents the extra artificial broadening found in an earlier linear-interpolation attempt. That attempt and the original log-grid and zero-bound attempts remain under `synthetic/attempt_*`. Neither the final generator nor its declared truth invokes pPXF.

The final run has 803 interior numerical solutions and seven lower-sigma bound hits, all in the low-S/N, low-dispersion regime. For each marginal stratum, the table gives median bias and robust scatter (1.4826 times the median absolute deviation), including bounded solutions. Full age, metallicity and DRP stratification, 90th-percentile errors and every individual trial are retained in CSV.

{table(strata,{'stratum':'Stratum','value':'Value','log_age_bias':'Δlog age','log_age_scatter':'Scatter','metal_bias':'Δ[M/H]','metal_scatter':'Scatter','velocity_bias':'Δv [km/s]','sigma_bias':'Δsigma [km/s]','sigma_scatter':'Scatter'})}

The conditional acceptance subset is native S/N 30--300, injected sigma 80--200 km/s, the tested 0.1--10 Gyr ages and −0.71 to +0.22 metallicities. Every age, metallicity, DRP, dispersion and S/N marginal within that subset passes the declared bias/scatter tolerances. This does not certify all intermediate SSP mixtures, arbitrary abundance patterns, arbitrary LSF curves or lower dispersions. Three realizations per exact cell are insufficient to calibrate detailed confidence-interval coverage. Sigma 20 km/s remains unsupported: positive-bound effects and residual sampling sensitivity are large relative to the injected value. Catalogue sigma is withheld outside the tested response/SNR regime and for any QC-rejected fit.

![Conditional stellar recovery]({O}/figures/synthetic_recovery.png)

*Median recovery error and the 16th--84th-percentile interval over the tested ages, metallicities, DRP profiles and realizations at each nominal S/N. Colours distinguish injected stellar dispersions. These are controlled errors relative to known inputs, not uncertainties measured for Sandra.*

## Regional mixtures and template mismatch

The mixture experiment contains 72 trials per treatment, combining both DRP profiles, young/old starting populations, sigma 20/80/200 km/s, three realizations and identical-population versus mixed-population spaxels. The mixed case combines the starting SSP with 1 and 10 Gyr SSPs with different responses. Population truth uses the common luminosity-normalization band. The median scalar approximation changes the spectral resolution model while retaining the same noisy regional input.

{table(mix,{'method':'Treatment','n':'Trials','bound':'Bound hits','sigma_bias':'Median Δsigma [km/s]','sigma_p90':'90th |Δsigma| [km/s]','age_bias':'Median Δlog age','metal_bias':'Median Δ[M/H]'})}

The scalar approximation drives 24 of 72 fits to the lower dispersion bound, whereas the homogenized and wavelength-dependent RMS approximations have no bound hits in this experiment. Homogenization gives a median sigma error of +1.058 km/s; the scalar treatment gives −11.400 km/s. The RMS approximation happens to perform comparably here. It remains an approximation to a response mixture, rather than an exact Gaussian response for arbitrary populations and spaxel weights. The experiment demonstrates a material scalar-resolution error without claiming that homogenization dominates every possible approximation.

A separate 48-fit perturbation experiment adds 4% localized absorption responses near 5175 and 5270 Angstrom to the generating model. It compares perturbed and unperturbed 1 and 10 Gyr, solar-metallicity spectra at sigma 80/200 km/s in both DRP profiles. This isolates one deliberate abundance-like spectral mismatch; it is not a different-library validation or a general stellar-population systematic-error budget. The perturbations and their paired parameter changes are retained under `template_mismatch/`.

## Five-pilot result and catalogue limits

All fitting attempts use the frozen rest-frame partitions. Numerical-success counts below include bound hits. The 28 insufficient-support regions remain in every inventory; they were not dropped or relabelled. Formal S/N and chi-square use the spatially independent native-noise premise and are not calibrated Sandra uncertainties.

{table(pilot,{'galaxy_id':'Galaxy','n_regions':'Regions','n_numerical_success':'Numerical fits','n_bound':'Bound hits','n_insufficient':'Insufficient support','median_snr':'Median formal S/N','median_chi2':'Median chi-square','n_population_usable':'Certified populations'})}

Bound hits affect 17.8% of all regions and 25.8% of numerical fits; insufficient support affects 31.1% of regions. There are no remaining numerical exceptions in the final pilot. Population certification is 0/90. Interior fits still fail the residual/noise-calibration gate, and some also place most luminosity weight on a template-grid edge. Catalogue age, metallicity, stellar velocity and dispersion columns are therefore missing; the raw fit values are explicitly retained as diagnostics. No placeholder uncertainties are assigned.

Safe outputs of this freeze are the fixed region identities and maps, native measured sums with their validity/count information, the separate controlled spectra, LSF and wavelength provenance, transformation operators/covariance under their stated assumptions, and fit/QC records. Stellar ages, metallicities, velocities, intrinsic dispersions, SFH detail, mass-weighted ages, masses, M/L and extinction are **not certified for the Sandra catalogue**. No intrinsic gas dispersion or deprojected dynamics is inferred.

![Pilot diagnostics]({O}/figures/pilot_diagnostics.png)

*All numerical pilot fits, including bound hits, are shown as QC-rejected diagnostics. Grey crosses are not certified population measurements. Fractions use all 18 regions per object; no-fit counts include insufficient support. The chi-square threshold is displayed under the stated independent-native-noise assumption.*

Each of the five `figures/*_spectroscopy.png` figures shows the categorical fixed partition with region identifiers/boundaries, a representative successful numerical fit, an insufficient-support region where present, native and controlled spectra, the stellar model and residuals on accepted pixels. Emission features outside the stellar masks remain visible in the data but do not enter the displayed residual sample.

## API, source and verification

CAPIVARA exposes `prepare_segment_spectra()` for the physically distinct products; it requires the contributing spaxel arrays because summed spectra cannot reconstruct their responses. The numerical operator is an internal Python helper using NumPy/SciPy through reticulate. The result is serializable in R, including sparse spectral covariance. capivaraPPXF exposes `fit_segment_spectra()` for the experimental stellar-only consumer. Its output retains numerical status, bounds, residuals, masks, template metadata, transformed grids, normalization and the uncertainty method. It does not turn numerical success into population certification. Research grids, target comparisons and pilot diagnostics remain under `tools/` and are excluded from package builds.

| Component | Version | V3 code commit |
| --- | --- | --- |
| CAPIVARA | 0.4.3.9000 | `{pins['capivara']}` |
| capivaraPPXF | 0.0.3.9000 | `{pins['capivaraPPXF']}` |
| SpectroPath, unchanged | 0.1.0 | `e1695030c50c09adf5eebdde56336c2342493c25` |

Changed runtime files are CAPIVARA `R/science_regional_spectra.R` and `inst/python/capivara_resolution.py`, and backend `R/stellar.R` and `inst/python/capivara_stellar.py`, plus their DESCRIPTION, NAMESPACE and generated help entries. Two R regression files and two Python regression scripts cover the added contract. The existing wavelength, medium, PRE/POST, gas guard and one-pixel grid regressions remain active.

The final CAPIVARA testthat run passes 481 expectations; capivaraPPXF passes 55, with its one intentional missing-pPXF-branch skip because pPXF is installed. The standalone suites pass eight response and four backend test cases, including controlled stellar recovery, bound retention, medium round trips, template incompatibility, covariance rank and invalid-support retention. Both source archives pass `R CMD build` and `R CMD check --no-manual` with zero errors, warnings or notes in the final check status. Normal and isolated installations are compared against all source function bodies, exports and bundled scripts; the exact results and final sessions are in `tests/` and `manifests/`.

The fitting runtime is Python 3.9, pPXF 9.4.5, NumPy 1.26.4 and SciPy 1.13.1; R is 4.5.2 on arm64 macOS. Independent FITS/figure/response tools use the separately recorded Python 3.12 environment. A NumPy deprecation warning inside the pinned pPXF `log_rebin` implementation appears in the standalone unittest log; it does not change the current calculation and is retained rather than hidden. Initial package-check failures caused by omitting the normal dependency library were corrected. Offline repository-index diagnostics do not change the final clean check status. Exact environments, seeds, archive/input hashes and source/install comparisons are retained.

Before the 46-object run, the spatial correlation tables and their wavelength mapping must be incorporated into a validated regional covariance calculation, including its interaction with different spaxel convolution operators. The fitting aperture must be specified for regions with permanently unusable contributors. Only after those changes should the residual/model adequacy, parameter bounds, template-family sensitivity and uncertainty calibration be reconsidered. V3 performs no Sandra-46 production, environment comparison, morphology revision, merge or push.

Evidence root: `{O}`. The final manifest records code/report pins, versions, all input and template hashes, validation artifacts and the unchanged V1/V2 ledger. `synthetic/configuration.json`, `synthetic/domain_gates.csv`, `variance/native_spatial_correlation_bands.csv`, `sandra_pilot/all_regions.csv` and `manifests/freeze_decisions.json` provide the machine-readable decisions. **Stop after V3; no population-baseline readiness flag is issued.**
'''
D.write_text(text)
print(D)
