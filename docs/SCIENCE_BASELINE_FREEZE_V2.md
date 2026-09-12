# CAPIVARA science-baseline freeze V2: physical wavelength and LSF contracts

Validated 2026-09-11. The [V1 freeze](SCIENCE_BASELINE_FREEZE.md) and all 766 historical evidence files are unchanged. This release validates observed/rest wavelength, wavelength medium, systemic-relative velocity and native wavelength-dependent LSF as separate contracts.

```
READY_FOR_RESTFRAME_SPECTRAL_SEGMENTATION
READY_FOR_SYSTEMIC_FRAME_KINEMATIC_SEGMENTATION
READY_FOR_SANDRA_SPATIAL_ATLAS_BASELINE
```

These statements cover fixed spectral partitions, systemic coordinates, observed kinematic/path descriptors and the spatial-atlas baseline. They do not certify intrinsic gas dispersions, stellar populations, deprojected dynamics or independent morphological accuracy. Native-LSF input is now explicitly rejected by the legacy joint pPXF fitter because its regional template-resolution contract remains incomplete. V1’s broader `NOT_READY_FOR_SANDRA_PILOT` inference decision is unchanged.

## Defect and retained historical interpretation

Before modification, both `segment()` and `segment_large()` passed the requested wavelength interval unchanged to the native-axis subsetting helper. Their `redshift` argument did not transform that interval. On a 4 × 4 cutout of 11004-12701, z = 0.018372238, both selected channels 1179–3779: 2601 samples at 4800–7400 Å observed, corresponding to 4713.40421595428–7266.49816626286 Å rest. A 4800–7400 Å rest request instead maps to observed bounds 4888.1867424–7535.9545612 Å and selects the available samples 4889–7535 Å.

The native FLUX header identifies MaNGA and `CTYPE3=WAVE`; its axis agrees with HDU WAVE, with 6732 channels from 3622 to 10353 Å at 1 Å spacing. These Sandra files contain the linear native DRP cubes inside MEGACUBE containers. WAVE is vacuum and heliocentric; it has not had the galaxy’s systemic redshift removed. [MaNGA DRP, Table B8](https://www.sdss4.org/wp-content/uploads/2016/07/law_arxiv.pdf).

V1’s description “4800–7400 Angstrom rest-frame feature window” was semantically incorrect for the full-spectrum partitions. Its statement that no scientific segmentation change was needed established reproducibility of the observed-window calculation after earlier coordinate repairs. It did not establish equivalence to a rest-window calculation. Those numerical partitions remain valid historical observed-frame results. The separately specified pPXF fitting interval was a different operation.

`pre_edit_proof.rds` and `pre_edit_selected_channels.csv` retain the original demonstration. `tools/verify_historical_wavelength_defect.R` replays it through the preserved V1 installation. No historical product was relabelled or overwritten.

## 1. Observed and rest wavelength

Both spectral APIs require an explicit `feature_wavelength_frame` for a bounded request. `wavelength_frame` can be an argument or input metadata; conflicting declarations fail. Rest-frame requests require a finite scalar redshift greater than −1, including for an already-rest-frame input. An explicit observed-frame request on native observed data is independent of redshift. Ambiguous historical calls fail with an explanation of the old native-axis behavior.

```r
seg <- capivara::segment_large(
  cube, Ncomp = 18L, redshift = z,
  wavelength_frame = "observed",
  feature_wavelength_range = c(4800, 7400),
  feature_wavelength_frame = "rest",
  knn_k = 40L, auto_k = TRUE, max_k = 100L,
  feature_scale = "robust_col", spatial_weight = 0.15,
  mask = fixed_galaxy_support, valid_mode = "signal"
)
seg$wavelength_provenance
```

On the native vacuum grid, rest bounds are multiplied by 1 + z and the inclusive native channels are selected. Neither the cube nor its flux samples are interpolated. Already-rest-frame input does not receive a second correction. An explicitly declared air axis is converted through vacuum when changing frame. The retained original cube supplies full native regional sums.

Every segmentation records `input_wavelength_frame`, `requested_feature_wavelength_range`, `requested_feature_wavelength_frame`, `systemic_redshift`, `selected_native_wavelength_min/max`, `selected_rest_wavelength_min/max`, `selected_channel_indices`, and `number_of_selected_channels`. Medium, one-based indexing and the absence of resampling are explicit. Variance selection uses the same channel indices. Summed spectra retain the selection provenance.

The table gives actual selected samples, not the continuous requested endpoints. Native wavelengths and rest wavelengths are in Å; channel indices are one-based. The complete numeric values are in [exact_windows.csv](../../../results/science_baseline_freeze_v2/physical/exact_windows.csv).

| Object | z | Requested frame | Native limits | Corresponding rest limits | Channels | N |
| --- | --- | --- | --- | --- | --- | --- |
| 11004-12701 | 0.018372238 | observed | 4800–7400 | 4713.404215954–7266.498166263 | 1179–3779 | 2601 |
| 11004-12701 | 0.018372238 | rest | 4889–7535 | 4800.798585792–7399.062659837 | 1268–3914 | 2647 |
| 11014-3704 | 0.024539573 | observed | 4800–7400 | 4685.031331630–7222.756636263 | 1179–3779 | 2601 |
| 11014-3704 | 0.024539573 | rest | 4918–7581 | 4800.205018533–7399.421359393 | 1297–3960 | 2664 |
| 8602-12705 | 0.0318464 | observed | 4800–7400 | 4651.855159838–7171.610038083 | 1179–3779 | 2601 |
| 8602-12705 | 0.0318464 | rest | 4953–7635 | 4800.133043058–7399.357113617 | 1332–4014 | 2683 |
| 8932-3701 | 0.0255681 | observed | 4800–7400 | 4680.332783362–7215.513041016 | 1179–3779 | 2601 |
| 8932-3701 | 0.0255681 | rest | 4923–7589 | 4800.266310935–7399.801144361 | 1302–3968 | 2667 |
| 9869-9102 | 0.0278487 | observed | 4800–7400 | 4669.948018614–7199.503195363 | 1179–3779 | 2601 |
| 9869-9102 | 0.0278487 | rest | 4934–7606 | 4800.317400800–7399.921797829 | 1313–3985 | 2673 |
| 8941-12704 | 0.0139907 | observed | 4800–7400 | 4733.771226896–7297.897308131 | 1179–3779 | 2601 |
| 8941-12704 | 0.0139907 | rest | 4868–7503 | 4800.832985944–7399.476149042 | 1247–3882 | 2636 |
| 9039-9101 | 0.1253 | observed | 4800–7400 | 4265.529192215–6576.024171332 | 1179–3779 | 2601 |
| 9039-9101 | 0.1253 | rest | 5402–8327 | 4800.497645072–7399.804496579 | 1781–4706 | 2926 |

## Observed-window versus rest-window partitions

Exactly five pilots and the two redshift extremes were compared. Reading the complete 46-object manifest identified 8941-12704 (z = 0.0139907) and 9039-9101 (z = 0.1253); no 46-object analysis was run. The spectral comparison held Ncomp = 18, graph settings, scaling, spatial weight and support fixed. Full native flux with the historical detector parameters defined one shared support per galaxy. All five observed-window partitions and masks reproduce V1.

| Object | ARI | VI [nats] | Boundary Jaccard | Spatial components, observed → rest |
| --- | --- | --- | --- | --- |
| 11004-12701 | 0.5710 | 1.0879 | 0.6201 | 34 → 35 |
| 11014-3704 | 0.8179 | 0.5669 | 0.6425 | 24 → 22 |
| 8602-12705 | 0.8449 | 0.5468 | 0.5278 | 32 → 34 |
| 8932-3701 | 0.8089 | 0.5103 | 0.6875 | 18 → 19 |
| 9869-9102 | 0.6659 | 1.0789 | 0.5756 | 35 → 33 |
| 8941-12704 | 0.8795 | 0.3139 | 0.7993 | 26 → 26 |
| 9039-9101 | 0.6816 | 0.7272 | 0.5520 | 49 → 45 |

VI uses natural logarithms. Boundary agreement is the Jaccard index of changed-label edges on the common four-neighbour grid. Fragmentation counts four-connected components across the 18 labels; diagonal contact does not connect components. Every comparison retains the same assigned footprint and exactly 18 labels. The differences are not monotonic in redshift: the smallest pilot ARI is 0.571, while the low-redshift leverage object has ARI 0.879. Galaxy spectra and spatial structure affect which boundaries move.

![Observed and rest partitions](../../../results/science_baseline_freeze_v2/validated/figures/11004-12701_partitions.png)

*Partitions of 11004-12701 under the historical observed interval and the requested rest interval. Both calculations use the same support and graph settings. Colour labels are matched by spatial overlap; agreement is not imposed.*

Regional comparisons use a maximum-overlap one-to-one assignment of labels. A forced pair with zero spatial overlap is flagged and must not be interpreted as the same physical region. The summed-spectrum discrepancy is the L1 difference divided by the historical sum’s L1 norm.

| Object | Median \|Δsize\| [spaxels] | Median summed-spectrum L1 | Median L1, overlapping pairs | Max \|Δbar fraction\|, overlapping | Zero-overlap pairs |
| --- | --- | --- | --- | --- | --- |
| 11004-12701 | 9.0 | 0.1186 | 0.0688 | 0.0823 | 1 |
| 11014-3704 | 2.5 | 0.1338 | 0.1338 | 0.2368 | 0 |
| 8602-12705 | 3.0 | 0.1381 | 0.0984 | 0.0126 | 1 |
| 8932-3701 | 1.5 | 0.0796 | 0.0796 | 0.1848 | 0 |
| 8941-12704 | 0.0 | 0.0000 | 0.0000 | 0.0453 | 0 |
| 9039-9101 | 4.0 | 0.3035 | 0.1737 | 0.1667 | 4 |
| 9869-9102 | 1.0 | 0.0469 | 0.0285 | 0.1080 | 1 |

The per-region sizes, normalized spectral-shape changes, bar fractions, overlap counts and component sizes are retained alongside the full native regional sums in `validated/comparisons/`. The highest-redshift object has several disconnected bright sources in the existing detector support; its comparison does not certify catalogue membership or a unique-galaxy mask. No support was tuned for agreement.

For all 14 partitions, the regional sums reproduce the direct sum of the same assigned native spaxels at every wavelength. Both maximum absolute and relative L1 conservation discrepancies are numerically zero. This conserves retained measured flux; it does not restore missing measurements or calibrate their uncertainties. All five bar-rotation checks retain IoU = 1, zero area change and PA discrepancy below 10⁻⁶ degrees.

![Native regional sums](../../../results/science_baseline_freeze_v2/validated/figures/11004-12701_summed_spectra.png)

*Full native summed spectra for the four largest historical regions and their spatial-overlap counterparts. Blue denotes observed-window segmentation and orange rest-window segmentation. Flux remains in native FLUX units; changes in region size and spectral mixture remain visible.*

## 2. Vacuum and air wavelength

All laboratory constants in CAPIVARA’s R code, native kinematic workflow and bundled full-workflow tutorial now come from one referenced registry. SpectroPath itself contains no laboratory line catalogue: its public functions consume the coordinate supplied by CAPIVARA. The registry exposes identifier, label, laboratory wavelength, medium, source reference, registry version and conversion provenance through `emission_lines()`.

The pinned reference is **SDSS DR17 DAP Table 2**. It lists Hα = **6564.608 Å** and [O III] = **5008.240 Å** in vacuum. The user’s approximate values, 6564.614 and 5008.239 Å, appear in the SDSS Classic table; the earlier Sandra configuration used the DR15 Hα value 6564.632 Å. These are distinct published reference choices. Relative to DR17, the Hα differences are about 0.274 and 1.096 km/s, respectively. [DR17 DAP definitions](https://www.sdss4.org/dr17/manga/manga-analysis-pipeline/), [SDSS Classic definitions](https://classic.sdss.org/dr5/products/spectra/vacwavelength.php), [DR15 DAP definitions](https://www.sdss4.org/dr15/manga/manga-analysis-pipeline/).

| Identifier | Laboratory vacuum wavelength [Å] |
| --- | --- |
| oii3726 | 3727.0920 |
| oii3729 | 3729.8750 |
| neiii3869 | 3869.8600 |
| hdelta | 4102.8922 |
| hgamma | 4341.6837 |
| hbeta | 4862.6830 |
| oiii4959 | 4960.2950 |
| oiii5007 | 5008.2400 |
| oi6300 | 6302.0460 |
| halpha | 6564.6080 |
| nii6548 | 6549.8600 |
| nii6583 | 6585.2700 |
| sii6716 | 6718.2950 |
| sii6731 | 6732.6740 |

The former [O II] blend entry was not a unique laboratory transition. The registry now separates 3726 and 3729; ambiguous `oii`/`oii3727` requests fail. The “strong” line group explicitly includes both transitions. No independent air table can drift away from the vacuum registry. `emission_lines("air")` explicitly converts the same entries with the SDSS-published Morton approximation; `convert_wavelength_medium()` uses its iterative inverse for air → vacuum and is tested for round trips. The Classic tabulated pairs differ from that page’s conversion approximation by up to 0.003 Å, so the tables are not treated as exact inverse pairs.

A laboratory wavelength must match its declared registry entry and medium. The physical helper rejects mismatched line/input media and rejects an air numerical value falsely labelled vacuum. Native MaNGA metadata also rejects an explicit air declaration. Regression profiles at z = 0.02 and 0.1253 fail if common air wavelengths are applied directly to vacuum WAVE: the induced systemic-zero bias is approximately 83 km/s. The air-coordinate tests explicitly convert through vacuum before evaluating redshift and optical velocity.

## 3. Systemic-relative velocity

For vacuum laboratory wavelength λ₀ and systemic redshift z, CAPIVARA evaluates the native samples using

\[
\lambda_{\rm sys}=\lambda_0(1+z),\qquad
v=c\left(\frac{\lambda_{\rm obs}}{\lambda_{\rm sys}}-1\right),\qquad
\lambda_{\rm obs}=\lambda_{\rm sys}(1+v/c),
\]

with c = 299792.458 km/s. This is the optical velocity convention. An already-rest-frame axis uses λ₀ in the denominator. The full cube is not rebinned. Native conventional features preserve systemic velocity; median-centred velocity and its subtracted offset are separate outputs. Scaling features does not redefine the physical origin.

The systemic profile is X(v) = (v, F(v)). `profile_centering_mode="local_centroid"` subtracts the measured positive-profile centroid from the same extracted coordinate. Extraction remains centred on the systemic wavelength. SpectroPath’s translation-invariant signatures are unchanged by translating those same samples, while the retained centroid encodes the removed velocity. These properties are tested through `as_spectral_path()`, `line_moments()`, `classical_features()` and `path_features()`.

Kinematic provenance records `line_rest_wavelength`, `wavelength_medium`, `input_wavelength_medium`, `systemic_redshift`, `observed_line_centre`, line identifier/reference, redshift source, optical velocity definition, extraction half-width and centring mode. Both systemic and local-centred representations retain these fields.

Every line below uses the registry’s laboratory vacuum wavelength and a requested velocity interval of **[−600, +600] km/s**. The sampled interval generally stops inside those limits. The last column checks wavelength → velocity → wavelength over the full native WAVE vector.

| Object | Line | Predicted observed λ [Å] | Extracted native λ [Å] | Sampled v [km/s] | Max round-trip \|Δλ\| [Å] |
| --- | --- | --- | --- | --- | --- |
| 11004-12701 | Hα | 6685.214541 | 6672–6698 | -592.594 to 573.352 | 9.09e-13 |
| 11004-12701 | [O III] | 5100.252577 | 5091–5110 | -543.866 to 572.953 | 1.82e-12 |
| 11014-3704 | Hα | 6725.700677 | 6713–6739 | -566.122 to 592.806 | 9.09e-13 |
| 11014-3704 | [O III] | 5131.140071 | 5121–5141 | -592.445 to 576.077 | 1.82e-12 |
| 8602-12705 | Hα | 6773.667132 | 6761–6787 | -560.628 to 590.093 | 9.09e-13 |
| 8602-12705 | [O III] | 5167.734414 | 5158–5178 | -564.716 to 595.531 | 1.82e-12 |
| 8932-3701 | Hα | 6732.452554 | 6719–6745 | -599.035 to 558.731 | 9.09e-13 |
| 8932-3701 | [O III] | 5136.291181 | 5127–5146 | -542.303 to 566.679 | 1.82e-12 |
| 9869-9102 | Hα | 6747.423799 | 6734–6760 | -596.428 to 558.769 | 9.09e-13 |
| 9869-9102 | [O III] | 5147.712973 | 5138–5158 | -565.664 to 599.096 | 1.82e-12 |
| 8941-12704 | Hα | 6656.451461 | 6644–6669 | -560.787 to 565.160 | 9.09e-13 |
| 8941-12704 | [O III] | 5078.308783 | 5069–5088 | -549.534 to 572.110 | 1.82e-12 |
| 9039-9101 | Hα | 7387.153382 | 7373–7401 | -574.386 to 561.937 | 1.36e-12 |
| 9039-9101 | [O III] | 5635.772472 | 5625–5647 | -573.037 to 597.243 | 9.09e-13 |

The public `segment_kinematics(..., segmentation_mode="path_signature")` smoke runs used two pilots from different DRP generations, both lines, 18 conventional-feature regions, 45 path-signature regions and the fixed 40-neighbour graph.

| Object | Line | Measured profiles | Assigned conventional spaxels | Path profiles | Median systemic v [km/s] |
| --- | --- | --- | --- | --- | --- |
| 11004-12701 | halpha | 1133 | 1133 | 1133 | 1.897 |
| 11004-12701 | oiii5007 | 1133 | 1133 | 1133 | 1.012 |
| 8602-12705 | halpha | 423 | 352 | 423 | -33.678 |
| 8602-12705 | oiii5007 | 720 | 716 | 720 | 36.728 |

Conventional features remain the existing flux proxy, centroid, observed width, asymmetry and third/fourth-moment proxies; path features remain `p2`, `p3u`, `p3F`, `p4F`, `p4T`, and `p_pm`. Undefined higher moments explain the smaller conventional assignment counts. Missing profiles are not zero-filled or spatially imputed. Equal physical offsets −300, −100, 0, +100 and +300 km/s at both pilot redshifts agree to 4.815 × 10⁻¹¹ km/s. The four representative local-centred profiles have centroid magnitude below 1.6 × 10⁻¹⁴ km/s.

Among jointly measured Hα spaxels, the median new-minus-V1 velocity is +1.928 km/s for 11004-12701 and −65.509 km/s for 8602-12705. V1 included 261 imputed Hα spaxels in the latter object. These changes combine the vacuum reference, retained systemic origin, channel selection and explicit measurement eligibility; they were not tuned to reproduce the old map.

## 4. Native wavelength-dependent LSF

The header inventory of all 46 containers found **40 DRP v2_7_1 cubes** with native PREDISP/DISP and **six v3_1_1 cubes** with LSFPRE/LSFPOST. This inventory reads metadata only. In both generations, native HDUs 4 and 5 are the POST and PRE LSF cubes, followed by WAVE in HDU 6. Every MEGACUBE also contains a later **DISP with 30 line planes**, which is a derived product and is excluded.

The native products contain Gaussian σλ in Å. PREDISP corresponds to LSFPRE; DISP corresponds to LSFPOST. PRE describes the LSF before detector-pixel integration and POST includes that integration. Semantic normalization does not assert identical calibration across DRP versions. The relation to resolving power is R = λ/[2√(2 ln 2) σλ], conventionally written with 2.355. [DR16 LSF documentation](https://www.sdss4.org/dr16/manga/manga-data/working-with-manga-data/), [DR17 LSF documentation](https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/), [Law et al. (2021)](https://arxiv.org/abs/2011.04675).

`read_manga_lsf()` validates the DRP version, native extension pair, dimensions and WAVE length, and checks surviving LSF WCS/unit cards against FLUX and the DRP definition. The Sandra containers omit LSF BUNIT and WCS cards: their units and pixelization semantics therefore come from the documented native DRP contract, not from nonexistent header evidence. The reader records that distinction. Both canonical arrays are in FITSio order **(x, y, wavelength)**, exactly like FLUX; Astropy uses **(wavelength, y, x)**. At 63 asymmetric spatial/wavelength samples, independent readers agree exactly, including explicit invalid-sample handling.

The canonical names are `lsf_sigma_angstrom_pre` and `lsf_sigma_angstrom_post`. Provenance retains native extension name/HDU, DRP version, WAVE sampling, dimensions, units source, medium and pixelization definition. Nonpositive/nonfinite LSF values remain missing; there is no interpolation or scalar replacement. All seven full native LSF cubes have finite positive values where valid, PRE ranges 0.902–2.473 Å and POST ranges 0.946–2.505 Å. POST² − PRE² is positive at every jointly valid sample. The output WAVE step is not assumed to equal the detector-pixel integration width.

These medians use the fixed spectral support at the native channel nearest the predicted line centre; they do not imply a single constant LSF over the spectrum. All extrema and sample counts are in [lsf_line_summary.csv](../../../results/science_baseline_freeze_v2/physical/lsf_line_summary.csv).

| Object | DRP | Line | PRE σλ [Å] | POST σλ [Å] | POST σ [km/s] | POST R |
| --- | --- | --- | --- | --- | --- | --- |
| 11004-12701 | v3_1_1 | halpha | 1.5207 | 1.5716 | 70.48 | 1806 |
| 11004-12701 | v3_1_1 | oiii5007 | 1.2362 | 1.2759 | 75.00 | 1697 |
| 11014-3704 | v3_1_1 | halpha | 1.4633 | 1.5157 | 67.56 | 1884 |
| 11014-3704 | v3_1_1 | oiii5007 | 1.1929 | 1.2323 | 72.00 | 1768 |
| 8602-12705 | v2_7_1 | halpha | 1.4658 | 1.5231 | 67.41 | 1889 |
| 8602-12705 | v2_7_1 | oiii5007 | 1.1375 | 1.1825 | 68.60 | 1856 |
| 8932-3701 | v2_7_1 | halpha | 1.4981 | 1.5548 | 69.24 | 1839 |
| 8932-3701 | v2_7_1 | oiii5007 | 1.2735 | 1.3140 | 76.70 | 1660 |
| 9869-9102 | v2_7_1 | halpha | 1.5220 | 1.5779 | 70.11 | 1816 |
| 9869-9102 | v2_7_1 | oiii5007 | 1.1501 | 1.1948 | 69.58 | 1830 |
| 8941-12704 | v2_7_1 | halpha | 1.4846 | 1.5413 | 69.42 | 1834 |
| 8941-12704 | v2_7_1 | oiii5007 | 1.1546 | 1.1986 | 70.76 | 1799 |
| 9039-9101 | v2_7_1 | halpha | 1.5928 | 1.6488 | 66.91 | 1903 |
| 9039-9101 | v2_7_1 | oiii5007 | 1.2335 | 1.2762 | 67.88 | 1875 |

`select_manga_lsf()` requires both the fitting method and its pixel-sampling treatment. A direct Gaussian evaluated at pixel centres selects POST. A template model that explicitly integrates pixels selects PRE; a point-sampled template selects POST. Template convolution by itself does not prove that pixel integration is included. These choices follow the SDSS method-dependent guidance and are tested for both native extension conventions. A controlled Gaussian-plus-pixel-integral test recovers its known intrinsic width with PRE and demonstrates the bias from counting pixel integration twice.

The native CAPIVARA moment and SpectroPath calculations fit neither of those instrumental models. They retain both native LSF arrays over each extracted line window and explicitly identify their widths as **observed**, with `correction_applied=FALSE`. A scalar FWHM is never substituted.

![Native LSF generations](../../../results/science_baseline_freeze_v2/physical/lsf_generations.png)

*Native PRE and POST σλ for one retained spatial sample in each of the two smoke-test galaxies. Solid blue is PRE and dashed orange POST. The wavelength dependence and localized structure are retained as read; the curves are not fitted or smoothed. The table above summarizes the fixed support, whereas these curves show individual spatial samples.*

The separate capivaraPPXF audit found a legacy joint SPS+gas implementation with one scalar FWHM shared by all regions. Its pPXF gas templates include pixel integration (`pixel=True`), implying PRE for that component. The stellar-template medium, template sampling and region-dependent LSF mixture require a separate resolution-matching calculation. The backend now preserves native LSF discovery metadata and **fails explicitly** on such input. It also removes the old fallback that could silently create default-air gas templates after a `vacuum=True` error. This is an enforced limitation, not a claim that native-LSF population fitting has been repaired.

The observed descriptors therefore remain unsuitable for unqualified intrinsic-dispersion inference. Line blending, continuum adequacy, IVAR/MASK science cuts and uncertainty calibration remain separate requirements. Integrating native Fλ along velocity has units of flux density times km/s; it is not automatically calibrated line flux.

## Validation, source and installation

CAPIVARA passes **32 test cases / 471 expectations**, with zero failures, warnings or skips. capivaraPPXF passes **52 expectations in 11 cases**, with zero failures and one intentional skip for the missing-pPXF branch because pPXF is installed. Both final source archives pass `R CMD check --no-manual`: **zero errors, warnings and notes**; CAPIVARA’s vignettes were rebuilt.

The regressions cover both spectral APIs and widely separated redshifts, exact channel indices, missing redshift, ambiguous metadata, no double correction, variance/SNR alignment, provenance propagation, line-medium mismatch, explicit conversion/inversion, systemic and local profiles, translation-invariant signatures, both DRP extension pairs, wrong units/dimensions, exclusion of derived DISP, invalid LSF samples and method-dependent pixel integration. An additional pre-existing small-cube SNR check was corrected so an unused default Ncomp cannot reject target-SNR selection; the fixed-18 comparisons are unaffected.

Normal and isolated installations match all **200 CAPIVARA functions, 21 exports and 12 bundled scripts** and all **18 capivaraPPXF functions, 13 exports and its Python helper**. SpectroPath’s installed 40 functions and 18 exports match its unchanged source.

| Component | Version | Validated source commit |
| --- | --- | --- |
| CAPIVARA | 0.4.2.9000 | `4e6263af3fd0bfc556d4a2273c76e4ce535f7829` |
| capivaraPPXF | 0.0.2.9000 | `cad2e5259f98f57570d39d0d515a0d9e35f1823c` |
| SpectroPath | 0.1.0 | `e1695030c50c09adf5eebdde56336c2342493c25` |
| Historical CAPIVARA implementation | 0.4.0.9000 | `2cdce2f54526ecc9b04eb3aa7f2a54b51f1e7100` |
| Observed/rest partition run | 0.4.1.9000 | `a243fdfe551a66c32e90ee39e04ce113918fb9aa` |

The seven partition comparisons were computed at the recorded 0.4.1 commit. The final 0.4.2 audit re-read all seven unchanged input SHA256s and reproduced their complete wavelength-selection provenance exactly. The clustering implementation, scaling, masks and selected flux samples are unchanged; the new registry and LSF metadata do not enter full-spectrum clustering. The four public kinematic runs and seven LSF audits use the final 0.4.2 runtime. The later bundled-tutorial correction changes no function or native smoke-workflow calculation.

The historical source/report head was `670a61857e74a840822ec026c659ba304de19c4d`; the retained historical capivaraPPXF pin is `9bfdea212e1704e80ca74dd715801df893e8ca07`. R 4.5.2 on arm64 macOS was used. The manifest records full package sessions, source and artifact SHA256s, normal/isolated install checks, all four physical contracts, the comparison metrics and readiness decisions.

## Evidence and scope

Evidence root: `results/science_baseline_freeze_v2/`.

- `validated/comparisons/` and `validated/figures/`: seven paired partitions, region sizes, spectra, overlap assignments, fragmentation and conservation.
- `physical/exact_windows.csv`, `line_coordinate_checks.csv`, `line_registry*.csv`: physical intervals, laboratory references and round trips.
- `lsf_header_inventory.json`, `physical/lsf_contracts.rds`, `lsf_line_summary.csv`, `orientation_verification.json`: native DRP semantics, sample values and independent orientation check.
- `physical/kinematic_smoke/`: four complete public-API runs, representative systemic/local profiles and both native LSF arrays over each line.
- `tests_physical.rds`, `backend_physical_tests.rds`, `physical_source_install_checks.txt`, `build_physical/`: tests, source/install comparisons and clean package checks.
- `freeze_manifest_v2.rds`: exact source pins, SHA256 evidence ledger and the physical readiness gates.

`tools/validate_wavelength_frames.R`, `validate_systemic_kinematics.R`, `validate_physical_coordinates.R`, `verify_lsf_orientation.py`, and `finalize_wavelength_freeze.R` reproduce the scoped checks. No 46-object production run, hierarchical segmentation, website/README narrative edit or push was performed.
