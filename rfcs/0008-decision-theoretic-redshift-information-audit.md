# RFC 0008 — Decision-aware partition audit under changing observation models

**Status:** research proposal with a *synthetic proof-of-concept*, not validated CAPIVARA code or a novel theorem.  
**Date:** 2026-11-20.  
**Architectural position:** read-only audit of existing R partition results first; exploratory method construction only after scientific validation gates.  
**Existing foundations:** [architecture](0001-platform-architecture.md), [scientific contracts](0002-scientific-api-specification.md), [visuals](0007-v7-colour-and-mobile-audit.md).

## Scientific problem

Ordinary spectral similarity does not determine which spatial differences matter for a particular astronomical estimand. Spatially summing two spectra can improve a global signal-to-noise ratio and simultaneously eliminate sensitivity to a physically interesting contrast. The choice of redshift, rest-frame coverage, spatial PSF, spectral LSF, sampling and instrumental noise changes which physical quantities are distinguishable.

**Candidate first target:** spatially varying emission-line ratios or line-profile contrasts, with a physically specified local line model and nuisance parameters. Do **not** translate any [N II]/H-alpha contrast directly into metallicity without separate excitation/calibration assumptions.

## Existing CAPIVARA components that must be preserved

- The **legacy all-pairs** L1/Ward.D2 hierarchy, the **sparse Euclidean graph Ward** hierarchy, and the **explicit observed-entry** spatial hierarchy are distinct; do not rebrand any one as another.
- \`capivara_representation\`, wavelength-frame/air-vacuum/LSF contracts, measured voxel validity and independent support are established source assets.
- \`regions_sf()\` provides exact region footprints, adjacency, holes and disconnected components; these should be the common visualization structure.
- The existing \`POST_V4_METHODOLOGICAL_ROADMAP.md\` explicitly says a fixed \`Ncomp\` has no intrinsic scientific privilege and asks for an evidence-based stopping rule. This proposal extends that programme by declaring the physical estimand and the information lost in each merge.

## 1. Observation-specific information

For a physical field parameter \`theta\`, nuisance parameters \`eta\`, intrinsic source \`F(theta,eta)\` and observation operator \`A_z\`:

\`\`\`math
Y_z = A_z[F(theta,eta)] + epsilon_z, \quad epsilon_z \sim N(0,Sigma_z).
\`\`\`

The mean vector \`mu_z\` has Jacobian \`J_z = d mu_z / d(theta,eta)\`. In the finite Gaussian model with parameter-independent covariance:

\`\`\`math
I_full(z)=J_z^T Sigma_z^{-1} J_z.
\`\`\`

For a **fixed** regional-summing operator \`B_Pi\`:

\`\`\`math
I_Pi(z)=J_z^T B_Pi^T(B_Pi Sigma_z B_Pi^T)^{dagger}B_Pi J_z.
\`\`\`

Nuisance-adjusted efficient information for a scalar physical target follows from the Schur complement of the nuisance block wherever estimable. The corresponding local Gaussian posterior-variance or decision-loss audit is **model-specific**: it is not generally equal to a Bayes risk for arbitrary scientific loss functions.

**Classical data-processing statement:** \`I_full-I_Pi\` is positive semidefinite for fixed linear compression and a fixed Gaussian model. One proof whitens \`J_z\` and observes that compression projects the whitened score onto a subspace. This is *established statistics* (closely related to MOPED); it is not a theorem original to CAPIVARA.

## 2. Ward identity in physical sensitivity space

In the independent, homoscedastic, full-coverage Gaussian model, define for region \`R\` the spatial mean of the **vector-valued Jacobian** over its pixels \`bar J_R(lambda,parameter)\`. The information lost by replacing pixels with regional sums is the within-region sum of squared Jacobian deviations, divided by noise variance.

For merging existing regions \`A,B\`, the **increase in loss of the full Fisher matrix** is

\`\`\`math
Delta I_loss(A,B)= [n_A n_B/(n_A+n_B)] sum_lambda
  (bar J_A(lambda)-bar J_B(lambda))^T
  (bar J_A(lambda)-bar J_B(lambda)) / sigma^2.
\`\`\`

This has the same algebraic form as Ward's merge criterion, but on physical-model *response vectors* rather than measured spectra. It gives a direct testable mechanism for a **task-aware spatial hierarchy**.

**Cautions:** this is a classical variance-decomposition identity, not a new theorem; greedy graph Ward is not a globally optimal partition. A scalar nuisance-adjusted Schur complement is not equal to a fixed trace of the full-Fisher loss; a specific decision functional must be declared. Correlated noise, missing channels, changing PSF/LSF, parameter-dependent covariance, learned response models, and selection of labels on the same noisy cube require nontrivial generalizations.

## 3. Re-observation across redshift

A genuine cross-z experiment must begin with the **same intrinsic physical scene**, then apply observation operators \`A_z\` with explicitly modeled bandpass, cosmological angular scale, luminosity/surface-brightness effects, spatial sampling, PSF, LSF, noise, throughput and selection. It should compare:

1. information lost to **instrumental observation** (intrinsic scene -> measured cube);
2. additional information lost to **regional compression** (measured cube -> partition sums);
3. uncertainty/decision loss due to **adaptive partition selection**.

When the relevant spectral feature leaves the passband, relative information-retention ratios may be undefined: report “not identifiable from these observations” instead of dividing by zero or inventing a confident estimate. For uncertain redshift, integrate over an explicitly specified redshift distribution or report sensitivity bounds—not simply a point-z calculation.

## 4. Executed synthetic proof of concept

An independent Python **research script**, not a port of R CAPIVARA and not part of the stable package, evaluated a 22 x 22 spaxel field with 826 observed-frame channels (4700–8000 A), Gaussian H-beta/[O III]/H-alpha/[N II] profiles, a spatial [N II] contrast target, two nuisance response directions (global [N II] amplitude and continuum slope), and a seeing Gaussian whose physical width changes with a fixed 2.5-arcsec angular FWHM in flat LCDM. The noise level per pixel/channel was held fixed to isolate spatial blur and bandpass effects; source luminosity dimming, throughput, realistically sampled pixels, wavelength-dependent PSF and observational covariance were **not** modeled. Labels for data-dependent methods were fitted on an independent simulated pilot, then evaluated through the known test-model Jacobian.

**Results at z=0.08, K=4** (nuisance-adjusted fraction of *full measured-cube Fisher information*):
- historical median-centred L1/Ward-style SciPy **analogue**: **1.50%**;
- four-neighbor **complete-data spatial Euclidean Ward surrogate**, not the V4 Rcpp implementation: **8.16%**;
- fixed square grid: **83.13%**;
- four known-orientation target-aligned stripes: **83.08%**;
- spatial Ward applied to the **true Jacobian response field** (oracle): **91.65%**.

The L1 and spatial-Ward results are poor in this *deliberately constructed mismatch*, not estimates of CAPIVARA's real-data performance. The derivative-based method is an oracle requiring forward-model knowledge; it is not deployable without parameter uncertainty, selection correction and strong tests. The sign and magnitude of the advantage are *not* established across source populations.

At \`z=0.215\` the toy [N II] sensitivity remained measurable, whereas near \`z=0.22\` the shifted line had passed beyond the deliberately limited 8000-A observed-frame bandpass; the full-data efficient information in this toy model fell to numerical zero. No partition can restore information absent from the measured data.

Six numerical unit checks passed: split-gain identity, information monotonicity under refinement, spatial connectivity, legacy analogue number of groups, redshift information disappearance and finite-grid information loss. These tests are **not** an observational validation or an R-parity certification.

The independent scientific plate, source script, fixtures and executable tests were generated as conversation artifacts for review. This RFC does **not** claim the source is committed in the CAPIVARA repository or that CI has run it.

## 5. Prior art and strongest comparisons

- [Heavens, Jimenez & Lahav, MOPED, MNRAS 317, 965 (2000)](https://doi.org/10.1046/j.1365-8711.2000.03692.x): model-derivative-weighted, locally Fisher-preserving spectral data compression.
- [Cappellari, PowerBin, MNRAS 544, 1432 (2025)](https://doi.org/10.1093/mnras/staf1726): adaptive spatial binning based on centroidal power diagrams, including capacity constraints that may already be tied to physically relevant uncertainty; strong practical comparison required.
- [RealSim-IFS (MNRAS 514, 2821, 2022)](https://academic.oup.com/mnras/article/514/2/2821/6603846): forward simulation of instrument/seeing/redshift effects on spatially resolved galaxy spectroscopy.
- [CAPIVARA, MNRAS (2025)](https://doi.org/10.1093/mnras/staf688): existing spectral segmentation and its declared variants.

**This is a real mathematical and software research problem, but not yet demonstrated literature-level novelty.** The most valuable first release could be a diagnostic, not a new segmentation method.

## 6. Architectural next step

**Stage A: read-only audit.** Implement an independent \`audit_partition_information\` experiment that accepts existing CAPIVARA region IDs, a specified estimand, a named physical forward model or derivative/score arrays, units, noise/LSF/PSF assumptions and redshift(s), and produces:
- full vs partition Fisher matrices and nuisance-adjusted target information;
- information-retention curves at existing hierarchy cuts, with “unobservable” states;
- region-resolved loss maps, exact polygons and linked spectral diagnostics;
- frozen-label/data-splitting provenance and warnings about model mismatch.

**Stage B: selection-adjusted science.** Evaluate coverage, posterior/decision risk and model sensitivity with held-out synthetic truth, fit nuisance uncertainty, realistic MaNGA/MUSE noise/PSF/LSF and actual CAPIVARA hierarchy outputs. Avoid treating selected labels as fixed if they were estimated from the same scientific cube.

**Stage C: only if Stage A/B survive, experimental score-space hierarchy.** Use physically validated Jacobian-response features and a clear scalar/matrix decision criterion. Benchmark against CAPIVARA V4, PowerBin, task-specific line fits and MOPED-like score compression, at matched spatial resolution and computational budget.

**Novel mathematical paper gate:** nontrivial *selection-aware* or *instrument/uncertainty-robust* guarantee, certified approximation or scientifically new identifiable functional, beyond the classical Fisher-data-processing and Ward variance identities. Failure means “useful audit tool”, not a new theorem or operator.

**Current decision:** keep this work out of the stable segmentation API. First run the read-only information audit on frozen native R region arrays and a scientifically defensible line-response model.
