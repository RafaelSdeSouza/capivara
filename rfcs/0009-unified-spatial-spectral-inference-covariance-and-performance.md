# RFC 0009 — Unified spatial–spectral inference: covariance, Fisher–Rao, kinematics and speed

**Status:** proposed scientific contract for Codex research implementation. No public interface changes authorized.  
**Baseline:** current main; related [RFC 0008](0008-decision-theoretic-redshift-information-audit.md).  
**Deliverable:** isolated experimental module, exact covariance tests and draft PR; do not merge or release.

## Central scientific decision

Analyze each galaxy in a **declared systemic rest-frame wavelength coordinate**, without making its systemic redshift a clustering feature. Internal LOS velocities remain nonzero relative to this frame. For **kinematics**, those line shifts are parameters of interest; for **other spectral representations**, they may be nuisances that must be explicitly modeled instead of silently removed. The redshift also remains important in converting measured data into rest-frame coordinates, observed-passband coverage, the instrument PSF/LSF, and kpc/arcsec geometry across objects.

Build **one shared observation model, one covariance engine, one likelihood/fit interface, one region geometry and one information audit**. Retain distinct task-specific spectral representations and distinct outputs for:
- spectral-continuum/absorption,
- emission-line line-strength/ratio/profile,
- stellar/gas kinematic fields.

### Existing source to preserve

`segment_emission_lines()` is *already implemented* in `R/segment_emission_lines.R` (line windows, moments, profiles and derivatives); `segment_kinematics()` and a native kinematic runner already exist. `capivara_representation()`, `R/observed_entry_hierarchy.R`, `R/science_wavelength_frames.R`, `R/spatial_geometry.R` supply scientific contracts and exact region polygons. `prepare_segment_spectra()` and `inst/python/capivara_resolution.py` provide experimental LSF handling and banded spectral covariance under **independent native spaxels**. `tools/audit_spatial_covariance.py` contains a four-band MaNGA covariance diagnostic but not a complete fitting correction. `docs/SCIENCE_BASELINE_FREEZE_V3.md` explicitly says **VARIANCE_CONTRACT FAILED** for unpropagated native spatial covariance, and the stellar population pilot remains **NOT_READY**. These limits must be preserved and raised as hard validation gates. `docs/CAPIVARA_2_COHERENCE_AUDIT.md` also identifies a sampled-vs-unsampled scaling divergence in the emission feature path.

## Joint forward model

For each intrinsic spaxel `p`, let its physical spectrum be
`F_p(lambda)=F_star,p(lambda;theta_star,p)+F_gas,p(lambda;theta_gas,p)`. Its continuum/absorption, stellar and gas velocities, gas line fluxes and (possibly multicomponent) dispersions are task-specific components of the same likelihood. Observed means are `mu(theta)=A_PSF,LSF,sampling[F(theta)]`. Signal response `A` and noise covariance `Sigma` must remain *separate*. A region formed from spectra of different velocities has a mixture profile; a single-Gaussian fit to the sum is not in general a measurement of its intrinsic dispersion.

For Gaussian, parameter-independent `Sigma`:
`-2 log L=(y-mu)^T Sigma^{-1}(y-mu)+logdet Sigma+constant`.
Optional spatial regularizers act on fitted parameter fields (graph-GMRF/Laplacian, edge-aware variants), explicitly as **priors** rather than observational evidence. Control their effect on sharp physical boundaries and identifiability. Use pPXF as a validated backend where scientifically compatible; do not rewrite it by default. Separate line-feature measurement and dynamical model hypotheses.

## Covariance must be correct before fitting

For each fixed regional sum/average `B`, propagate **data, response and covariance with the same operator**:
- `mu_B=B mu`, `J_B=B J`, and **`Sigma_B=B Sigma B^T`** including covariance *between different regions*.
- Native diagonal variance/IVAR; masks and support must remain missing, never silently imputed as zero.
- Real MaNGA DR17 sparse spatial matrices `GCORREL`, `RCORREL`, `ICORREL`, `ZCORREL` carry genuine spatial noise correlations; respect native index order, symmetry, band metadata and positive-definite checks.
- LSF / rest-grid resampling creates **spectral** correlations, which must be propagated by a declared linear operator `R Sigma R^T`.
- Joint spatial–spectral covariance may be nonseparable. If only four spatial correlation anchors and a spectral kernel are available, **do not mislabel an approximation as the exact full covariance**. Carry explicit approximation scope, error tests and provenance. Never interpret PSF signal mixing as a substitute for noise covariance.
- Implement matrix-free covariance operators and bounded/sparse alternatives instead of dense covariance of dimension `(n_spaxel*n_wave)^2`.

Official MaNGA reference: https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/ ; MaNGA DAP covariance guidance: https://sdss-mangadap.readthedocs.io/en/latest/spatialcovariance.html .

## Fisher–Rao / decision-aware audit

For a **fixed** partition and Gaussian `Sigma` independent of parameters:
`I_full = J^T Sigma^{-1} J`, `I_B = J^T B^T (B Sigma B^T)^{-1} B J`, and `I_full - I_B >= 0` (PSD, up to numerical tolerance). Compute nuisance-adjusted efficient information with a correctly scoped Schur complement. For ill-conditioned or singular cases, return an identifiability warning rather than a spurious interval. If noise covariance depends on the model parameters, the usual score Fisher contains an additional trace term and must not be approximated by the mean-Jacobian term without an explicit declaration.

This **first stage is a read-only diagnostic of existing region partitions**, not a new clustering algorithm. Treat labels selected from the same noisy observations as data-dependent; evaluate selection effects with independent simulations or data splitting. Distinguish an unsupervised learned spectral representation from parameter-specific Fisher–Rao geometry.

## Optimize without changing science

- Time and profile each stage: feature extraction, covariance, spectral summaries, backend fitting, Fisher solves and R/Python crossings.
- Use native sparse CSR/CSC blocks for spatial covariance; banded/sparse spectral operators; cached sparse binning `B`; batched `Sigma^{-1} J` solves (few RHS) without global inverse; local spectral windows for line fits; FFT or direct convolution as justified by the actual PSF/LSF boundary conditions.
- Vectorize or batch regional fitting/preparation to avoid per-region `reticulate` crossings when performance tests justify refactoring.
- Use the existing C++ graph code and exact polygon indexing; do not introduce dense pairwise spaxel distances for large IFUs.
- Profile CPU/Rcpp/NumPy/SciPy first. Consider optional JAX/torch/Metal only with measured end-to-end benefit and R/Python numerical parity. Report wall time, peak memory and accuracy, not an unverified speedup.
- For every approximation, specify when it is valid, error tolerance and the independent dense reference test.

## Mandatory acceptance checks

1. Rest-coordinate equivalence, consistent flux-density Jacobian and covariance transformation, preservation of relative local LOS velocities.
2. Exact two-spaxel Fisher benchmark: common-mode information retained by binning, contrast information lost; independent and spatial/spectral correlated noise.
3. `B Sigma B^T` equals small dense reference, *including inter-bin covariance*; masks not treated as zero observations.
4. Correct reconstruction of real MaNGA covariance sparse triangles/band metadata with no double-counted off-diagonal entries.
5. Fisher data-processing PSD and monotonicity under refinement; nuisance degeneracy returns explicitly nonidentifiable.
6. No unsupported single-Gaussian interpretation of merged multi-velocity profiles.
7. Sampled/unsampled emission-line features: documented and tested normalisation and missingness behavior; no unintended change in historical outputs.
8. Synthetic no-gradient, two-zone, weak-line, bad-support and low-S/N tests; compare methods on held-out cases, not only post-selection conditional estimates.
9. Stable existing tests (especially `tests/testthat/test-semantic-representations.R`, `test-spatial-foundation.R`, `test-line-medium-lsf.R`, `test-science-products.R`) must keep passing. No changes to public functions in stage 1.
10. Reproducible benchmarks on a tiny dense exact case and MaNGA-scale cubes. A computational shortcut is only acceptable if it meets a declared scientific accuracy tolerance.

## Codex execution order

**Phase 0:** inspect current branch and source; audit all current contracts/known failures.  
**Phase 1:** implement experimental covariance operator and exact numerical tests.  
**Phase 2:** read-only Fisher/decision-information audit for fixed existing CAPIVARA partitions.  
**Phase 3:** independently validated joint star+gas fitting demonstration with optional spatial regularization; compare with pPXF/MaNGA DAP, and calibrate uncertainty.  
**Phase 4:** consider task-aware stopping or new hierarchy *only after* Phases 1–3 pass.

Work on a new `codex/unified-inference-covariance-audit` worktree branched from current main. Read this RFC and RFC 0008 from the separate documentation design branch; do not merge the documentation branch. Produce a draft PR documenting exact tests, hardware, speed/memory measurements, scientific limitations and PASS/FAIL gates. **Never merge, publish or change the released website without separate approval.**
