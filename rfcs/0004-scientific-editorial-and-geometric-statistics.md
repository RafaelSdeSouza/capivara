# RFC 0004 — Scientific editorial design and geometric statistics

**Status:** proposed design and research programme. **No website has been deployed and no new mathematical operator has been implemented.**  
**Date:** 2026-11-20  
**Companions:** [Architecture](0001-platform-architecture.md), [Scientific API](0002-scientific-api-specification.md), [Validation and website](0003-validation-documentation-release.md).

## 1. Editorial direction: as precise as the best scientific packages

**Recommended direction:** *Scientific Editorial* — predominantly light, exceptionally clear scientific figures, minimal ornament, readable typography, understated accent colours and compact mathematics embedded at the point where it explains a result. A darker *Spectral Atlas* treatment may be appropriate for an opening page or animated physical-field demonstration, but must never displace legible light reference pages.

**Reference patterns, not visual copies:**
- [Astropy docs](https://docs.astropy.org/en/stable/index.html): clear role of getting started, concepts, guides and API.
- [Photutils](https://docs.astropy.org/projects/photutils/en/stable/index.html): task-based scientific reference, meaningful astronomical examples and citation guidance.
- [SunPy](https://docs.sunpy.org/en/stable/): tutorial, how-to, example gallery, topic guides and references as *different user journeys*.
- [SunPy example authoring](https://docs.sunpy.org/en/latest/dev_guide/contents/example_gallery.html): each example motivated by an astrophysical question, linked to scientific papers and illustrated with a reproducible figure.
- [Matplotlib gallery](https://matplotlib.org/stable/gallery/index.html) / [scikit-learn visualization API](https://scikit-learn.org/stable/visualizations.html): composable figures that render existing calculated results without secretly rerunning inference.

### Writing voice

- Scientific, concise, natural English as in a strong methods paper or well-maintained astronomical software manual. Avoid promotional abstractions: “unlock”, “seamless”, “cutting-edge”, “revolutionary”, “powerful toolkit”, “new paradigm”, “AI-driven insights”, “physically meaningful” without defining the estimand.
- Speak in **specific astronomical quantities**: spaxels, observed/rest wavelengths, measured flux, finite channel coverage, variance, line-of-sight velocities, PSF/LSF, line profiles, detected support, regional sums.
- Headings should name a question or outcome, not advertise the mechanism: “Which spaxels have similar spectra?”, “Does the result depend on the wavelength interval?”, “When does a symmetry survive observation?”
- The first paragraph states the purpose, the quantity being compared and a limitation. Follow with a small equation or appropriate plot, not a list of internal function names.
- Every technical term used in a headline is defined in the page; every equation states its mathematical domain and the units/frame of the relevant symbols.
- Do not use institutional-sounding slogans or generic marketing copy in scientific articles. Readability comes from explanation, examples and good section sequencing, not lowering mathematical standards.
- Distinguish **software correctness**, **mathematical guarantee**, **instrument-scoped scientific validation** and **observational conclusion** in the same sentence where necessary.

**Reference rewrite example**

Weak: “Cluster a 2D representation of a data cube. Processes, scales and computes distances.”

Preferred:

> **Segment an IFU cube by spectral similarity.** Groups eligible spaxels according to a declared spectral representation and a specified hierarchical merge criterion. The resulting regions summarize measured spectral similarities; they are not automatically distinct stellar populations, gas phases or connected structures. The historical `segment()` path uses median-centred Manhattan distances followed by Ward.D2 updates. Use the explicit observed-entry representation when missing-channel semantics and validated eligibility rules are essential.

**Reference rewrite for a math page**

Weak: “A transformative equivariant operator that preserves the geometry of spectra.”

Preferred:

> **When does a symmetry survive the instrument?** A Doppler shift has a well-defined action on an intrinsic spectrum. A broadband filter discards spectral information, so an exact corresponding action on the measured flux need not exist. We use the measurement null space to test this obstruction before constructing a model.

**Reference rewrite for a visual caption**

> *Synthetic demonstration:* the left panel shows a simulated spatially resolved source; the centre panel is a categorical partition of selected pixels; the right panel shows representative region spectra. The regions encode an algorithmic grouping, not physically identified components. Real MaNGA results should replace this fixture for the published gallery.

### Editorial safeguards

- Always label simulation, real observation, fitted model, residual, validation benchmark and schematic separately.
- Never claim that all current R functions have identical scientific validation; scope evidence to each method/instrument.
- Do not merge “R + Python” in installation guidance until Python core has its own released and validated scope. Proposed navigation may state “Python in development”.
- Preserve current package logo, copyrights, author and publication attribution unless expressly approved otherwise. Preview marks are temporary placeholders.
- Existing documentation URLs and cite instructions must remain discoverable and, where possible, backward-compatible.

## 2. Visual language

**A. Scientific Editorial (preferred for public documentation).** White canvas, refined Inter-like sans, very restrained teal/blue accents, large multi-panel maps/line plots, fine lines, generous whitespace and publication-style captions. Category colors are consistent across map and spectra. Equations are typographically precise but never visually larger than the primary data figure. The page should be understandable in full-page view and at mobile width.

**B. Spectral Atlas (for science storytelling).** Dark navy backgrounds and luminous—but correctly scaled—physical-field imagery for the homepage and special demonstrators. Switch to a light surface for long mathematics and API documentation. Physical interpretation must come from labels and explicit diagnostics, not cosmic aesthetics or arbitrary contours.

**C. Methods Notebook (mathematics and reproducibility).** Persistent conceptual navigation, numbered assumptions, compact propositions, plots immediately under the equations, coded counterexamples, links to reference papers and validation/QC notes in a separate panel. The operator's **known** mathematics and the research hypothesis must have visibly different labels.

**Visualization rules:** figure-first; 3:2 and 16:9 display figures, with journal/vector export; physically correct axes and WCS conventions; spectral units and line identity; categorical palettes for labels (no numerical gradient), continuous scales for actual physical scalars, blue/orange or accessible contrasts for paired transformations, and no decorative footer on figures.

The screenshots produced for this design discussion are **local mockups** with synthetic data. They are neither an installation preview of the current package nor a proposed scientific result. Replace simulation with verified real-data figures only after validation and provenance checks.

## 3. Where the mathematics enters the codebase

Place the new scientific reasoning at **three contracts**, not in an oversized `segment()` argument list.

### Layer I — What is the representation?

Intrinsic spectral-spatial field `F(r,u)`, where `u = log λ`. For a declared group `G` with action `U_g` on input fields and `ρ(g)` on output features, a representation `T` is equivariant if

```math
T(U_g F)=\rho(g)\,T(F).
```

**Scientifically justified choices:**
- Constant Doppler/redshift shift, `u→u+\delta` with `\delta=log(1+z)`, corresponds to the additive translation group in log-wavelength. Account for Jacobians/flux-density convention rather than translating arbitrary arrays as though they were the same physical quantity.
- Idealized rotations/translations of a planar IFU field admit `SE(2)` actions, but detector PSF anisotropy, footprint, sky coordinates and a galaxy's actual physical anisotropy mean **invariance is not necessarily desired**.
- Positive amplitude scaling defines a multiplicative action, with conditional `F/A(F)` invariance where `A(F)>0`. If absolute brightness is the target, retain amplitude as a separate observable.
- Axial directions (filaments, projected magnetic orientations) inhabit a projective space (e.g. `RP¹` on the sky), not an oriented-vector space. Spin-weighted fields require their own representation laws; do not infer that every spectral cube is a spin field.

**Not all transformations are groups:** wavelength filtering, PSF/LSF smoothing, finite pixel sampling and other non-invertible operators do not form a conventional invertible group action. Fixed-law attenuation can have a compositional semigroup interpretation under suitable assumptions. The correct composition law matters more than the label “equivariant”.

### Layer II — What survives a measurement?

Let `A:V→W` be a linear idealized observation operator (e.g., bandpass integration or PSF/LSF convolution on finite discretizations). A linear symmetry action `U_g:V→V` descends to `im(A)` exactly if an induced linear map `V_g` exists with

```math
A U_g = V_g A.
```

**Classical exact criterion:**

```math
\exists V_g\ \text{on im}(A)\text{ with } A U_g=V_g A
\iff U_g(\ker A)\subseteq\ker A.
```

**Elementary proof:** if `Ax=0`, factorization forces `AU_gx=0`. Conversely define `V_g(Ax)=AU_gx`; this is well-defined precisely because any two representatives of the same measurement differ by a null-space element whose transform also lies in the null space. If `G` includes inverses, the full family of such maps satisfies the group laws on `im(A)` when the condition holds for all `g∈G`.

This is **not** a new theorem. It identifies where naive spectral equivariance is impossible after measurement. In particular, two intrinsic spectra may give the same single filter flux, but their redshifted versions different fluxes. A model that demands an exact induced redshift action on that one measured scalar is asserting unavailable information.

When `A` and `U_g` fail to intertwine, quantify an **observable, norm-defined defect** under a physically specified instrument and restricted source class; avoid unqualified global operator norms when the spaces are not bounded/comparable. For spatially varying velocity `δ(r)` and PSF `B`, `[B,T_δ]` is generally nonzero; beam smearing is established physics, not a discovery.

### Layer III — Which statistical estimands remain determined?

For a target physical functional `Q(F)` and compatible model class `M`, exact identification from a noiseless measurement `y=AF` means `Q` is constant on each nonempty observational fiber `{F∈M:AF=y}`. Under uncertainty, define an admissible set/range with explicit noise and physical assumptions. That definition is classical; the research goal is a **tractable theorem or operator** that provides useful sharp or certified approximate bounds for scientifically relevant spectro-spatial covariance/statistics.

**No new methods should be claimed until:** (i) the strongest relevant prior art is identified, (ii) the structural obstruction is demonstrated, (iii) a nontrivial construction/guarantee is derived, and (iv) independently held-out tests establish an astrophysical advantage.

### Research API boundary (proposal, not implemented)

Stable:
- `representation`: explicit mathematical meaning/invariances, support/coverage and retained amplitudes;
- `segmentation`: stable objective and validated regions;
- `regions`: scientifically honest spectra, covariance/variance and reconstruction;
- `visualization`: views of measured/reconstructed data.

Experimental:
- `experimental.geometry`: transformation declarations, action/response compatibility and counterexamples;
- `experimental.spectrospatial`: cross-line transport/covariance estimands, instrumental commutator/observability diagnostics;
- `experimental.identification`: conditional identified sets and sharpness/approximation certification.

No exported methods named `equivariant_segment`, `gauge_spectral_operator` or `information_preserving_transform` until precise operational semantics and meaningful empirical baselines exist.

## 4. First falsification gate before coding

1. **Filter counterexample:** build two nonnegative spectra with identical flux through a given filter but different flux after a known redshift, demonstrating no induced exact scalar action. Then determine if multiple physically feasible filters or SED constraints yield informative bounds.
2. **IFU noncommutativity:** simulate two lines with matched intrinsic kinematics but different emissivity distributions; convolve with known spatial PSF and wavelength-dependent LSF. Compare measured line-profile correlation under shift-before-PSF versus PSF-before-shift. Report actual physical test quantities and uncertainty.
3. **Baselines:** compare against spectral correlation function, standard velocity alignment, instrument-forward modelling, SPENDER-like physically transformed latent approaches and classical set-valued inverse-problem constructions. An attractive figure, by itself, is not a novel statistical operator.

**Decision:** Prioritize an *observability-compatibility* mathematical investigation while keeping standard CAPIVARA segmentation stable. Borrow from gauge and group-equivariant theory with precise attribution, rather than presenting group equivariance as new.

## 5. Scientific source anchors

- [Cohen, Geiger & Weiler (2019), general theory of equivariant CNNs](https://proceedings.neurips.cc/paper/2019/hash/b9cfe8b6042cf759dc4c0cccb27a6737-Abstract.html) — classification of group-equivariant kernels on homogeneous spaces.
- [Cohen et al. (2019), Gauge Equivariant CNNs](https://proceedings.mlr.press/v97/cohen19d.html) — intrinsic local frame transformations.
- [Melchior et al., Autoencoding Galaxy Spectra I](https://arxiv.org/abs/2211.07890) — redshift/resolution and observational transformations already modelled.
- [Liang et al., Autoencoding Galaxy Spectra II](https://arxiv.org/abs/2302.02496) — redshift invariance already investigated; do not recreate it.
- [CAPIVARA original MNRAS paper](https://doi.org/10.1093/mnras/staf688) — published scientific baseline.
- [Current CAPIVARA semantic-representations vignette](../vignettes/semantic-representations.Rmd) — existing conditional scale invariance and explicit validation domain.

**Important:** The mathematical obstruction is useful design guidance; *a new observationally compatible statistical operator is not yet established*. This RFC proposes where research should sit in the platform, not that any theorem or validated method has been delivered.