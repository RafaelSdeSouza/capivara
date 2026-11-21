# RFC 0007 — V7 colour specification and screen-by-screen mobile audit

**Status:** visual-design specification and audit, **not** a code or science change.  
**Reviewed:** 20 November 2026.  
**Source inspected:** CAPIVARA v6 mobile/desktop HTML preview, real MaNGA 7443-12703 figure assets and linked-aperture example 10218-12703; current online pkgdown information.  
**External comparisons:** [Astropy](https://docs.astropy.org/en/stable/), [SunPy](https://docs.sunpy.org/en/stable/), [S2FFT](https://astro-informatics.github.io/s2fft/), [CARTA](https://cartavis.org/), [Cubeviz](https://jdaviz.readthedocs.io/en/latest/cubeviz/index.html).

## A. The scientific colour grammar

| Role | V7 rule | Reference |
| --- | --- | --- |
| Page | #FFFFFF paper, #F3F7FA secondary scientific panel | White, not dark UI by default |
| Main text | #17242E graphite; secondary text #526471 | Never low-contrast labels |
| Interface selection | #1959D1 cobalt, reserved for links, focus, selected tabs | Do not repeat selected-state blue as a continuous data scale without a legend |
| Borders | #DCE5EA, no ornamental frames inside a scientific plot | Thin separators |
| Region identity | Categorical only: initial eight swatches #1959D1 #E76532 #00A6AD #8C54CE #D39B12 #CE4678 #36946B #399ADE | **Not sufficient alone for a 45-region partition.** |
| 9–45 regions | Generate a wide-separation palette (Glasbey-like, tested for visual accessibility), then use the **actual region-adjacency graph** to reduce near-identical colours on neighbouring polygons. Add stable region IDs, topology-preserving edges and selection highlighting | A colour alone cannot make 45 categories universally distinguishable |
| Nonnegative intensity | Matplotlib magma or inferno, with specified physical scale and normalization | A signed continuum-subtracted map should not use a nonnegative-only intensity rule |
| Signed velocity | RdBu_r centred at a declared physical zero (e.g. systemic-relative km/s), with scientifically defended limits | Not interchangeable with categorical regions |
| Uncertainty/variance | cividis or a verified monotone sequential map; missing uncertainty is masked, not minimum colour | Use correct units and uncertainty type |
| Selected region | Full-opacity fill, contrast-outline (white + dark), searchable ID; unselected region de-emphasized without loss of label identity | Same colour for corresponding line/map/legend only **within the same partition** |

**Preserve the source science:** The v6 MaNGA region raster is an actual analysis figure. It must *not* be recoloured by pixel hue or interpreted as continuous intensity. Regenerate the map from the original region-ID matrix and a validated geometry object before applying new categorical colours. The original figure's colour choices are historical data visualization, not a V7 scientific colour reference.

**Plot rules:** every colorbar supplies quantity and units, every truncation and nonlinear stretch is declared, and a region-ID map uses no continuous gradient key. When 45 nominal groups make the map unreadable, show a readable selected-region state and an overview with ID lookup; never silently collapse or invent labels.

## B. Screen-by-screen V6 critique

### 1. Header and first mobile screen

**Observed:** a decorative rainbow rail, a multi-colour square logomark, a large blue Documentation button, a large Georgia headline and smaller caption. These compete with the data.  
**Relative to:** Astropy — the main navigation quickly identifies Getting Started, Guide and reference; SunPy — scientific visuals sit within a task-oriented tutorial/gallery.  
**V7 fix:** retain CAPIVARA's original logo; remove rainbow rail and temporary motif; one accent for navigation; reduce headline by roughly one-third. Make a real scientific figure visible in the first meaningful viewport while preserving readable text. On mobile the first action should be **Explore MaNGA example** or **Get started**, with documentation reachable without becoming the only prominent control.

### 2. Hero pair

**Observed:** the optical context image is a small or blown-up crop with source overlay, while the region image is a low-resolution raster derived from the manuscript. The pair lacks a scientific common coordinate/scale explanation. The small four-colour indicator below looks like a key but does not explain the 45 regions.  
**Relative to:** SunPy's scientific gallery and S2FFT's diagrams, where the image or computational property itself is the story.  
**V7 fix:** one physically registered, native-science plate on the same target; large real MaNGA field with calibrated intensity and **exact CAPIVARA region polygons** where the spatial contracts permit. Display optical context as separate labelled inset if image grids are not registered. Prominent panel title and metadata should state representation, region count, wavelength frame and source; avoid unlabelled colour dots.

### 3. VV 705 science map

**Observed:** the current source image contains patchy coloured raster cells with no individually selectable polygon layer, yet 45 regions are claimed. Controls switch between white light, spectral regions, H-alpha velocity and path regions, but there are **no per-region spectrum links** for this object. A gradient-like legend ornament is applied even to the categorical region view.  
**Relative to:** CARTA/Cubeviz already link regions and spectral profiles.  
**V7 fix:** rebuild polygons via existing R `regions_sf()` and exact label matrices, preserving holes and disconnected components, mapped to measured regional spectra from the **same galaxy**. Use `region_id` not colour positions as the selection key; selection should update spectra and region metadata. Different map modes get their own physically appropriate colourbar, units and captions; the region view never gets a continuous gradient.

### 4. Second example: MaNGA 10218-12703

**Observed:** measured spectra and aperture selection work, but they belong to a different MaNGA source and use aperture circles, not CAPIVARA regions. As the hero says “One galaxy. Many spectra,” this creates a conceptual discontinuity between a region figure for 7443 and linked spectra for 10218.  
**Relative to:** Photutils/Scikit-learn worked examples ordinarily keep the problem, input, method, result and verification together.  
**V7 fix:** establish a primary single-object viewer from 7443 that links its **own** segmentation to its **own** measured spectra. Retain the 10218 apertures in a separate gallery article, labelled as independent observational aperture extraction.

### 5. Mathematical section

**Observed:** the scalar positive-amplitude invariance formula is mathematically valid when the amplitude functional is positive homogeneous, but appears after the gallery with limited connection to a measured decision; no linked figure demonstrates what the choice changes.  
**Relative to:** S2FFT introduces an explicit mathematical problem, a diagram, an algorithmic consequence and a tradeoff.  
**V7 fix:** place `F/A(F)` beside a real before/after spectral-shape example, with retained amplitude and changed region membership. Separate a reference page on established equivariance from an experimental page on instrument incompatibility; include counterexamples, baseline citations, numerical tests and failure cases.

### 6. Mobile interaction, typography and information architecture

**Observed:** oversized display serif headings, repeated all-caps kickers and large blocks of explanatory text overtake plots. Four selector buttons wrap on mobile; important plot information sits far from the active controls. The live site exposes a lengthy synthetic setup prominently and does not yet provide a coherent R/Python entry route.  
**Relative to:** Astropy's stable user/reference separation, SunPy's example gallery, CARTA's linked views.  
**V7 fix:** a restrained serif headline only, Inter for all scientific labels and UI, 13–15 px readable mobile data text, grouped tabs with visible active selection, and a close-by linked plot. Give the first screen one real-data science promise; move longer methodology into science guide with source-level equations and tests. Navigation: Explore · Learn · Mathematics · API · Cite; R live, Python in development until validated.

## C. Final flagship layout decision

**Desktop:** white header, short natural title (max 2 lines), massive single MaNGA scientific plate taking a clear majority of the width, exact selectable region polygons in the map on the left, linked spectrum and compact measurement metadata on the right. The measured region ID and spectral axis are the shared anchor. Below it: separate modes for H-alpha velocity, kinematics, volume/slice and spectral-feature PCA, each labelled and constrained.

**Phone:** compact header, one-sentence scientific introduction, large legible map, region selector/tooltip and full-width spectrum immediately below. Do not force a dashboard into two tiny columns; interaction changes the single science plate. Keep a static export and an accessible alternative.

**Non-negotiable gate:** The first *live* demo needs actual region label matrices, spectral arrays and source data/provenance for the same target. This audit and colour specimen are not evidence that those products or new group-equivariant operators have been implemented or validated.

## D. Decision

Proceed with **white V7 editorial canvas + science-specific colour encoding + real polygon/spectrum linkage**. Do not publish recoloured manuscript pixels as new labels. Keep design and benchmarking on the draft RFC branch; no stable package code, live website or scientific release is changed by this document.
