# Capivara website content map

This map records the Phase B migration from the former Quarto site and drifting
README sources to one pkgdown site. No scientific claim is declared obsolete.

| Source content | Destination | Treatment | Parity notes |
|---|---|---|---|
| `index.qmd` hero, purpose, primary MaNGA mosaic | `README.qmd` / homepage | Retained and shortened | Semantic H1, compact logo, corrected arXiv identifier, no broken badge wall; mosaic is fully visible and has a provenance limitation. |
| `index.qmd` four feature cards | Homepage open workflow plus backend, support, and products guides | Redistributed | Cube, support, regions, and regional spectra remain distinct; nested cards removed. |
| Installation and optional Torch note | Homepage and Getting started | Retained | Homepage keeps the standard install command; first run does not require Torch. |
| Local-file-only minimal workflow | Getting started; observed-data template in Examples | Replaced and retained by role | A fixed-seed synthetic cube is the executable first run. The FITS-path workflow remains an explicitly non-executable user-data template. |
| Full-spectrum and feature-window explanation | Choosing a segmentation backend | Moved | Preserves the difference between `feature_wavelength_range` and S/N `wavelength_range`. |
| Segmentation plotting, NA/zero labels, FITS writing, spatial WCS | FITS and DS9 export | Moved | Code, label semantics, WCS caveat, and map-to-table correspondence retained. |
| Exact-backend memory growth and `segment_large()` | Choosing a segmentation backend | Moved | Memory estimate, sparse graph parameters, backend metadata, and limitations retained. |
| Starlet support controls | Support masks | Moved | All homepage controls retained; adaptive support added from the documented API. |
| Core API list | Generated Reference index | Replaced by generated navigation | Every established documented export is grouped by scientific role. Pending emission-line exports are deliberately not advertised. |
| Companion packages and scope | Homepage continuation links and relevant guides | Condensed | Capivara remains the segmentation layer; specialised fitting remains separate. |
| `manga_8443_6102_compare_current.png` | Examples | Moved | Caption and purpose retained; exact provenance remains unresolved. |
| Homepage BibTeX | `inst/CITATION`, pkgdown citation metadata, Paper nav | Canonicalised | Complete ten-author journal record; DOI and arXiv identifier verified. |
| Newer kinematic material present only in `README.md` | Existing Kinematic analysis and Bisymmetric bar models articles | Preserved | Native kinematic, panel, and explicit bar-hypothesis concepts remain available without dominating first use. |
| `about.qmd` project description, repository, paper | Homepage `#about`, sidebar, Paper and GitHub nav | Moved and redirected | `/about.html` targets `index.html#about`. |
| `news.qmd` include of `NEWS.md` | Generated pkgdown news | Replaced and redirected | `/news.html` targets `/news/index.html`. |
| `styles.css` site styling | `pkgdown/extra.css` | Reimplemented | Retains package blue/orange identity with narrower Bootstrap corrections and native theme support. |
| Existing `kinematic-analysis.Rmd` | Same article route | Retained | Only navigation and interpretive cross-links added. |
| Existing `bisymmetric-bar-model.Rmd` | Same article route | Retained | Only navigation and the no-detection caveat added. |

## Redirect and anchor parity

- `/about.html` → `/index.html#about`, preserving incoming query strings and
  fragments where supplied.
- `/news.html` → `/news/index.html`, preserving incoming query strings and
  fragments.
- Homepage anchors retained: `#installation`, `#quick-start`, and `#core-api`.

The former Quarto source files may be retired only after the pkgdown build,
redirect generation, link validation, and rendered content review pass.
