# Capivara website visual audit

Status: **local implementation complete**. Publication is ready after the
updated workflow dependency is committed and pushed. The previous deployment
failed because `rsconnect` was absent from the CI website library; it is now
declared in `.github/workflows/publish.yml`.

Local preview: `python3 -m http.server 8765 --bind 127.0.0.1 --directory docs`

## Rendered QA

| Page | Viewport | Theme | Screenshot | H1 | Overflow | Image alt | Result |
|---|---:|---|---|---:|---|---|---|
| Home | 1440 × 900 | light | `docs/website-audit/screenshots/capivara-home-desktop-light.png` | 1 | none | complete | pass |
| Home | 1440 × 900 | dark | `docs/website-audit/screenshots/capivara-home-desktop-dark.png` | 1 | none | complete | pass |
| Home | 390 × 844 | light | `docs/website-audit/screenshots/capivara-home-mobile-light.png` | 1 | none | complete | pass |
| Getting Started | 1440 × 900 | light | `docs/website-audit/screenshots/capivara-getting-started-desktop-light.png` | 1 | none | complete | pass |
| Getting Started | 390 × 844 | light | `docs/website-audit/screenshots/capivara-getting-started-mobile-light.png` | 1 | none | complete | pass |
| Examples | 1440 × 900 | light | `docs/website-audit/screenshots/capivara-examples-desktop-light.png` | 1 | none | complete | pass |
| Reference | 1440 × 900 | light | `docs/website-audit/screenshots/capivara-reference-desktop-light.png` | 1 | none | complete | pass |
| Support Masks | 1440 × 900 | light | `docs/website-audit/screenshots/capivara-support-masks-desktop-light.png` | 1 | none | complete | pass |
| Support Masks | 768 × 1024 | light | browser inspection | 1 | none | complete | pass |

## Functional checks

- Exact primary navigation order: Get started, User guide, Examples, Reference,
  Paper, GitHub.
- Mobile navigation opens and exposes all top-level links; User guide remains a
  labelled submenu.
- Native light and dark theme controls both render without overflow.
- Search for `segment_large` returns the function reference and related guides.
- Keyboard focus on search has a visible 3 px orange outline.
- Browser console recorded no errors or warnings.
- Internal-link checker passed 37 HTML files, including fragments and redirect
  targets.
- `/about.html` and `/news.html` preserve query strings and incoming fragments.

## Scientific visuals

- The real MaNGA mosaic renders at its natural aspect ratio without cropping.
- The fixed-seed first run uses a 24 × 24 × 96 line-rich IFS simulation and six
  categorical regions. Its map and regional spectra use a restrained
  publication palette and no ordered-label interpretation.
- File-level provenance is recorded in `docs/website_image_provenance.csv`.
  Two retained real panels remain explicitly marked `incomplete-holding`; their
  producer/configuration mapping has not been invented.

## Known repository condition

A non-mutating documentation audit found pre-existing drift from pending
emission-line roxygen blocks: running `devtools::document()` in a temporary copy
would reorder one existing NAMESPACE export and create `emission_lines.Rd` and
`segment_emission_lines.Rd`. Those pending exports were present before Phase B,
are omitted from the website reference, and were not documented or advertised
merely to make the gate pass.
