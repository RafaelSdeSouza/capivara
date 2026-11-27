# Atlas example development notes

Internal planning and validation checklist. Candidate archive paths, links and targets must be verified before publication.

# CAPIVARA example gallery

The installed [Get Started](https://rafaelsdesouza.com.br/capivara/articles/getting-started.html) guide is the first stop. Existing observed examples are in the [MaNGA examples vignette](../../vignettes/examples.Rmd), including the published [segmentation mosaic](../../man/figures/mosaic_segmented_sagui.png) and [8443-6102 backend comparison](../../man/figures/manga_8443_6102_compare_current.png). The [V4 validation report](../../SPATIAL_FOUNDATION_V4_FREEZE.md) describes the five real MaNGA replay targets. These are existing segmentation results, not new multi-survey physical fits.

## Reproduce the examples

See [Analyse your own data](../../vignettes/your-own-data.Rmd) and the installed scripts in [inst/tutorials](../../inst/tutorials). The synthetic demo does not need an archive download. The observational starter accepts a local MaNGA DRP cube or a prepared RDS from any instrument. No survey FITS cubes are committed.

## Planned flagship observational atlas

| Target | Programme / instrument | Science demonstration | Data source | Status |
| --- | --- | --- | --- | --- |
| 8602-12705 | MaNGA / SDSS | Disc kinematics and spatial segmentation | [SDSS DR17](https://www.sdss4.org/dr17/manga/manga-data/data-access/) | Eight frozen iFUN regions with measured native spectra in the Atlas; no new physical fit |
| NGC 628 | PHANGS / MUSE | H II regions and spiral disc | [PHANGS data releases](https://www.phangs.org/data) | Candidate; exact files, crop and validation pending |
| NGC 5485 | CALIFA / PMAS-PPAK | Early-type stellar kinematics | [CALIFA DR3](https://califa.caha.es/DR3/) | Candidate; exact files and validation pending |
| JO206 | GASP / MUSE | Ram-pressure stripping, extraplanar emission | [ESO GASP release](https://www.eso.org/cms/eso-archive-news/second-data-release-of-the-large-programme--dissecting-gas-strip.html) | Candidate; exact files and validation pending |
| NGC 5643 | MAGNUM / MUSE | Nuclear ionized gas and outflows | [ESO Science Archive](https://archive.eso.org/cms.html), programme 60.A-9339; [science paper](https://doi.org/10.1051/0004-6361/201526581) | Candidate; verify final reduction and access |
| CATID 272822 | SAMI / AAT | Stellar–gas kinematic misalignment | [Data Central](https://datacentral.org.au/) | Candidate; catalogue cross-match, exact files and validation pending |

## Numerical observed MaNGA collection

The Atlas contains a [15-galaxy numerical MaNGA collection](../../inst/atlas/objects/collection.json), including [8602-12705](../../inst/atlas/objects/8602-12705.json). Each entry reuses eight frozen iFUN full-flux regions and their regional spectra from verified DR17 linear DRP HDUs. The compact entries retain 2,201 native observed vacuum channels from 4800 to 7000 Å, native pixel labels, per-channel valid counts, diagonal variance and the saved spatial covariance-adjusted variance. A regional flux is displayed as a fixed-aperture sum only when every region spaxel passes finite flux, positive inverse variance and the MaNGA DONOTUSE mask bit in that channel. Other DRP mask bits follow the frozen iFUN selection rule. The shaded band in the viewer uses the saved spatial covariance-adjusted variance; spectral covariance was not reconstructed. These entries contain no new velocity, dispersion, line fit or physical classification.

Run `python research/atlas_examples/build_ifun_collection.py IFUN_CAPIVARA_EAT_MANGA NEW_OUTPUT_DIRECTORY` to rebuild all 15 entries, or `python research/atlas_examples/reuse_ifun_manga.py TARGET IFUN_CAPIVARA_EAT_MANGA NEW_OUTPUT_DIRECTORY` for one target. The exporter checks archived file hashes, target identity, spatial and wavelength grids, and flux units. It independently recomputes every selected channel's regional valid count, flux sum and diagonal variance from the source cube. Each output includes Atlas JSON, a native-pixel SVG partition and two descriptive catalogue CSV files. Source cubes stay outside this repository. The separate `prepare_observed_pilot.py` path starts from a MaNGA LOGCUBE or CALIFA cube and has different spectral sampling and region definitions; its results have not been substituted for the frozen iFUN collection.

For the MaNGA target, the DR17 archive naming convention gives the following paths (availability of the exact files has not been checked in this environment):

- [DRP LOGCUBE](https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/8602/stack/manga-8602-12705-LOGCUBE.fits.gz)
- [DAP MAPS HYB10-MILESHC-MASTARSSP](https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/HYB10-MILESHC-MASTARSSP/8602/12705/manga-8602-12705-MAPS-HYB10-MILESHC-MASTARSSP.fits.gz)

Do not interpret DAP measurements as CAPIVARA fits. Do not transfer existing MEGACUBE labels to LOGCUBE/MAPS unless pixel alignment and WCS have been verified.

## Publishing policy

Each finished example should include an archive identifier, release, citation, exact file hashes, configuration, software versions, wavelength/flux conventions, masks and WCS validation, and a description of which physical products passed quality gates. Commit small example scripts, provenance and validated web-ready products only; large native cubes stay in the survey archives. MUSE mosaics should initially use a documented science-driven crop, retaining the WCS and separate coverage mask. Large-cube processing and every non-MaNGA instrument require validation before claims of support. Unavailable products must be visibly unavailable, not simulated or zero-filled.

The multi-survey interactive atlas, unified joint fitting, calibrated anomaly detection and semantic catalogue are development goals, not features of this example index.
