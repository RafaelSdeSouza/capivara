# CAPIVARA for Python — experimental foundation

A genuine Python implementation is being developed alongside the mature R
package, not as a translation that silently changes the scientific objective.

**Status:** exploratory 0.0.1 development package. Only the algorithmic
historical segmentation baseline and elementary region-summary contracts are
present. **R/Python numerical parity has NOT been certified.** Do not use this
package for a published science result without validating the required path.

## Install for development

    python -m pip install -e "./python[test,visual]"
    python -m pytest python/tests

## Minimal use

    import numpy as np
    from capivara import IFUCube, segment_legacy, region_spectra
    from capivara.plotting import plot_diagnostic

    rng = np.random.default_rng(12)
    flux = rng.uniform(0.1, 1.0, size=(8, 8, 48))
    cube = IFUCube(flux=flux, wavelength=np.linspace(4800, 7400, 48))
    result = segment_legacy(cube, n_regions=4)
    summary = region_spectra(result)
    fig, axes = plot_diagnostic(result)
    fig.savefig("capivara_diagnostic.svg")

## Mathematical honesty

The historical R algorithm centres each spaxel by its median, computes
pairwise Manhattan distances and updates those dissimilarities using Ward.D2.
This does not generally minimize Euclidean within-cluster sums of squares.
SciPy and R can differ at ties or due to numerical conventions. A separate
Euclidean-variance objective must have its own explicit API and tests.

The mask is voxel-wise. No-contributor channels return missing (NaN), not
invented zero flux. Diagonal variance summation is provided only when every
contributing flux has a valid variance. Spatial covariance is not propagated.

## Before releasing

1. Golden R/Python partitions from the frozen R implementation, compared
   modulo label permutation, with strict wavelength and coordinate conventions.
2. Stable sample/voxel masks and missingness; exact spectral and spatial axes.
3. Sparse Ward port of the Rcpp algorithm with graph and tie diagnostics.
4. Real FITS/WCS/LSF adapters, reference MaNGA/MUSE observations, benchmarks.
5. Publication figures, validation documentation and a versioned reference site.

The R package remains the source of truth until the respective parity gate
passes. This directory does not replace the R implementation.
