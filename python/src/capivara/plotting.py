"""Publication-oriented visual diagnostics without hidden scientific inference."""
from __future__ import annotations

import numpy as np

from .core import SegmentationResult, region_spectra


def plot_diagnostic(result: SegmentationResult, *, max_spectra: int = 7):
    """Four-panel audit of support, partitions, summed spectra and reconstruction.

    This reports computed data only; region colours are categorical, never
    physical ordinals. The residual compares observed flux with a within-region
    *mean* spectral model; it is not a new measurement or independent validation.
    Returns (figure, axes) for full caller control and vector export.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    cube = result.cube
    f = cube.flux
    valid = cube.valid
    white = np.where(valid, f, 0.0).sum(axis=-1)
    white[~valid.any(axis=-1)] = np.nan
    regions = result.labels
    summary = region_spectra(result)
    ids = summary["region_ids"]
    recon = np.full(f.shape, np.nan)
    for i, k in enumerate(ids):
        recon[regions == k] = summary["mean_flux"][i]

    residual = np.where(valid, f - recon, np.nan)
    residual_rms = np.sqrt(np.nanmean(residual**2, axis=-1))
    residual_rms[~valid.any(axis=-1)] = np.nan

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    ax = axes.ravel()
    ax[0].imshow(white, origin="lower", cmap="magma", interpolation="nearest")
    ax[0].set_title("Measured flux integrated over channels")
    n = max(len(ids), 1)
    cmap = ListedColormap(plt.get_cmap("tab20")(np.linspace(0, 1, n)))
    idx = np.full(regions.shape, np.nan)
    for i, k in enumerate(ids):
        idx[regions == k] = i
    ax[1].imshow(idx, origin="lower", cmap=cmap, vmin=-0.5, vmax=n-0.5,
                 interpolation="nearest")
    ax[1].set_title("Spectral region labels (categorical)")
    for i in range(min(len(ids), max_spectra)):
        ax[2].plot(cube.wavelength, summary["sum_flux"][i], lw=1.6,
                   label=f"Region {ids[i]}")
    if len(ids):
        ax[2].legend(frameon=False, fontsize=8, ncol=2)
    ax[2].set(xlabel=f"Wavelength ({cube.wavelength_unit}; {cube.frame})",
              ylabel=f"Summed flux ({cube.flux_unit})")
    ax[2].set_title("Flux-preserving regional spectra")
    ax[3].imshow(residual_rms, origin="lower", cmap="cividis",
                 interpolation="nearest")
    ax[3].set_title("RMS residual to regional mean spectrum")
    for k in (0, 1, 3):
        ax[k].set(xlabel="Column (pixel)", ylabel="Row (pixel)")
    fig.suptitle("CAPIVARA | Scientific diagnostic", fontsize=16)
    return fig, axes
