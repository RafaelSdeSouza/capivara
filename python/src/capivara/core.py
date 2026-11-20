"""Explicit observational contracts for a small IFU data cube.

This is a Python-port prototype, not a claim of R/Python numerical parity.
Axis convention is [row, column, spectral_channel]. Wavelengths are native;
no redshift correction, air/vacuum transformation or resampling is implicit.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class IFUCube:
    """An IFU observation with explicitly declared coordinates and validity.

    flux: measured flux [row, col, wavelength]; not necessarily calibrated to
        a physical flux unit unless flux_unit says so.
    wavelength: strictly increasing native-axis values in wavelength_unit.
    variance: optional diagonal elementwise variance; spatial correlations
        are NOT represented and cannot be silently presumed absent.
    valid: optional Boolean voxel mask; True means eligible measured sample.
        Invalid voxels are retained in flux but not used for summaries.
    frame: either 'observed' or 'rest'. Conversion is never automatic.
    """

    flux: np.ndarray
    wavelength: np.ndarray
    variance: np.ndarray | None = None
    valid: np.ndarray | None = None
    frame: str = "observed"
    wavelength_unit: str = "Angstrom"
    wavelength_medium: str = "vacuum"
    flux_unit: str = "unspecified"

    def __post_init__(self) -> None:
        flux = np.asarray(self.flux, dtype=float)
        wave = np.asarray(self.wavelength, dtype=float)
        if flux.ndim != 3 or min(flux.shape) == 0:
            raise ValueError("flux must be nonempty [row, col, wavelength] array")
        if wave.ndim != 1 or len(wave) != flux.shape[-1]:
            raise ValueError("wavelength must have one value per spectral channel")
        if not np.all(np.isfinite(wave)) or not np.all(np.diff(wave) > 0):
            raise ValueError("wavelength must be finite and strictly increasing")
        if self.frame not in {"observed", "rest"}:
            raise ValueError("frame must be 'observed' or 'rest'")
        if self.wavelength_medium not in {"vacuum", "air"}:
            raise ValueError("wavelength_medium must be 'vacuum' or 'air'")
        if self.valid is None:
            valid = np.isfinite(flux)
        else:
            valid = np.asarray(self.valid, dtype=bool)
            if valid.shape != flux.shape:
                raise ValueError("valid mask must match flux shape")
            valid = valid & np.isfinite(flux)
        variance = None
        if self.variance is not None:
            variance = np.asarray(self.variance, dtype=float)
            if variance.shape != flux.shape:
                raise ValueError("variance must match flux shape")
        object.__setattr__(self, "flux", flux)
        object.__setattr__(self, "wavelength", wave)
        object.__setattr__(self, "valid", valid)
        object.__setattr__(self, "variance", variance)


@dataclass(frozen=True)
class SegmentationResult:
    """Region labels and provenance. A label of -1 means unassigned."""

    labels: np.ndarray
    cube: IFUCube
    method: str
    notes: str

    def __post_init__(self) -> None:
        labels = np.asarray(self.labels, dtype=np.int64)
        if labels.shape != self.cube.flux.shape[:2]:
            raise ValueError("labels must match the spatial dimensions")
        object.__setattr__(self, "labels", labels)


def region_spectra(result: SegmentationResult) -> dict[str, np.ndarray]:
    """Flux summaries with per-channel measurement counts.

    Given regional membership C and valid measured samples V, compute
    S_C(lambda) = sum_{p in C and V} flux[p, lambda].
    Channels with no contribution remain NaN, not zero.
    If diagonal variances exist, sum variances only if *every* measured flux
    included in that sum has a valid positive variance. Otherwise NaN.
    This is not a full spatial-covariance propagation.
    """
    cube = result.cube
    ids = np.unique(result.labels[result.labels >= 0])
    n = len(cube.wavelength)
    sums = np.full((len(ids), n), np.nan)
    means = np.full_like(sums, np.nan)
    counts = np.zeros((len(ids), n), dtype=np.int64)
    var_sum = np.full_like(sums, np.nan) if cube.variance is not None else None

    for i, label in enumerate(ids):
        pixels = result.labels == label
        f = cube.flux[pixels]
        v = cube.valid[pixels]
        counts[i] = v.sum(axis=0)
        s = np.where(v, f, 0.0).sum(axis=0)
        ok = counts[i] > 0
        sums[i, ok] = s[ok]
        means[i, ok] = s[ok] / counts[i, ok]
        if var_sum is not None:
            vv = cube.variance[pixels]
            has_var = v & np.isfinite(vv) & (vv > 0)
            complete = ok & (has_var.sum(axis=0) == counts[i])
            candidate = np.where(has_var, vv, 0.0).sum(axis=0)
            var_sum[i, complete] = candidate[complete]

    return {
        "region_ids": ids,
        "wavelength": cube.wavelength.copy(),
        "sum_flux": sums,
        "mean_flux": means,
        "finite_counts": counts,
        "sum_variance_diagonal": var_sum,
    }
