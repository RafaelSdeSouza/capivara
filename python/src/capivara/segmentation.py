"""Historical CAPIVARA segmentation baseline for parity investigations.

Important: R's legacy path applies Ward.D2 updates to Manhattan distances.
This is NOT minimization of ordinary Euclidean within-cluster SSE.
Do not silently replace the distance metric or claim the resulting hierarchy
is numerically identical across SciPy, fastcluster and different tie orders.
"""
from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import cut_tree, linkage
from scipy.spatial.distance import pdist

from .core import IFUCube, SegmentationResult


def segment_legacy(cube: IFUCube, n_regions: int = 15) -> SegmentationResult:
    """Reproduce the *algorithmic choices* of R's historical segment() baseline.

    Scientific question: which measured spaxel spectra are similar after their
    additive continuum level is removed?
    Features: subtract per-spaxel median, set unmeasured channels to zero.
    Distance: L1 / Manhattan, followed by Ward-style Lance-Williams updates.
    Eligibility: positive sum over finite measured channels.
    Limitations: zero fill encodes a missing-data assumption; no spatial
    connectivity guarantee; no physical stellar-population interpretation.
    R/Python numerical parity has not yet been certified with gold fixtures.
    """
    f = cube.flux
    v = cube.valid
    n_regions = int(n_regions)
    signal = np.where(v, f, 0.0).sum(axis=-1)
    eligible = v.any(axis=-1) & (signal > 0)
    rows = np.flatnonzero(eligible.reshape(-1))
    if n_regions < 1 or n_regions > len(rows):
        raise ValueError("n_regions must lie between 1 and eligible spaxels")
    labels = np.full(f.shape[:2], -1, dtype=np.int64)
    if len(rows) == 1:
        labels.ravel()[rows] = 0
    else:
        X = f.reshape(-1, f.shape[-1])[rows]
        V = v.reshape(-1, f.shape[-1])[rows]
        med = np.array([np.median(x[mask]) for x, mask in zip(X, V)])
        X = np.where(V, X - med[:, None], 0.0)
        if n_regions == len(rows):
            flat_labels = np.arange(len(rows))
        else:
            d = pdist(X, metric="cityblock")
            tree = linkage(d, method="ward")
            flat_labels = cut_tree(tree, n_clusters=n_regions).ravel()
        labels.ravel()[rows] = flat_labels

    return SegmentationResult(
        labels=labels, cube=cube, method="historical_L1_WardD2",
        notes=("Algorithmic baseline only; R/Python parity, ties, missingness "
               "and MaNGA science products require cross-language fixtures."),
    )
