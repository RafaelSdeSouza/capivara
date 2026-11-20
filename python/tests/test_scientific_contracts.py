import numpy as np
import pytest

from capivara import IFUCube, SegmentationResult, region_spectra, segment_legacy


def cube(flux, variance=None, valid=None):
    a = np.asarray(flux, dtype=float)
    return IFUCube(a, np.arange(1, a.shape[-1] + 1), variance=variance,
                   valid=valid)


def test_wavelength_is_strictly_increasing():
    with pytest.raises(ValueError):
        IFUCube(np.ones((2, 2, 3)), [1, 3, 2])


def test_missing_flux_does_not_become_zero_measurement():
    c = cube(np.array([[[10., np.nan]], [[5., np.nan]]]))
    r = SegmentationResult(np.zeros((2, 1), dtype=int), c, "fixture", "")
    s = region_spectra(r)
    assert s["sum_flux"][0, 0] == 15
    assert np.isnan(s["sum_flux"][0, 1])
    assert s["finite_counts"][0, 1] == 0


def test_partial_variance_does_not_understate_error():
    c = cube(np.array([[[10.]], [[5.]]]),
             variance=np.array([[[4.]], [[np.nan]]]))
    r = SegmentationResult(np.zeros((2, 1), dtype=int), c, "fixture", "")
    s = region_spectra(r)
    assert s["sum_flux"][0, 0] == 15
    assert np.isnan(s["sum_variance_diagonal"][0, 0])


def test_full_variance_sum_is_correct_under_diagonal_assumption():
    c = cube(np.array([[[10.]], [[5.]]]),
             variance=np.array([[[4.]], [[9.]]]))
    r = SegmentationResult(np.zeros((2, 1), dtype=int), c, "fixture", "")
    assert region_spectra(r)["sum_variance_diagonal"][0, 0] == 13


def test_legacy_preserves_pixel_geometry_and_support():
    f = np.array([[[1., 2., 3.], [2., 1., 3.]],
                  [[2., 3., 1.], [-1., -2., -3.]]])
    s = segment_legacy(cube(f), n_regions=2)
    assert s.labels.shape == (2, 2)
    assert s.labels[1, 1] == -1
    assert len(np.unique(s.labels[s.labels >= 0])) == 2


def test_legacy_one_spaxel():
    f = np.array([[[1., 2.]], [[-1., -3.]]])
    assert segment_legacy(cube(f), n_regions=1).labels.tolist() == [[0], [-1]]
