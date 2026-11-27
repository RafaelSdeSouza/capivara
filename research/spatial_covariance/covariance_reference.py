"""Reference spatial covariance propagation for CAPIVARA V3 regional operators.

Experimental; not a production fitting backend. Native wavelength samples are
assumed independent, but spaxels can be correlated at each wavelength. The
caller supplies a covariance (not correlation) for each native wavelength,
in the exact spaxel ordering of the retained LSF operators. No interpolation
of MaNGA's four broadband correlation extensions is implied. A full input
spectral covariance would require a more general operator contract.
"""
import numpy as np
from scipy.sparse import csr_matrix, hstack, diags


def regional_covariance(operators, native_covariances, *, valid_native=None,
                        valid_output=None, check_psd=False):
    """Compute sum_l H_l C_l H_l.T with sparse operators.

    operators: one (n_output, n_native) matrix per spaxel; e.g. retained
      Gaussian LSF convolution matrices from prepare_region.
    native_covariances: iterable of n_native (n_spaxel, n_spaxel) matrices;
      covariance at each native wavelength, before validity masking.
    valid_native: optional boolean (n_spaxel, n_native) eligibility array.
    valid_output: optional boolean vector to select valid output channels;
      returns the corresponding principal submatrix, not zero-variance channels.
    check_psd: optionally check the input eigenvalues (dense, small cases only).
      Caller must establish positive semidefiniteness for production inputs.
    """
    ops = [csr_matrix(a, dtype=float) for a in operators]
    if not ops or len({a.shape for a in ops}) != 1:
        raise ValueError('Require nonempty same-shape spaxel operators')
    n_out, n_wave = ops[0].shape
    n_spax = len(ops)
    if not n_out or not n_wave or any(not np.isfinite(a.data).all() for a in ops):
        raise ValueError('Operators must be finite and nonempty')
    if valid_native is None:
        m = np.ones((n_spax, n_wave), dtype=bool)
    else:
        m = np.asarray(valid_native)
        if m.dtype != bool or m.shape != (n_spax, n_wave):
            raise ValueError('Invalid native eligibility')
    if valid_output is None:
        keep = np.arange(n_out)
    else:
        mask = np.asarray(valid_output)
        if mask.dtype != bool or mask.shape != (n_out,):
            raise ValueError('Invalid output eligibility')
        keep = np.flatnonzero(mask)
    ops = [a[keep, :].tocsc() for a in ops]
    out = csr_matrix((keep.size, keep.size), dtype=float)
    covs = iter(native_covariances)
    for l in range(n_wave):
        try:
            c = csr_matrix(next(covs), dtype=float)
        except StopIteration as exc:
            raise ValueError('Too few native covariance matrices') from exc
        if c.shape != (n_spax, n_spax) or not np.isfinite(c.data).all():
            raise ValueError(f'Invalid covariance at native channel {l}')
        scale = max(float(np.max(np.abs(c.data), initial=0.0)), np.finfo(float).tiny)
        delta = c - c.T
        if np.max(np.abs(delta.data), initial=0.0) > 1e-10 * scale:
            raise ValueError(f'Asymmetric covariance at native channel {l}')
        c = (c + c.T) * 0.5
        if np.any(c.diagonal() < 0):
            raise ValueError(f'Negative variance at native channel {l}')
        if check_psd and np.min(np.linalg.eigvalsh(c.toarray())) < -1e-10 * scale:
            raise ValueError(f'Non-PSD covariance at native channel {l}')
        # H_l contains the l-th column of each spaxel response. Zero the
        # column for an invalid native measurement before propagation.
        h = hstack([a[:, l] for a in ops], format='csr') @ diags(m[:, l].astype(float))
        out = out + h @ c @ h.T
    sentinel = object()
    if next(covs, sentinel) is not sentinel:
        raise ValueError('Too many native covariance matrices')
    return ((out + out.T) * 0.5).tocsr()


def self_test():
    # Different per-spaxel LSF operators with correlated input spaxels.
    a = np.array([[1., .2, 0.], [0., .8, .3], [0., 0., .7]])
    b = np.array([[.8, 0., 0.], [.2, .7, .1], [0., .3, .9]])
    cs = [np.array([[1+l, .2], [.2, 2+l]]) for l in range(3)]
    ref = sum(np.column_stack([a[:, l], b[:, l]]) @ cs[l] @
              np.column_stack([a[:, l], b[:, l]]).T for l in range(3))
    actual = regional_covariance([a, b], cs, check_psd=True).toarray()
    np.testing.assert_allclose(actual, ref, atol=1e-12)
    # Independent-spaxel limit must reproduce current V3 diagonal-native formula.
    diag = [np.diag(np.diag(c)) for c in cs]
    expected = a @ np.diag([c[0, 0] for c in cs]) @ a.T + b @ np.diag([c[1, 1] for c in cs]) @ b.T
    np.testing.assert_allclose(regional_covariance([a, b], diag).toarray(), expected, atol=1e-12)
    # Spatially correlated identity response, variance of two-spaxel sum.
    ident = np.eye(3)
    np.testing.assert_allclose(regional_covariance([ident, ident], cs).toarray(),
                               np.diag([c.sum() for c in cs]), atol=1e-12)
    # Eligibility and principal-submatrix output, not zero-filled outputs.
    m = np.ones((2, 3), bool); m[1, 1] = False
    bmasked = b.copy(); bmasked[:, 1] = 0
    expected = sum(np.column_stack([a[:, l], bmasked[:, l]]) @ cs[l] @
                   np.column_stack([a[:, l], bmasked[:, l]]).T for l in range(3))
    np.testing.assert_allclose(regional_covariance([a, b], cs, valid_native=m,
        valid_output=np.array([True, False, True])).toarray(), expected[np.ix_([0, 2], [0, 2])], atol=1e-12)
    # Reject an invalid covariance, rather than silently regularizing it.
    try:
        regional_covariance([a,b], [np.array([[1.,2.],[2.,1.]])]*3, check_psd=True)
    except ValueError:
        pass
    else:
        raise AssertionError('Indefinite input was accepted')
    print('PASS: 5 reference checks (dense equality, diagonal limit, sum, masks, PSD rejection)')


if __name__ == '__main__':
    self_test()
