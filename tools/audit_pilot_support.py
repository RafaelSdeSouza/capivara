"""Diagnose inherited spatial support without changing any V2 partition.

Continuum intensity, DRP eligibility and numerical fitting support are different
quantities. None supplies a certified galaxy-membership boundary here.
"""
from pathlib import Path
import csv
import hashlib
import json
import warnings

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap

P = Path(__file__).resolve().parents[3]
O = P / 'results/science_baseline_freeze_v3'
D = O / 'sandra_pilot/support_audit'
D.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 170})


def write_csv(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


region_rows, spaxel_rows, summaries, provenance = [], [], [], []
for path in sorted((O / 'physical').glob('*_provenance.json')):
    prov = json.loads(path.read_text())
    galaxy = prov['galaxy_id']
    with np.load(O / 'physical' / f'{galaxy}_spaxels.npz') as a:
        partition = a['cluster_map']
    records_path = O / 'sandra_pilot' / galaxy / 'regions.csv'
    records = {int(r['region_id']): r for r in csv.DictReader(records_path.open())}
    available = sorted(k for k, r in records.items() if r['fit_status'].startswith('SUCCESS'))
    unsupported = sorted(set(records) - set(available))
    inherited = np.isfinite(partition)
    fit_support = np.isin(partition, available)
    with fits.open(prov['source_path'], memmap=True, uint=False) as hdus:
        wave = np.asarray(hdus['WAVE'].data, dtype=float)
        rest = wave / (1 + prov['redshift'])
        indices = np.flatnonzero((rest >= 5050) & (rest <= 5500))
        selection = slice(indices[0], indices[-1] + 1)
        flux = np.asarray(hdus['FLUX'].data[selection], dtype=float)
        ivar = np.asarray(hdus['IVAR'].data[selection], dtype=float)
        mask = np.asarray(hdus['MASK'].data[selection], dtype=np.int64)
        assert flux.shape[1:] == partition.shape
        valid = np.isfinite(flux) & np.isfinite(ivar) & (ivar > 0) & (mask == 0)
        dontuse = (mask & 1024) != 0  # SDSS MANGA_DRP3PIXMASK bit 10
        lowcov = (mask & 2) != 0     # bit 1
        forestar = (mask & 8) != 0   # bit 3
        valid_fraction = valid.mean(axis=0)
        dontuse_fraction = dontuse.mean(axis=0)
        # All-invalid pixels remain missing; no zero-filled sky or interpolation.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='All-NaN slice encountered')
            continuum = np.nanmedian(np.where(valid, flux, np.nan), axis=0)
            formal_snr = np.nanmedian(np.where(valid, flux * np.sqrt(np.maximum(ivar, 0)), np.nan), axis=0)
        for region, row in records.items():
            spatial = partition == region
            n = int(spatial.sum())
            assert n == int(row['n_spaxels'])
            region_rows.append(dict(
                galaxy_id=galaxy, region_id=region, n_spaxels=n,
                numerical_fit_available=region in available,
                fit_status=row['fit_status'],
                display_class='FITTING_SUPPORT' if region in available else 'INSUFFICIENT_SUPPORT',
                galaxy_membership='NOT_VALIDATED',
                valid_continuum_sample_fraction=float(valid[:, spatial].mean()),
                donotuse_continuum_sample_fraction=float(dontuse[:, spatial].mean()),
                lowcov_continuum_sample_fraction=float(lowcov[:, spatial].mean()),
                forestar_continuum_sample_fraction=float(forestar[:, spatial].mean()),
                median_valid_continuum=float(np.nanmedian(continuum[spatial])) if np.isfinite(continuum[spatial]).any() else np.nan,
                median_formal_spaxel_snr=float(np.nanmedian(formal_snr[spatial])) if np.isfinite(formal_snr[spatial]).any() else np.nan,
                n_spaxels_with_any_donotuse=int(np.any(dontuse[:, spatial], axis=0).sum()),
                n_spaxels_with_complete_continuum_support=int(np.all(valid[:, spatial], axis=0).sum()),
                native_min=float(wave[indices[0]]), native_max=float(wave[indices[-1]]),
                rest_vacuum_min=float(rest[indices[0]]), rest_vacuum_max=float(rest[indices[-1]])))
        for y, x in zip(*np.where(inherited)):
            spaxel_rows.append(dict(galaxy_id=galaxy, x=int(x), y=int(y),
                region_id=int(partition[y, x]), numerical_fit_available=bool(fit_support[y, x]),
                valid_fraction=float(valid_fraction[y, x]), donotuse_fraction=float(dontuse_fraction[y, x]),
                median_valid_continuum=float(continuum[y, x]), median_formal_snr=float(formal_snr[y, x])))
        summaries.append(dict(galaxy_id=galaxy, inherited_regions=len(records),
            fitting_support_regions=len(available), insufficient_support_regions=len(unsupported),
            inherited_spaxels=int(inherited.sum()), fitting_support_spaxels=int(fit_support.sum()),
            insufficient_support_spaxels=int((inherited & ~fit_support).sum()),
            valid_sample_fraction_fitting_support=float(valid[:, fit_support].mean()),
            valid_sample_fraction_insufficient_support=float(valid[:, inherited & ~fit_support].mean()),
            fitting_support_ids=';'.join(map(str, available)), insufficient_support_ids=';'.join(map(str, unsupported))))
        provenance.append(dict(galaxy_id=galaxy, source_sha256=prov['source_sha256'],
            segmentation_sha256=prov['segmentation_sha256'],
            fit_records_sha256=hashlib.sha256(records_path.read_bytes()).hexdigest(),
            redshift=prov['redshift'], rest_vacuum_requested=[5050, 5500],
            native_selected=[float(wave[indices[0]]), float(wave[indices[-1]])],
            rest_vacuum_selected=[float(rest[indices[0]]), float(rest[indices[-1]])],
            n_channels=len(indices), fits_mask_ext=hdus['MASK'].name,
            drp_version=prov['drp_version']))

    np.savez_compressed(D / f'{galaxy}_support_maps.npz', inherited_partition=partition,
        fitting_support=fit_support, continuum_median_valid=continuum,
        continuum_valid_fraction=valid_fraction, continuum_donotuse_fraction=dontuse_fraction)
    fig, axs = plt.subplots(1, 4, figsize=(13, 3.9), layout='constrained', sharex=True, sharey=True)
    finite = continuum[np.isfinite(continuum) & (continuum > 0)]
    lo, hi = np.quantile(finite, [.05, .995])
    im = axs[0].imshow(np.ma.masked_less_equal(continuum, 0), origin='lower', cmap='magma',
                       norm=LogNorm(vmin=lo, vmax=hi), interpolation='nearest')
    fig.colorbar(im, ax=axs[0], shrink=.67, label='Median eligible flux density')
    axs[0].set_title('Eligible continuum')
    im = axs[1].imshow(np.where(np.isfinite(continuum) | inherited, dontuse_fraction, np.nan),
                       origin='lower', cmap='Greys', vmin=0, vmax=1, interpolation='nearest')
    fig.colorbar(im, ax=axs[1], shrink=.67, label='DONOTUSE sample fraction')
    axs[1].set_title('Native DRP quality')
    axs[2].imshow(partition, origin='lower', cmap='tab20', vmin=.5, vmax=20.5, interpolation='nearest')
    axs[2].set_title('Inherited V2 partition')
    axs[3].imshow(np.where(inherited & ~fit_support, 1., np.nan), origin='lower',
                  cmap=ListedColormap(['#e3e3e3']), vmin=0, vmax=1, interpolation='nearest')
    axs[3].imshow(np.where(fit_support, partition, np.nan), origin='lower',
                  cmap='tab20', vmin=.5, vmax=20.5, interpolation='nearest')
    axs[3].set_title('Numerical fitting support')
    for region in records:
        y, x = np.where(partition == region)
        centre = np.argmin((x - np.median(x))**2 + (y - np.median(y))**2)
        for ax in [axs[2]] + ([axs[3]] if region in available else []):
            ax.contour(partition == region, levels=[.5], colors='white', linewidths=.3)
            ax.text(x[centre], y[centre], str(region), ha='center', va='center', fontsize=6)
    for ax in axs:
        ax.set_xlabel('Native x pixel')
    axs[0].set_ylabel('Native y pixel')
    fig.savefig(O / 'figures' / f'{galaxy}_support_audit.png')
    plt.close(fig)
    print(galaxy, available, 'insufficient:', unsupported, flush=True)

assert len(region_rows) == 90 and len(summaries) == 5
assert sum(r['numerical_fit_available'] for r in region_rows) == 62
write_csv(D / 'region_support.csv', region_rows)
write_csv(D / 'spaxel_support.csv', spaxel_rows)
write_csv(D / 'summary.csv', summaries)
(D / 'provenance.json').write_text(json.dumps(dict(
    purpose='Diagnose inherited support; no segmentation or galaxy-membership inference.',
    quality_policy='Finite FLUX, finite positive IVAR, MASK==0; per-pixel continuum median over eligible samples.',
    fitting_display_policy='Colour only fixed regions with SUCCESS_INTERIOR or SUCCESS_BOUND_HIT; grey marks insufficient wavelength support. All population QC fails.',
    spatial_support_status='NOT_VALIDATED_AS_GALAXY_MEMBERSHIP',
    mask_reference='https://www.sdss4.org/dr14/algorithms/bitmasks/',
    continuum_units='1E-17 erg/s/cm^2/Angstrom/spaxel',
    formal_snr_limitation='Per-spaxel diagonal diagnostic; no regional or spectral covariance correction.',
    inputs=provenance), indent=2))
