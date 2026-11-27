"""Download and prepare two public observational pilot cubes for CAPIVARA.
Only native observed channels 4800--7000 Angstrom are retained.
No interpolation, physical fitting or automatic host-membership claims.
Requires numpy, astropy, requests. Large source files remain outside Git.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

SDSS = "https://data.sdss.org/sas/dr17/manga/spectro"
CALIFA = "https://califa.caha.es/FTP-PUB/reduced/V500/reduced_v2.2"
MANGA_TARGETS = {"8602-12705", "8443-6102"}


def download(url, path):
    if path.exists():
        return
    partial = path.with_suffix(path.suffix + ".partial")
    with requests.get(url, stream=True, timeout=(30, 180)) as response:
        response.raise_for_status()
        n = 0
        with partial.open("wb") as f:
            for block in response.iter_content(chunk_size=1024 * 1024):
                n += len(block)
                if n > 2_000_000_000:
                    raise ValueError("Unexpected download size; stop.")
                f.write(block)
    partial.replace(path)


def sha256(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def verify_identity(h, identity):
    if str(h.get("PLATEIFU", "")).strip() != identity:
        raise ValueError("PLATEIFU missing or differs from requested target")


def manga(target, root):
    plate, ifu = target.split("-")
    drp = f"manga-{target}-LOGCUBE.fits.gz"
    daptype = "HYB10-MILESHC-MASTARSSP"
    dap = f"manga-{target}-MAPS-{daptype}.fits.gz"
    url = f"{SDSS}/redux/v3_1_1/{plate}/stack/{drp}"
    ref = f"{SDSS}/analysis/v3_1_1/3.1.0/{daptype}/{plate}/{ifu}/{dap}"
    p, q = root/drp, root/dap
    download(url, p)
    download(ref, q)
    with fits.open(q) as h:
        verify_identity(h[0].header, target)
        z = float(h[0].header["SCINPVEL"]) / 299792.458
    with fits.open(p, memmap=False) as h:
        verify_identity(h[0].header, target)
        if int(h[0].header.get("DRP3QUAL", 0)) & (1 << 30):
            raise ValueError("MaNGA DRP CRITICAL quality flag")
        wave = np.asarray(h["WAVE"].data, dtype=float)
        use = (wave >= 4800) & (wave <= 7000)
        f = np.asarray(h["FLUX"].data[use], dtype=np.float32)
        iv = np.asarray(h["IVAR"].data[use], dtype=np.float32)
        dq = h["MASK"].data[use]
        valid = np.isfinite(f) & np.isfinite(iv) & (iv > 0) & (dq == 0)
        var = np.full_like(f, np.nan)
        np.divide(1., iv, out=var, where=valid)
        unit = h["FLUX"].header["BUNIT"]
    meta = dict(object_id=target, survey="SDSS-IV MaNGA DR17",
                redshift=z, redshift_source="DAP MAPS PRIMARY SCINPVEL / c",
                wavelength_medium="vacuum", flux_unit=unit,
                cube_url=url, cube_sha256=sha256(p), metadata_url=ref,
                metadata_sha256=sha256(q),
                mask_policy="MASK==0 and finite positive IVAR and finite FLUX",
                archive_release="DR17 v3_1_1 / DAP 3.1.0")
    return wave[use], f, var, valid, meta


def califa(root):
    p = root/"NGC5485.V500.rscube.fits.gz"
    url = f"{CALIFA}/{p.name}"
    download(url, p)
    with fits.open(p, memmap=False) as h:
        header = h[0].header
        obj = str(header.get("OBJECT", "")).upper().replace(" ", "")
        if obj != "NGC5485":
            raise ValueError(f"CALIFA OBJECT not confirmed: {obj!r}")
        ctype = str(header.get("CTYPE3", "")).upper().strip()
        if ctype not in ("WAVE", "AWAV"):
            raise ValueError(f"CALIFA wavelength medium unverified ({ctype!r}); inspect source metadata before continuing.")
        medium = "air" if ctype == "AWAV" else "vacuum"
        # Require explicitly declared spectral units and separable spectral WCS.
        if not header.get("CUNIT3"):
            raise ValueError("CALIFA wavelength unit absent; inspect header.")
        sw = WCS(header).spectral
        raw = np.asarray(h[0].data)
        wave = (sw.pixel_to_world_values(np.arange(raw.shape[0])) *
                u.Unit(sw.world_axis_units[0])).to_value(u.AA)
        use = (wave >= 4800) & (wave <= 7000)
        f = np.asarray(raw[use], dtype=np.float32)
        err = np.asarray(h["ERROR"].data[use], dtype=np.float32)
        bad = h["BADPIX"].data[use]
        valid = np.isfinite(f) & np.isfinite(err) & (err > 0) & (bad == 0)
        var = np.where(valid, err**2, np.nan)
        unit = header["BUNIT"]
    meta = dict(object_id="NGC5485", survey="CALIFA V500 DR3 v2.2",
                redshift=0.00637901,
                redshift_source="https://califa.caha.es/Observed_Objects.html",
                wavelength_medium=medium, flux_unit=unit,
                cube_url=url, cube_sha256=sha256(p),
                mask_policy="BADPIX==0 and finite positive ERROR and finite flux",
                archive_release="CALIFA V500 reduced_v2.2")
    return wave[use], f, var, valid, meta


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("target", choices=sorted(MANGA_TARGETS | {"NGC5485"}))
    ap.add_argument("directory", type=Path)
    a = ap.parse_args()
    a.directory.mkdir(parents=True, exist_ok=True)
    wave, f, var, valid, meta = (manga(a.target, a.directory)
        if a.target in MANGA_TARGETS else califa(a.directory))
    assert f.ndim == 3 and f.shape == var.shape == valid.shape
    assert len(wave) == f.shape[0] and len(wave) >= 100
    assert np.all(np.isfinite(wave)) and np.all(np.diff(wave)>0)
    assert np.isfinite(meta["redshift"]) and meta["redshift"] > -1
    valid &= np.isfinite(var) & (var > 0)
    var[~valid] = np.nan
    # A continuum-selected analysis footprint, NOT a galaxy membership mask.
    anchor = (wave >= 5400) & (wave <= 5500)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = np.nanmedian(np.where(valid[anchor], f[anchor] / np.sqrt(var[anchor]), np.nan), axis=0)
    support = (valid.mean(axis=0) >= .7) & np.isfinite(snr) & (snr >= 3)
    if support.sum() < 16:
        raise ValueError("Insufficient continuum-selected spatial support")
    meta.update(wavelength_frame="observed", wavelength_range=[float(wave[0]),float(wave[-1])],
                wavelength_selection="Native observed channels 4800--7000 Angstrom; no resampling",
                spatial_selection="At least 70% eligible samples; median per-channel S/N >=3 in observed 5400--5500 A. Not a validated host boundary.",
                native_grid_shape=list(f.shape), support_count=int(support.sum()),
                covariance="Not propagated; physical fits withheld")
    fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(f,name="FLUX"),
        fits.ImageHDU(var,name="VARIANCE"),fits.ImageHDU(valid.astype(np.uint8),name="VALID"),
        fits.ImageHDU(wave,name="WAVE"),fits.ImageHDU(support.astype(np.uint8),name="SUPPORT")
        ]).writeto(a.directory/"prepared.fits",overwrite=False,checksum=True)
    (a.directory/"source.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    print(json.dumps(meta,indent=2))

if __name__ == "__main__":
    main()
