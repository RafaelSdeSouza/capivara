#!/usr/bin/env python3
"""DR17 pilot preflight. No fits, segmentation, or inferred galaxy membership."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ID = "8602-12705"
DAPTYPE = "HYB10-MILESHC-MASTARSSP"
C = 299792.458
DRP = "https://data.sdss.org/sas/dr17/manga/spectro/redux/v3_1_1/8602/stack/"
DAP = ("https://data.sdss.org/sas/dr17/manga/spectro/analysis/v3_1_1/3.1.0/"
       + DAPTYPE + "/8602/12705/")
CUBE_NAME = f"manga-{ID}-LOGCUBE.fits"
MAPS_NAME = f"manga-{ID}-MAPS-{DAPTYPE}.fits.gz"


def digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def clean(a):
    a = np.asarray(a, dtype=float)
    return np.where(np.isfinite(a), a, None).tolist()


def masked(value, mask, ivar):
    value, mask, ivar = map(np.asarray, (value, mask, ivar))
    if not value.shape == mask.shape == ivar.shape:
        raise ValueError("MAPS value/mask/IVAR shapes differ.")
    good = np.isfinite(value) & (mask == 0) & np.isfinite(ivar) & (ivar > 0)
    return np.where(good, value, np.nan)


def intrinsic(raw, correction):
    raw, correction = map(np.asarray, (raw, correction))
    q = raw**2 - correction**2
    ok = (raw >= 0) & (correction >= 0) & (q > 0)
    return np.where(ok, np.sqrt(np.maximum(q, 0)), np.nan)


def main():
    from astropy import __version__ as astropy_version
    from astropy.io import fits
    from astropy.wcs import WCS

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cube", type=Path, help="Uncompressed DR17 DRP LOGCUBE")
    p.add_argument("maps", type=Path, help="Matching DR17 DAP MAPS")
    p.add_argument("output", type=Path)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    with fits.open(args.cube, memmap=True) as c, fits.open(args.maps, memmap=False) as m:
        for h in (c[0].header, m[0].header):
            if str(h.get("PLATEIFU", "")).strip() != ID:
                raise ValueError("Unexpected target in FITS header")
        if str(m[0].header.get("DAPTYPE", "")).strip() != DAPTYPE:
            raise ValueError("Unexpected DAPTYPE")
        for key, h in (("DRP3QUAL", c[0].header), ("DAPQUAL", m[0].header)):
            if key not in h or int(h[key]) & (1 << 30):
                raise ValueError(f"Missing or critical {key}")
        z = float(m[0].header["SCINPVEL"]) / C
        if not np.isfinite(z) or not 0 < z < 1:
            raise ValueError("Invalid DAP input redshift")
        flux, ivar, mask = (c[k].data for k in ("FLUX", "IVAR", "MASK"))
        wave = np.asarray(c["WAVE"].data, dtype=float).ravel()
        if not (flux.shape == ivar.shape == mask.shape) or flux.ndim != 3:
            raise ValueError("Incompatible native cubes")
        if flux.shape[0] != len(wave) or not np.all(np.isfinite(wave)) or not np.all(np.diff(wave) > 0):
            raise ValueError("Invalid WAVE extension")
        ny, nx = flux.shape[1:]
        if m["BINID"].data.shape[1:] != (ny, nx):
            raise ValueError("DRP/DAP shapes differ")
        w1, w2 = WCS(c["FLUX"].header).celestial, WCS(m["BINID"].header).celestial
        if not (w1.has_celestial and w2.has_celestial):
            raise ValueError("Celestial WCS missing")
        yy, xx = np.indices((ny, nx), dtype=float)
        sep = w1.pixel_to_world(xx, yy).separation(w2.pixel_to_world(xx, yy)).arcsec
        if not np.isfinite(sep).all() or float(sep.max()) > 1e-3:
            raise ValueError("DRP/DAP WCS mismatch; no reprojection")
        # BINID channel 1: analysed binning footprint, not proof of galaxy membership.
        aperture = m["BINID"].data[0] >= 0
        if aperture.sum() < 16:
            raise ValueError("Insufficient DAP binning support")
        with (args.output / "aperture.csv").open("w", newline="") as stream:
            w = csv.writer(stream)
            w.writerow(["axis1", "axis2", "included"])
            w.writerows((x+1, y+1, int(aperture[y, x])) for y in range(ny) for x in range(nx))
        # Asymmetric flux/IVAR witnesses detect transposition in the R FITS reader.
        # Use three well separated channels and a deterministic set of sampled pixels.
        witness = []
        pixels = np.argwhere(aperture)
        for k in (len(wave)//3, len(wave)//2, 2*len(wave)//3):
            for y, x in pixels[np.linspace(0, len(pixels)-1, 12, dtype=int)]:
                f, v, q = float(flux[k, y, x]), float(ivar[k, y, x]), int(mask[k, y, x])
                if np.isfinite(f) and np.isfinite(v) and v > 0 and q == 0:
                    witness.append((int(x)+1, int(y)+1, k+1, float(wave[k]), f, v))
        if len(witness) < 6 or not any(x != y for x, y, *_ in witness):
            raise ValueError("Insufficient orientation witnesses")
        with (args.output / "witness.csv").open("w", newline="") as stream:
            w = csv.writer(stream)
            w.writerow(["axis1", "axis2", "channel", "wavelength", "flux", "ivar"])
            w.writerows(witness)
        meta = dict(object=ID, release="SDSS DR17", daptype=DAPTYPE,
                    redshift=z, redshift_source="DAP PRIMARY SCINPVEL / 299792.458 (input reference)",
                    drp_url=DRP+CUBE_NAME+".gz", dap_url=DAP+MAPS_NAME,
                    cube_sha256=digest(args.cube), maps_sha256=digest(args.maps),
                    cube_hash_representation="uncompressed FITS",
                    nx=nx, ny=ny, n_wave=len(wave),
                    wcs_max_separation_arcsec=float(sep.max()),
                    aperture_pixels=int(aperture.sum()), aperture="DAP BINID channel 1 >= 0; analysis footprint only",
                    drp3qual=int(c[0].header["DRP3QUAL"]), dapqual=int(m[0].header["DAPQUAL"]),
                    numpy=np.__version__, astropy=astropy_version)
        (args.output/"preflight.json").write_text(json.dumps(meta, indent=2, allow_nan=False)+"\n")
        # Independent survey reference maps. Never label these as CAPIVARA fits.
        result = []
        def add(name,title,group,kind,unit,values,note):
            if values.shape != (ny,nx):
                raise ValueError("Reference map dimensions differ")
            result.append(dict(id=name,title=title,group=group,kind=kind,unit=unit,
                               note=note,values=clean(values)))
        def q(name,channel=None):
            arrays=[np.asarray(m[name+s].data) for s in ("","_MASK","_IVAR")]
            return masked(*(a if channel is None else a[channel] for a in arrays))
        add("dap_stellar_velocity","DAP stellar velocity","kinematics","diverging","km/s",
            q("STELLAR_VEL"),"Survey reference; stellar measurements repeated across DAP bins, relative to DAP input redshift.")
        add("dap_stellar_sigma","DAP stellar dispersion","kinematics","continuous","km/s",
            intrinsic(q("STELLAR_SIGMA"),m["STELLAR_SIGMACORR"].data[0]),
            "Survey reference; template-resolution corrected in quadrature; unresolved values withheld.")
        h=m["EMLINE_GFLUX"].header
        names=[str(h.get(f"C{i+1:02d}","")).strip() for i in range(m["EMLINE_GFLUX"].data.shape[0])]
        # Include all named channels within the viewer's cap; DAP gas velocities are tied.
        for i,name in enumerate(names):
            if not name:
                raise ValueError("Missing DAP line channel label")
            add(f"dap_flux_{i:02d}",f"DAP {name} flux","emission","continuous","10^-17 erg/s/cm²/spaxel",
                q("EMLINE_GFLUX",i),"Survey Gaussian line flux; Galactic extinction corrected; strict MASK==0 and positive IVAR.")
        oiii=[i for i,name in enumerate(names) if "OIII" in name.upper() and ("5007" in name or "5008" in name)]
        if len(oiii)==1:
            i=oiii[0]
            add("dap_gas_velocity","DAP gas velocity","kinematics","diverging","km/s",
                q("EMLINE_GVEL",i),"Survey reference, [O III] channel; velocities tied across lines, not independent per-line fits.")
            add("dap_oiii_sigma","DAP [O III] dispersion","kinematics","continuous","km/s",
                intrinsic(q("EMLINE_GSIGMA",i),m["EMLINE_INSTSIGMA"].data[i]),
                "Survey reference; instrumental correction in quadrature; unresolved widths withheld.")
        (args.output/"dap_reference.json").write_text(json.dumps(result,allow_nan=False,separators=(",",":")))
        print(json.dumps(meta,indent=2))


if __name__ == "__main__":
    main()
