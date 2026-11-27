"""Publish a compact Atlas entry from frozen iFUN MaNGA regional products.

This does not segment or fit a cube. It checks the archived source cube,
region map and regional spectra before writing native-channel Atlas data.
Usage: python reuse_ifun_manga.py TARGET IFUN_CAPIVARA_EAT_MANGA OUTPUT_DIRECTORY
"""
import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from astropy.io import fits


PRODUCT = Path("results/sandra46_regional_master_v1")


def digest(path):
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def manifest_hash(manifest, relative):
    with manifest.open(newline="") as stream:
        matches = [r["sha256"] for r in csv.DictReader(stream)
                   if r["path"] == str(relative)]
    if len(matches) != 1:
        raise ValueError(f"Missing or ambiguous manifest entry: {relative}")
    return matches[0]


def finite(x):
    return float(x) if np.isfinite(x) else None


def record_metadata(path):
    code = ('x<-readRDS(commandArgs(TRUE)[1]);cat(x$plateifu,x$input_path,'
            'x$input_sha256,x$spatial_contract$redshift,'
            'x$spatial_contract$redshift_source,x$capivara_version,'
            'x$representation$type,sep="\\t")')
    fields = subprocess.run(["Rscript", "-e", code, str(path)],
                            check=True, capture_output=True, text=True).stdout.split("\t")
    if len(fields) != 7:
        raise ValueError("Unexpected frozen segmentation record")
    return fields


def adjacency(labels, spectra, ids):
    contacts = {}
    for left, right in ((labels[:, :-1], labels[:, 1:]),
                        (labels[:-1, :], labels[1:, :])):
        for a, b in zip(left.flat, right.flat):
            if a <= 0 or b <= 0 or a == b:
                continue
            pair = tuple(sorted((int(a), int(b))))
            contacts[pair] = contacts.get(pair, 0) + 1
    edges = []
    for (a, b), shared in sorted(contacts.items()):
        x, y = spectra[ids.index(a)], spectra[ids.index(b)]
        common = np.isfinite(x) & np.isfinite(y)
        fraction = float(common.mean())
        distance = None
        if common.sum() >= 3 and fraction >= .7:
            u, v = x[common], y[common]
            if np.max(np.abs(u)) > 0 and np.max(np.abs(v)) > 0:
                u = u / np.max(np.abs(u))
                v = v / np.max(np.abs(v))
                distance = float(np.clip(
                    1 - np.dot(u, v) / np.sqrt(np.dot(u, u) * np.dot(v, v)),
                    0, 2))
        edges.append(dict(region_a=a, region_b=b, shared_edges=shared,
                          common_fraction=fraction, cosine_distance=distance))
    return edges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target")
    parser.add_argument("ifun_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    root = args.ifun_root.resolve()
    out = args.output_dir
    target = args.target
    if not re.fullmatch(r"[0-9]{4,5}-[0-9]{4,5}", target):
        raise ValueError("Require a MaNGA plate-IFU identifier")
    labels_path = PRODUCT / f"segmentation/{target}/full_flux_k08_labels.fits"
    spectra_path = PRODUCT / f"spectra/{target}_regional_spectra.npz"
    record_path = PRODUCT / f"segmentation/{target}/full_flux_k08_segmentation_record.rds"
    source_rows = [r for r in csv.DictReader(
        (root / PRODUCT / "provenance/cube_sources.csv").open(newline=""))
        if r["plateifu"] == target]
    if len(source_rows) != 1:
        raise ValueError("Missing or ambiguous frozen source cube record")
    source = source_rows[0]
    cube_path = Path(source["path"])
    if out.exists():
        raise FileExistsError(f"Output exists: {out}")
    paths = [root / p for p in (cube_path, labels_path, spectra_path, record_path)]
    if not all(p.is_file() for p in paths):
        raise FileNotFoundError("The frozen iFUN cube, labels, spectra and record are required")
    source_manifest = root / PRODUCT / "provenance/source_manifest.csv"
    artifact_manifest = root / PRODUCT / "provenance/artifact_manifest.csv"
    plateifu, record_input, record_hash, redshift, redshift_source, version, representation = record_metadata(root / record_path)
    if (plateifu != target or record_input != str(cube_path) or
            record_hash != source["sha256"] or representation != "full_flux"):
        raise ValueError("Frozen segmentation record differs from source manifest")
    actual_hashes = {}
    for rel, path in zip((cube_path, labels_path, spectra_path, record_path), paths):
        expected = manifest_hash(
            source_manifest if rel == cube_path else artifact_manifest, rel)
        actual = digest(path)
        if actual != expected:
            raise ValueError(f"Frozen source hash changed: {rel}")
        actual_hashes[rel.name] = actual

    with fits.open(root / labels_path, memmap=True) as handle:
        labels = np.asarray(handle[0].data, dtype=np.int64)
    data = np.load(root / spectra_path)
    ids = [int(i) for i in data["region_id"]]
    sizes = data["n_spaxels"].astype(int)
    if ids != list(range(1, 9)) or len(labels.shape) != 2:
        raise ValueError("Unexpected region IDs or spatial grid")
    if any(np.count_nonzero(labels == i) != n for i, n in zip(ids, sizes)):
        raise ValueError("Region map and regional spectra disagree on membership")
    if not np.all((labels < 0) | np.isin(labels, ids)):
        raise ValueError("Unexpected spatial label")
    wave = data["wave_vacuum"]
    select = np.flatnonzero((wave >= 4800) & (wave <= 7000))
    if select.size < 100 or not np.all(np.diff(wave) > 0):
        raise ValueError("Unexpected native wavelength grid")

    with fits.open(root / cube_path, memmap=True) as handle:
        if handle[0].header.get("PLATEIFU") != target:
            raise ValueError("Source cube target differs")
        if int(handle[0].header.get("DRP3QUAL", 0)) & (1 << 30):
            raise ValueError("Critical DRP quality flag")
        if handle["FLUX"].data.shape != (len(wave), *labels.shape):
            raise ValueError("Region map does not match the native cube grid")
        if not np.array_equal(handle["WAVE"].data, wave):
            raise ValueError("Saved spectra and source wavelength grids differ")
        unit = handle["FLUX"].header.get("BUNIT")
        if unit != "1E-17 erg/s/cm^2/Angstrom/spaxel":
            raise ValueError(f"Unverified flux unit: {unit}")
        # Recompute only native per-channel sums and diagonal variance.
        # This is an independent numerical check, not a new segmentation.
        for channel in select:
            flux = handle["FLUX"].data[channel]
            ivar = handle["IVAR"].data[channel]
            mask = handle["MASK"].data[channel]
            # Match the frozen iFUN DRP DONOTUSE-bit rule. Other mask bits
            # remain recorded in the source cube and are not silently excluded.
            good = (np.isfinite(flux) & np.isfinite(ivar) &
                    (ivar > 0) & ((mask & 1024) == 0))
            for row, region in enumerate(ids):
                eligible = good & (labels == region)
                count = int(eligible.sum())
                if count != data["n_valid"][row, channel]:
                    raise ValueError(f"Valid-count mismatch at {region}, {channel}")
                if count:
                    observed = np.sum(flux[eligible], dtype=float)
                    diagonal = np.sum(1 / ivar[eligible], dtype=float)
                    if not np.isclose(observed, data["flux"][row, channel],
                                      rtol=1e-10, atol=1e-10):
                        raise ValueError(f"Flux-sum mismatch at {region}, {channel}")
                    if not np.isclose(diagonal,
                                      data["variance_diagonal"][row, channel],
                                      rtol=2e-7, atol=1e-9):
                        raise ValueError(f"Variance mismatch at {region}, {channel}")

    counts = data["n_valid"][:, select].astype(int)
    complete = counts == sizes[:, None]
    fixed = np.where(complete, data["flux"][:, select], np.nan)
    diagonal = np.where(complete, data["variance_diagonal"][:, select], np.nan)
    spatial = np.where(complete, data["variance_spatial"][:, select], np.nan)
    if np.any(np.isfinite(spatial) & (spatial < 0)):
        raise ValueError("Negative spatial covariance-adjusted variance")
    edges = adjacency(labels, fixed, ids)
    regions = []
    features = []
    for row, region in enumerate(ids):
        neighbours = [e for e in edges if region in (e["region_a"], e["region_b"])]
        distances = [e["cosine_distance"] for e in neighbours
                     if e["cosine_distance"] is not None]
        vals = fixed[row] / sizes[row]
        feature = dict(region=region, n_spaxels=int(sizes[row]),
                       complete_fraction=float(np.isfinite(fixed[row]).mean()),
                       median_flux_density_per_spaxel=finite(np.nanmedian(vals)),
                       n_neighbours=len(neighbours),
                       max_neighbour_shape_distance=max(distances) if distances else None)
        features.append(feature)
        regions.append(dict(**feature,
                            spectrum=[finite(v) for v in fixed[row]],
                            valid_counts=counts[row].tolist(),
                            variance_diagonal=[finite(v) for v in diagonal[row]],
                            variance_spatial=[finite(v) for v in spatial[row]],
                            median_lsf_sigma_angstrom=finite(
                                np.nanmedian(data["lsf_sigma"][row, select]))))
    payload = dict(
        schema="capivara-semantic-atlas-v1", object=target,
        source="Observed · frozen iFUN MaNGA regional product",
        nx=labels.shape[1], ny=labels.shape[0], representation="full_flux",
        wavelength_frame="observed", wavelength_medium="vacuum",
        flux_unit=unit, redshift=float(redshift), requested_k=8, actual_k=8,
        wavelength=wave[select].tolist(),
        labels=[[int(v) if v > 0 else None for v in row] for row in labels],
        regions=regions, edges=edges,
        provenance=dict(
            survey="SDSS-IV MaNGA DR17", archive_url=source["url"],
            cube_kind=source["acquisition"],
            input_sha256=actual_hashes,
            original_record=str(record_path),
            redshift_source=redshift_source,
            spatial_coordinates=f"Native {labels.shape[1]} x {labels.shape[0]} pixel grid; no celestial WCS exported",
            wavelength_selection="Native observed vacuum channels 4800--7000 Angstrom; no interpolation",
            flux_policy="Fixed-aperture sum shown only when all region spaxels pass finite FLUX, positive IVAR and MaNGA MASK DONOTUSE bit 1024 clear",
            uncertainty="Per-channel diagonal IVAR sum and frozen iFUN spatial covariance-adjusted variance; spectral covariance not reconstructed",
            lsf="Median regional post-pixel DRP sigma from the frozen iFUN product; no deconvolution",
            interpretation="Descriptive spectral regions; no new physical fit or calibrated anomaly probability",
            independent_check="Every selected native-channel flux sum, valid count and diagonal variance recomputed from frozen source HDUs",
            software=f"Frozen iFUN CAPIVARA {version} record; this exporter does not rerun segmentation"))
    out.mkdir(parents=True)
    (out / f"{target}.json").write_text(
        json.dumps(payload, ensure_ascii=False, allow_nan=False,
                   separators=(",", ":")) + "\n", encoding="utf-8")
    for name, rows in ((f"{target}-region_features.csv", features),
                       (f"{target}-region_adjacency.csv", edges)):
        with (out / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else
                                    ["region_a", "region_b", "shared_edges",
                                     "common_fraction", "cosine_distance"],
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    colours = ("#1E295B", "#E3D75A", "#4E74A6", "#B58B2B",
               "#68A69C", "#C96442", "#80538A", "#44312A")
    tiles = [f'<rect x="{x}" y="{y}" width="1" height="1" fill="{colours[v-1]}"/>'
             for y, row in enumerate(labels[::-1]) for x, v in enumerate(row) if v > 0]
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {labels.shape[1]} {labels.shape[0]}" '
           'role="img" aria-label="Native-pixel eight-region MaNGA partition">'
           f'<rect width="{labels.shape[1]}" height="{labels.shape[0]}" fill="#FCFCFA"/>' +
           ''.join(tiles) + '</svg>\n')
    (out / f"{target}-regions.svg").write_text(svg, encoding="utf-8")
    print(f"{target}: {len(ids)} regions, {len(select)} verified native channels, "
          f"{sum(sizes)} labelled spaxels -> {out}")


if __name__ == "__main__":
    main()
