"""Rebuild the compact 15-galaxy MaNGA collection from frozen iFUN products.

Usage: python build_ifun_collection.py IFUN_CAPIVARA_EAT_MANGA NEW_OUTPUT_DIRECTORY
The source cubes and all temporary products remain outside Git.
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ifun_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    root = args.ifun_root.resolve()
    out = args.output_dir
    if out.exists():
        raise FileExistsError(f"Output exists: {out}")
    manifest = (root / "results/sandra46_regional_master_v1/provenance/"
                "cube_sources.csv")
    with manifest.open(newline="") as stream:
        targets = sorted(r["plateifu"] for r in csv.DictReader(stream))
    if len(targets) != 15 or len(set(targets)) != 15:
        raise ValueError("Expected the 15 frozen iFUN MaNGA sources")
    out.mkdir(parents=True)
    stage = out / ".staging"
    stage.mkdir()
    collection = []
    for target in targets:
        target_out = stage / target
        subprocess.run([sys.executable, str(Path(__file__).with_name("reuse_ifun_manga.py")),
                        target, str(root), str(target_out)], check=True)
        payload = json.loads((target_out / f"{target}.json").read_text())
        collection.append(dict(object=target, n_regions=payload["actual_k"],
                               n_channels=len(payload["wavelength"]),
                               n_spaxels=sum(r["n_spaxels"] for r in payload["regions"]),
                               spectrum_status="measured native fixed-aperture sums",
                               physical_fit_status="not run in this Atlas"))
        for path in target_out.iterdir():
            path.rename(out / path.name)
        target_out.rmdir()
    stage.rmdir()
    (out / "collection.json").write_text(
        json.dumps(dict(schema="capivara-observed-collection-v1",
                        source="Frozen iFUN MaNGA DR17 regional products",
                        entries=collection), ensure_ascii=False,
                   allow_nan=False, separators=(",", ":")) + "\n", encoding="utf-8")
    print(f"Published {len(collection)} checked observed MaNGA entries in {out}")


if __name__ == "__main__":
    main()
