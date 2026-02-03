#!/usr/bin/env python3

import argparse
import json
import pickle
from pathlib import Path
from typing import Dict


def find_pkl_files(input_dir: Path) -> Dict[str, Path]:
    pkls = {}
    for pkl in input_dir.glob("*.fa.cds-only.pkl"):
        sample = pkl.name.split(".fa.cds-only.pkl")[0]
        pkls[sample] = pkl
    return pkls


def load_pkl(path: Path) -> dict:
    with path.open("rb") as fh:
        return pickle.load(fh)


def iter_cdss(data: dict):
    for feat in data.get("features", []):
        if feat.get("type") == "cds":
            yield feat


def build_index(input_dir: Path) -> dict:
    index = {}
    pkls = find_pkl_files(input_dir)

    for sample, pkl_path in pkls.items():
        data = load_pkl(pkl_path)

        for cds in iter_cdss(data):
            aa_seq = cds.get("aa")
            aa_hash = cds.get("aa_hexdigest")

            if not aa_seq or not aa_hash:
                continue

            entry = index.setdefault(
                aa_hash,
                {
                    "aa_seq": aa_seq,
                    "members": [],
                },
            )

            entry["members"].append(
                {
                    "sample": sample,
                    "contig": cds["sequence"],
                    "locus": cds.get("locus"),
                    "start": cds.get("start"),
                    "end": cds.get("stop") or cds.get("end"),
                    "strand": cds.get("strand"),
                }
            )

    return index


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build CDS-coords index from cds-only bakta results")
    parser.add_argument("input_dir", type=Path, help="Folder with *.cds-only.pkl files")
    parser.add_argument("output_json", type=Path, help="Path to write resulting JSON file")
    args = parser.parse_args()

    print(f"Building CDS-coords index from files in {args.input_dir}")
    index = build_index(args.input_dir)

    print(f"Writing results to {args.output_json}")
    args.output_json.write_text(json.dumps(index, indent=2))
