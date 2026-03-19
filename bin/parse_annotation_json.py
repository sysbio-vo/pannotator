#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def extract_protein_annotations(annotations_json: dict) -> dict:
    assert "features" in annotations_json.keys()

    result = {protein["aa_hexdigest"]: protein for protein in annotations_json["features"]}

    return result


def load_json(infile: str) -> dict:
    with open(infile, "r") as f:
        res_dict = json.load(f)
    return res_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build CDS-coords index from cds-only bakta results")
    parser.add_argument(
        "-a", "--bakta_annotation_json", type=Path, help="Folder with *.cds-only.gff3 and *.cds-only.faa files"
    )
    parser.add_argument("-o", "--proteins_json", type=Path, help="Path to write resulting JSON file")
    args = parser.parse_args()

    protein_annotations_dict = extract_protein_annotations(load_json(args.bakta_annotation_json))
    args.proteins_json.write_text(json.dumps(protein_annotations_dict, indent=2))
