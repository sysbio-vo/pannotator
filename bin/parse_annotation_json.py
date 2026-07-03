#!/usr/bin/env python3

import argparse
from pathlib import Path

import utils as ut


def extract_protein_annotations(annotations_json: dict) -> dict:
    assert "features" in annotations_json

    result = {protein["aa_hexdigest"]: protein for protein in annotations_json["features"]}

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein annotations from a Bakta JSON and write an index keyed by aa_hexdigest"
    )
    parser.add_argument(
        "-a",
        "--bakta_annotation_json",
        type=Path,
        required=True,
        help="Path to Bakta annotation JSON file",
    )
    parser.add_argument(
        "-o",
        "--proteins_json",
        type=Path,
        required=True,
        help="Output path for the proteins index JSON (keyed by aa_hexdigest)",
    )
    args = parser.parse_args()

    protein_annotations_dict = extract_protein_annotations(ut.load_json(args.bakta_annotation_json))
    ut.dump_json(protein_annotations_dict, args.proteins_json)
