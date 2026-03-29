#!/usr/bin/env python3

import argparse
import json
import operator
from functools import reduce
from pathlib import Path
from typing import Iterable


def load_json(infile: str) -> dict:
    with open(infile, "r") as f:
        res_dict = json.load(f)
    return res_dict


def dump_json(data: dict, outpath: Path) -> None:
    with open(outpath, "w") as f:
        json.dump(data, f, indent=4)


def merge_bulk_annotation_jsons(json_paths: Iterable[Path]):
    annotation_dicts = [load_json(json_path) for json_path in json_paths]
    return reduce(operator.ior, annotation_dicts, {})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein annotations from a Bakta JSON and write an index keyed by aa_hexdigest"
    )
    parser.add_argument(
        "-i",
        "--input_annotation_jsons",
        type=Path,
        nargs="+",
        required=True,
        help="Path to Bakta annotation JSON file",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        required=True,
        help="Output path for the merged protein annotation JSON (keyed by aa_hexdigest)",
    )

    args = parser.parse_args()

    merged_annotation = merge_bulk_annotation_jsons(args.input_annotation_jsons)
    print(f"Saving merged annotation to {args.out}")
    dump_json(merged_annotation, args.out)
