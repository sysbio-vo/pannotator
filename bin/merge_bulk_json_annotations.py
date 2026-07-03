#!/usr/bin/env python3

import argparse
import operator
from functools import reduce
from pathlib import Path
from typing import Iterable

import utils as ut


def merge_jsons(json_paths: Iterable[Path]):
    jsons_list = [ut.load_json(json_path) for json_path in json_paths]
    return reduce(operator.ior, jsons_list, {})


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

    merged_annotation = merge_jsons(args.input_annotation_jsons)
    print(f"Saving merged annotation to {args.out}")
    ut.dump_json(merged_annotation, args.out)
