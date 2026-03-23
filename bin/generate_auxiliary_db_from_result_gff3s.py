#!/usr/bin/env python3

import argparse
import json

# import pickle
# import sqlite3
from pathlib import Path
from typing import Any, Dict


# TODO: move to a helper module and import from there
def load_annotations(json_path: Path) -> Dict[str, Dict[str, Any]]:
    with open(json_path, "r") as f:
        return json.load(f)


def split_annotations_into_lookup_and_alignment(annotations: Dict[str, Dict[str, Any]]) -> None:
    pass


def create_db() -> None:
    pass


def main():
    p = argparse.ArgumentParser(description="Annotate CDS pickles with information from JSON annotations")
    p.add_argument("--annotations", required=True, help="bulk_protein_annotations.json file")
    p.add_argument("--out", default="auxiliary_db.db", help="Output auxiliary DB path")
    args = p.parse_args()

    json_path = Path(args.annotations)
    # output_path = Path(args.out)

    if not json_path.exists():
        print(f"Annotations JSON file does not exist: {json_path}")
        return

    split_annotations_into_lookup_and_alignment()


if __name__ == "__main__":
    main()
