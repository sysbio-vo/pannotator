#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from typing import Any, Dict

from bakta.constants import CDS_NOT_TRUNCATED
from bakta.utils import calc_aa_hash
from Bio import SeqIO


def load_annotations(json_path: Path) -> Dict[str, Dict[str, Any]]:
    with open(json_path, "r") as f:
        return json.load(f)


def dump_json(data: dict, outpath: Path) -> None:
    with open(outpath, "w") as f:
        json.dump(data, f, indent=4)


def annotate_bulk_proteins_with_auxDB(fasta_path: Path, auxiliary_db: dict) -> dict:

    looked_up_annotations = {}
    remaining_records = []

    # Open the file and parse it
    for record in SeqIO.parse(fasta_path, "fasta"):

        aa_hexdigest = calc_aa_hash(str(record.seq))[1]

        protein_annotation = auxiliary_db.get(aa_hexdigest, None)

        # print(f"{aa_hexdigest = } DEBUG {record.description = }")

        # TODO: take truncation information into account !!!
        # for now if a protein is not explicitly defined as not-truncated we don't fetch annotation for it
        if protein_annotation is None or CDS_NOT_TRUNCATED not in record.description:
            remaining_records.append(record)
            continue

        looked_up_annotations[aa_hexdigest] = protein_annotation

    return looked_up_annotations, remaining_records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Lookup annotations in an auxiliary JSON database and return a sliced JSON object"
    )
    parser.add_argument(
        "--auxiliary_db",
        type=Path,
        required=True,
        help="Path to auxiliary annotation JSON file. Please Refer to \
              Pannotator documentation for the auxDB format requirements",
    )

    parser.add_argument("--proteins_fa", type=Path, required=True, help="Input proteins to annotate in FASTA format")

    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output annotation JSON compatible with Pannotator bulk JSON annotation format",
    )

    parser.add_argument(
        "--remaining_proteins_filename",
        type=Path,
        required=True,
        help="Output FASTA for the proteins not found in the auxDB and to be passed to Bakta for annotation",
    )

    args = parser.parse_args()

    auxDB = load_annotations(args.auxiliary_db)
    looked_up_annotations, remaining_records = annotate_bulk_proteins_with_auxDB(args.proteins_fa, auxDB)

    # dump looked-up annotations in a Pannotator-compatible JSON
    dump_json(looked_up_annotations, args.out)

    # write remaining FASTA records to a new file
    with open(args.remaining_proteins_filename, "w") as output_handle:
        count = SeqIO.write(remaining_records, output_handle, "fasta")
    print(f"Saved {count} records to {args.remaining_proteins_filename} for Bakta annotation")
