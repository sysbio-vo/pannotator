#!/usr/bin/env python3

import argparse
import json
import os
import pickle
from pathlib import Path
from typing import Iterable

import bakta.constants as bc

"""
All bakta feature fields for CDS, pseudogene and sORF features
derived from a test run output
{'aa',
 'aa_digest',
 'aa_hexdigest',
 'db_xrefs',
 'end',
 'frame',
 'gene',
 'genes',
 'hypothetical',
 'id',
 'ips',
 'locus',
 'nt',
 'pfams',
 'product',
 'psc',
 'pscc',
 'pseudo-inference',
 'pseudogene',
 'rbs_motif',
 'seq_stats',
 'sequence',
 'start',
 'start_type',
 'stop',
 'strand',
 'truncated',
 'type',
 'ups'}


 All fields from Bakta bulk protein anno (CDS only, no sORF, no pseudogenes)
 {'aa',
 'aa_hexdigest',
 'db_xrefs',
 'description', # contains information abount truncation, not utilised currently for auxDB
 'expert',
 'gene',
 'genes',
 'hypothetical',
 'id',
 'ips',
 'length',
 'locus',
 'pfams',
 'product',
 'psc',
 'pscc',
 'seq_stats',
 'type',
 'ups'}
"""

CACHED_FEATURE_TYPES = (bc.FEATURE_CDS, bc.FEATURE_SORF, bc.PSEUDOGENE)

MANDATORY_BAKTA_FEATURE_FIELDS = ("aa", "type", "db_xrefs")
BAKTA_FEATURE_FIELDS = (
    "expert",
    "gene",
    "genes",
    "hypothetical",
    "ups",
    "ips",
    "psc",
    "pscc",
    "length",
    "pfams",
    "product",
    "seq_stats",
    "pseudogene",
)  # CDS + postprocessing (pseudogene, pfams, expert) features only for now


def read_pickle(pickle_path: str) -> dict:
    with open(pickle_path, "rb") as fh:
        pickled_obj = pickle.load(fh)
    return pickled_obj


def load_json(json_path: Path) -> dict:
    with open(json_path, "r") as f:
        return json.load(f)


def dump_json(data: dict, outpath: Path) -> None:
    with open(outpath, "w") as f:
        json.dump(data, f, indent=4)


def sample_feature_to_annotation_entry(bakta_feature: dict) -> dict:

    protein_annotation = {}

    for mandatory_field in MANDATORY_BAKTA_FEATURE_FIELDS:
        protein_annotation[mandatory_field] = bakta_feature[mandatory_field]

    for feature_field in BAKTA_FEATURE_FIELDS:
        if feature_field in bakta_feature:
            protein_annotation[feature_field] = bakta_feature[feature_field]

    return protein_annotation


def collect_annotations(pickle_paths: Iterable[Path]) -> dict:
    unique_features = {}
    for pickle_path in pickle_paths:
        sample_data = read_pickle(pickle_path)
        for feature in sample_data["features"]:
            if feature["type"] not in CACHED_FEATURE_TYPES:
                continue
            if feature["aa_hexdigest"] not in unique_features:
                protein_anno = sample_feature_to_annotation_entry(feature)
                unique_features[feature["aa_hexdigest"]] = protein_anno
    return unique_features


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein annotations from a Bakta JSON and write an index keyed by aa_hexdigest"
    )
    parser.add_argument(
        "-i",
        "--annotation_pickles",
        type=Path,
        nargs="+",
        required=True,
        help="List of annotation pickle objects to extract CDS features from and generate auxiliary DB",
    )

    parser.add_argument(
        "-a",
        "--auxiliary_db",
        type=Path,
        required=True,
        help="Path to existing JSON DB to update with new proteins",
    )

    parser.add_argument(
        "-o",
        "--updated_db_out",
        type=Path,
        required=True,
        help="Output path for updated auxiliary JSON DB",
    )

    args = parser.parse_args()

    bulk_annotations = collect_annotations(args.annotation_pickles)
    if os.path.exists(args.auxiliary_db):
        print(f"Updating existing auxiliary DB {args.auxiliary_db}")
        existing_aux_db = load_json(args.auxiliary_db)
        bulk_annotations |= existing_aux_db
    else:
        print(f"Saving auxiliary DB annotation to {args.updated_db_out}")
    dump_json(bulk_annotations, args.updated_db_out)
