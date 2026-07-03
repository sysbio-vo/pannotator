#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from typing import Iterable

import bakta.constants as bc
import utils as ut

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
    bc.HYPOTHETICAL_PROTEIN_NOT_PSEUDOGENE,
)  # CDS + postprocessing (pseudogene, pfams, expert) features only for now


def sample_feature_to_annotation_entry(bakta_feature: dict) -> dict:

    protein_annotation = {}

    for mandatory_field in MANDATORY_BAKTA_FEATURE_FIELDS:
        protein_annotation[mandatory_field] = bakta_feature[mandatory_field]

    for feature_field in BAKTA_FEATURE_FIELDS:
        if feature_field in bakta_feature:
            protein_annotation[feature_field] = bakta_feature[feature_field]

    return protein_annotation


def update_cds_annotation(pickle_paths: Iterable[Path], cds_annotation_before_filtering: dict) -> dict:
    for pickle_path in pickle_paths:
        sample_data = ut.load_pickle(pickle_path)  # TODO: change if other serialization formats will be added
        for feature in sample_data["features"]:
            if feature["type"] not in CACHED_FEATURE_TYPES:
                continue

            protein_anno = sample_feature_to_annotation_entry(feature)

            if feature["aa_hexdigest"] not in cds_annotation_before_filtering:
                cds_annotation_before_filtering[feature["aa_hexdigest"]] = protein_anno
            else:
                cds_annotation_before_filtering[feature["aa_hexdigest"]] |= protein_anno


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract protein annotations from a Bakta JSON and write an index keyed by aa_hexdigest"
    )

    parser.add_argument(
        "-c",
        "--cds_annotation_pickles",
        type=Path,
        nargs="+",
        required=True,
        help="List of CDS annotation pickles with annotated pseudogenes",
    )

    parser.add_argument(
        "-i",
        "--sorf_annotation_pickles",
        type=Path,
        nargs="+",
        required=True,
        help="List of annotation pickle objects to extract sORF features from and generate auxiliary DB",
    )

    parser.add_argument(
        "-a",
        "--auxiliary_db",
        type=Path,
        required=True,
        help="Path to JSON to generate or update with new proteins",
    )

    parser.add_argument(
        "-b",
        "--bulk_annotation_before_filtering",
        type=Path,
        required=True,
        help="Path to bulk CDS annotation JSON",
    )

    parser.add_argument(
        "-o",
        "--updated_db_out",
        type=Path,
        required=True,
        help="Output path for updated auxiliary JSON DB",
    )

    args = parser.parse_args()

    bulk_annotations = ut.load_json(args.bulk_annotation_before_filtering)
    # update bulk annotations with pseudogene search results
    update_cds_annotation(args.cds_annotation_pickles, bulk_annotations)
    # update bulk annotations with sORF annotations
    update_cds_annotation(args.sorf_annotation_pickles, bulk_annotations)

    if os.path.exists(args.auxiliary_db):
        print(f"Updating existing auxiliary DB {args.auxiliary_db}")
        existing_aux_db = ut.load_pangenome(args.auxiliary_db)
        bulk_annotations = existing_aux_db | bulk_annotations
    else:
        print(f"Saving auxiliary DB annotation to {args.updated_db_out}")
    ut.dump_pangenome(bulk_annotations, args.updated_db_out)
