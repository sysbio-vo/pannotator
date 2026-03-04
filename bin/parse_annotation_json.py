#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def extract_protein_annotations(annotations_json: dict) -> dict:
    assert "features" in annotations_json

    result = dict()

    for protein in annotations_json["features"]:
        aa = protein["aa"]
        aa_hexdigest = protein["aa_hexdigest"]
        product = protein["product"]
        protein_type = protein["type"]
        protein_locus = protein["locus"]
        protein_db_xrefs = protein["db_xrefs"]
        protein_length = protein["length"]
        protein_hypothetical = protein.get("hypothetical", False)
        protein_description = protein["description"]
        protein_gene = protein["gene"]
        protein_genes = protein["genes"]

        record = {
            "aa": aa,
            "description": protein_description,
            "length": protein_length,
            "type": protein_type,
            "locus": protein_locus,
            "product": product,
            "hypothetical": protein_hypothetical,
            "gene": protein_gene,
            "genes": protein_genes,
            "db_xrefs": protein_db_xrefs,
        }

        if "pscc" in protein:
            record["pscc"] = protein["pscc"]
        if "seq_stats" in protein:
            record["seq_stats"] = protein["seq_stats"]

        result[aa_hexdigest] = record

    return result


def load_json(infile: str) -> dict:
    with open(infile, "r") as f:
        res_dict = json.load(f)
    return res_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract protein annotations from a Bakta JSON and write an index keyed by aa_hexdigest")
    parser.add_argument(
        "-a", "--bakta_annotation_json", type=Path, required=True, 
        help="Path to Bakta annotation JSON file",
        )
    parser.add_argument(
        "-o", "--proteins_json", type=Path, required=True, 
        help="Output path for the proteins index JSON (keyed by aa_hexdigest)"
        )
    args = parser.parse_args()

    protein_annotations_dict = extract_protein_annotations(load_json(args.bakta_annotation_json))
    args.proteins_json.write_text(json.dumps(protein_annotations_dict, indent=2))
