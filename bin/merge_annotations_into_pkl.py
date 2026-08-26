#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import utils as ut


def load_annotations(json_path: Path) -> Dict[str, Dict[str, Any]]:
    with open(json_path, "r") as f:
        return json.load(f)


def find_annotation(feature: Dict[str, Any], annotations: Dict[str, Dict[str, Any]]) -> bool:
    aa_hexdigest = feature.get("aa_hexdigest")
    if not aa_hexdigest:
        return False

    annotation = annotations.get(aa_hexdigest)
    if not annotation:
        return False

    if "gene" in annotation:
        feature["gene"] = annotation["gene"]

    if "product" in annotation:
        feature["product"] = annotation["product"]

    if "db_xrefs" in annotation:
        feature["db_xrefs"] = annotation["db_xrefs"]

    if "pfams" in annotation:
        feature["pfams"] = annotation["pfams"]

    if "pscc" in annotation:
        feature["pscc"] = annotation["pscc"]

    if "genes" in annotation:
        feature["genes"] = annotation["genes"]

    if "pseudogene" in annotation:
        feature["pseudogene"] = annotation["pseudogene"]

    # this is the bc.HYPOTHETICAL_PROTEIN_NOT_PSEUDOGENE constant
    if "hypothetical_but_not_pseudogene" in annotation:
        feature["hypothetical_but_not_pseudogene"] = annotation["hypothetical_but_not_pseudogene"]

    # if "expert" in annotation:
    #     feature["expert"] = annotation["expert"]

    if annotation.get("hypothetical", False):
        feature["hypothetical"] = True

        if "seq_stats" in annotation:
            feature["seq_stats"] = annotation["seq_stats"]

        if "ips" in annotation:
            feature["ips"] = annotation["ips"]

        if "ups" in annotation:
            feature["ups"] = annotation["ups"]

        if "psc" in annotation:
            feature["psc"] = annotation["psc"]

    else:
        if "hypothetical" in feature:
            del feature["hypothetical"]

    return True

def annotate_sample_features(sample_id: str, data: Dict[str, Any], annotations: Dict[str, Dict[str, Any]]) -> None:
    features = data.get("features", [])
    if not features:
        print(f"No features found for sample {sample_id}")
        return

    total = annotated = 0
    missing = []
    for feature in features:
        if feature.get("type") == "cds":
            total += 1
            if find_annotation(feature, annotations):
                annotated += 1
            else:
                missing.append(feature.get("aa_hexdigest", "no_hexdigest"))

    print(f"[{sample_id}] CDS features: {total}, annotated: {annotated}, missing: {len(missing)}")

def process_batch(pkl_path: Path, annotations: dict, output_path: Path) -> None:
    print(f"Processing batch: {pkl_path}")
    batch = ut.load_pickle(pkl_path)  # {sample_id: {"features": [...], ...}, ...}

    if not batch:
        print(f"No samples found in {pkl_path}")
        return

    for sample_id, data in batch.items():
        annotate_sample_features(sample_id, data, annotations)  # unchanged, mutates in place

    ut.dump_pickle(batch, output_path)
    print(f"Annotated batch pickle saved to: {output_path}")


def process_input_pickle(input_pickle: Path, json_path: Path, output_folder: Path) -> None:
    annotations = load_annotations(json_path)
    new_name = input_pickle.name.replace(".cds-only.pkl", ".cds-annotated.pkl")
    output_path = output_folder / new_name
    process_batch(input_pickle, annotations, output_path)

# def process_input_folder(input_folder: Path, json_path: Path, output_folder: Path) -> None:
#     annotations = load_annotations(json_path)
#     pickle_files = list(input_folder.glob("*.cds-only.pkl"))

#     if not pickle_files:
#         print(f"No batch pickle files found in {input_folder}")
#         return

#     print(f"Found {len(pickle_files)} batch pickle files to propagate annotations to.")
#     output_folder.mkdir(parents=True, exist_ok=True)

#     for pkl_path in pickle_files:
#         new_name = pkl_path.name.replace(".cds-only.pkl", ".cds-annotated.pkl")
#         output_path = output_folder / new_name
#         process_batch(pkl_path, annotations, output_path)

#     print(f"All {len(pickle_files)} batches processed.")

def main():
    p = argparse.ArgumentParser(description="Annotate CDS pickles (batch format) with information from JSON annotations")
    p.add_argument("--pickle_in", required=True, help="batch-level .cds-only.pkl file")
    p.add_argument("--annotations", required=True, help="bulk_protein_annotations.json file")
    p.add_argument("--pickle_out", default="annotated_pkl", help="Folder for output annotated batch pickle files")
    args = p.parse_args()

    input_pickle = Path(args.pickle_in)
    json_path = Path(args.annotations)
    output_folder = Path(args.pickle_out)

    if not input_pickle.exists():
        print(f"Input file does not exist: {input_pickle}")
        return
    if not json_path.exists():
        print(f"Annotations JSON file does not exist: {json_path}")
        return

    if not output_folder.exists():
        output_folder.mkdir(parents=True, exist_ok=True)

    process_input_pickle(input_pickle, json_path, output_folder)



if __name__ == "__main__":
    main()
