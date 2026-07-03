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


def process_sample(pkl_path: Path, annotations: Dict[str, Dict[str, Any]], output_path: Path) -> None:
    print(f"Processing sample: {pkl_path}")

    total = 0
    annotated = 0
    missing = []

    data = ut.load_pickle(pkl_path)  # TODO: change if more serialization formats will be added

    features = data.get("features", [])
    if not features:
        print(f"No features found in {pkl_path}")
        return

    for feature in features:
        if feature.get("type") == "cds":
            total += 1

            if find_annotation(feature, annotations):
                annotated += 1
            else:
                missing.append(feature.get("aa_hexdigest", "no_hexdigest"))

    print(f"Total CDS features: {total}")
    print(f"Annotated: {annotated}")
    print(f"Missing: {len(missing)}")

    ut.dump_pickle(data, output_path)  # TODO: change if more serialization formats will be added

    print(f"Annotated pickle saved to: {output_path}")


def process_input_folder(input_folder: Path, json_path: Path, output_folder: Path) -> None:
    annotations = load_annotations(json_path)

    pickle_files = list(input_folder.glob("*.cds-only.pkl"))  # TODO: change if more serialization formats will be added

    if not pickle_files:
        print(f"No sample pickle files found in {input_folder}")
        return

    print(f"Found {len(pickle_files)} sample pickle files to propagate annotations to.")

    output_folder.mkdir(parents=True, exist_ok=True)

    for pkl_path in pickle_files:
        new_name = pkl_path.name.replace(".cds-only.pkl", ".cds-annotated.pkl")
        output_path = output_folder / new_name

        process_sample(pkl_path, annotations, output_path)

    print(f"All {len(pickle_files)} samples processed.")


def main():
    p = argparse.ArgumentParser(description="Annotate CDS pickles with information from JSON annotations")
    p.add_argument("--pickle_folder", required=True, help="Folder containing .cds-only.pkl files")
    p.add_argument("--annotations", required=True, help="bulk_protein_annotations.json file")
    p.add_argument("--out", default="annotated_pkl", help="Output directory")
    args = p.parse_args()

    input_folder = Path(args.pickle_folder)
    json_path = Path(args.annotations)
    output_folder = Path(args.out)

    if not input_folder.exists():
        print(f"Input folder does not exist: {input_folder}")
        return

    if not json_path.exists():
        print(f"Annotations JSON file does not exist: {json_path}")
        return

    process_input_folder(input_folder, json_path, output_folder)


if __name__ == "__main__":
    main()
