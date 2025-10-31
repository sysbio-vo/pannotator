#!/usr/bin/env python3

import json 
import argparse
from pathlib import Path

def map_annotations(cds_index_path: Path, annotations_path: Path) -> dict:
    with open(cds_index_path, 'r') as f:
        cds_index = json.load(f)
    with open(annotations_path, 'r') as f:
        annotations = json.load(f)
    
    merged = {}
    for aa_hexdigest, cds in cds_index.items():
        protein = annotations.get(aa_hexdigest)
        if protein:
            for k, v in protein.items():
                if k not in cds:
                    cds[k] = v
            merged[aa_hexdigest] = cds
        else:
            merged[aa_hexdigest] = cds
            merged[aa_hexdigest]['annotation_missing'] = True

    return merged


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge CDS index with protein annotations')
    parser.add_argument("cds_index", type=Path, help='Path to cds_index.json')
    parser.add_argument("annotations", type=Path, help='Path to bulk_protein_annotations.json')
    parser.add_argument("-o", "--output", type=Path, help='Path to write resulting JSON file')
    args = parser.parse_args()

    merged_index = map_annotations(args.cds_index, args.annotations)

    with open(args.output, 'w') as f:
        json.dump(merged_index, f, indent=2)