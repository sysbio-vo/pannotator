#!/usr/bin/env python3
"""Merge Bakta output pickle files for a batch of genomes of any size into a single batch-level pickle file"""

import argparse
from pathlib import Path

import utils as ut


def merge_batch(assembly_ids: list[str], pkl_files: list[Path]) -> dict:
    if len(assembly_ids) != len(pkl_files):
        raise ValueError(
            f"Mismatch: {len(assembly_ids)} sample IDs vs {len(pkl_files)} pickle files"
        )

    merged = {}
    for assembly_id, pkl_path in zip(assembly_ids, pkl_files):
        merged[assembly_id] = ut.load_pickle(pkl_path)

    return merged

def main():
    p = argparse.ArgumentParser(description="Merge pickle files within batches.")

    p.add_argument(
        "--assembly_ids", 
        required=True, 
        help="IDs for assemblies within the batch, comma-separated and positionally identical to pkls"
        )
    p.add_argument(
        "--out",
        required=True,
        help="Output batch pickle file path"
        )
    p.add_argument(
        "pkl_files",    # positional arg
        nargs='+',
        help="Pickle file paths, positionally identical to supplied IDs"
        )
    
    args = p.parse_args()

    asm_ids = args.assembly_ids.split(',')
    merged = merge_batch(asm_ids, args.pkl_files)
    ut.dump_pickle(merged, Path(args.out))
    print(f"Merged {len(merged)} samples into {args.out}")


if __name__ == "__main__":
    main()
