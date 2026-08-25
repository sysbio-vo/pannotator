#!/usr/bin/env python3
"""Batch and unbatch Bakta pickle outputs (CDS/RNA).

batch:   merge N per-sample pickles into one batch-level pickle, {sample_id: data, ...}
unbatch: split one batch-level pickle back into N per-sample pickle files,
         plus a manifest listing sample_id -> written path (for tools like
         bakta_pseudo_bulk that expect a manifest of individual files).
"""

import argparse
from pathlib import Path

import utils as ut


def batch_pickles(sample_ids: list[str], pkl_files: list[Path]) -> dict:
    if len(sample_ids) != len(pkl_files):
        raise ValueError(
            f"Mismatch: {len(sample_ids)} sample IDs vs {len(pkl_files)} pickle files"
        )
    return {sid: ut.load_pickle(path) for sid, path in zip(sample_ids, pkl_files)}


def unbatch_pickles(batch: dict, out_dir: Path, suffix: str) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    written = {}
    for sample_id, data in batch.items():
        out_path = out_dir / f"{sample_id}{suffix}"
        ut.dump_pickle(data, out_path)
        written[sample_id] = out_path
    return written


def write_manifest(paths_by_id: dict[str, Path], manifest_path: Path) -> None:
    with open(manifest_path, "w") as mf:
        mf.write("sample_id\tpath\n")
        mf.writelines(f"{sample_id}\t{path}\n" for sample_id, path in paths_by_id.items())


def cmd_batch(args):
    sample_ids = args.sample_ids.split(",")
    merged = batch_pickles(sample_ids, args.pkl_files)
    ut.dump_pickle(merged, Path(args.output))
    print(f"Batched {len(merged)} samples into {args.output}")


def cmd_unbatch(args):
    batch = ut.load_pickle(Path(args.input))
    written = unbatch_pickles(batch, Path(args.out_dir), args.suffix)
    print(f"Unbatched {len(written)} samples into {args.out_dir}")
    if args.manifest:
        write_manifest(written, Path(args.manifest))
        print(f"Manifest written to {args.manifest}")


def main():
    p = argparse.ArgumentParser(description="Batch/unbatch Bakta pickle files")
    sub = p.add_subparsers(dest="command", required=True)

    batch_p = sub.add_parser("batch", help="Merge per-sample pickles into one batch pickle")
    batch_p.add_argument("--sample-ids", required=True, help="Comma-separated sample IDs, positionally identical to pkl_files")
    batch_p.add_argument("--output", required=True, help="Output batch pickle path")
    batch_p.add_argument("pkl_files", nargs="+", type=Path, help="Per-sample pickle files, in sample-id order")
    batch_p.set_defaults(func=cmd_batch)

    unbatch_p = sub.add_parser("unbatch", help="Split a batch pickle back into per-sample files")
    unbatch_p.add_argument("--input", required=True, help="Batch-level pickle file to split")
    unbatch_p.add_argument("--out-dir", required=True, help="Directory to write per-sample pickle files")
    unbatch_p.add_argument("--suffix", default=".pkl", help="Filename suffix for written files, e.g. '.cds-only.pkl'")
    unbatch_p.add_argument("--manifest", default=None, help="Optional manifest TSV path (sample_id -> written path)")
    unbatch_p.set_defaults(func=cmd_unbatch)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
