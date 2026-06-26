#!/usr/bin/env python3

import argparse
import hashlib
import json
from pathlib import Path
from typing import Dict, Iterable, Set, Tuple


# TODO: replace with an existing function from Bakta
def md5_hexdigest(seq: str) -> str:
    return hashlib.md5(seq.encode("utf-8")).hexdigest()


def iter_tsv_pairs(tsv_path: Path) -> Iterable[Tuple[str, str]]:
    """
    Iterate over a TSV file and yield (rep_locus, member_locus) pairs.
    """
    with tsv_path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep_locus, member_locus = parts[0], parts[1]
            yield rep_locus, member_locus


# TODO: can we use a standard parser instead? Like the one from the Biopython package?
def fasta_parser(fasta_path: Path, loci_to_extend: Set[str]) -> Dict[str, Tuple[str, str, int]]:
    """
    Parse a FASTA file and return mappings:
        locus: (aa_seq, aa_hexdigest, length)
    for entries in `loci_to_extend`.
    """
    result = {}
    current_id = None
    seq_chunks = []

    with fasta_path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                header = line[1:].strip()
                locus = header.split()[0] if header else ""

                if current_id is None:
                    current_id = locus
                    seq_chunks = []
                else:
                    if not seq_chunks:
                        current_id = locus
                    else:
                        if current_id in loci_to_extend:
                            aa_seq = "".join(seq_chunks)
                            aa_hexdigest = md5_hexdigest(aa_seq)
                            result[current_id] = (aa_seq, aa_hexdigest, len(aa_seq))

                        current_id = locus
                        seq_chunks = []

                continue

            if current_id is None:
                continue

            seq_chunks.append(line)

    if current_id is not None and seq_chunks and current_id in loci_to_extend:
        aa_seq = "".join(seq_chunks)
        aa_hexdigest = md5_hexdigest(aa_seq)
        result[current_id] = (aa_seq, aa_hexdigest, len(aa_seq))

    return result


def clone_rep_annotation(rep_annotation: dict) -> dict:
    """
    Clone a representative annotation dict, removing fields that will be updated for the member
    """
    new_entry = rep_annotation.copy()
    for field in ["locus", "aa", "length"]:
        new_entry.pop(field, None)
    return new_entry


def main():
    parser = argparse.ArgumentParser(description="Propagate representative annotations to all MMseqs cluster members")
    parser.add_argument("--tsv", type=Path, required=True, help="mmseqs clustering TSV")
    parser.add_argument("--fasta", type=Path, required=True, help="FASTA with all sequences")
    parser.add_argument(
        "--json-in", type=Path, required=True, help="bulk_protein_annotations.json (representatives annotated)"
    )
    parser.add_argument("--json-out", type=Path, required=True, help="output JSON with all cluster members added")
    parser.add_argument(
        "--add-trace",
        action="store_true",
        help="If set, adds trace fields: cluster_rep_locus and cluster_rep_aa_hexdigest to propagated entries",
    )
    args = parser.parse_args()

    # load rep - member pairs
    loci_to_extend = set()
    rep_loci = set()
    n_pairs = 0
    for rep, mem in iter_tsv_pairs(args.tsv):
        loci_to_extend.add(rep)
        loci_to_extend.add(mem)
        rep_loci.add(rep)
        n_pairs += 1

    # load sequences for loci to extend
    locus_to_seq = fasta_parser(args.fasta, loci_to_extend)

    # load existing annotations
    with args.json_in.open("r") as f:
        annotations = json.load(f)

    # propagate annotations from reps to members
    added_count = 0
    processed_count = 0

    for rep, mem in iter_tsv_pairs(args.tsv):
        processed_count += 1

        rep_entry = locus_to_seq.get(rep)
        mem_entry = locus_to_seq.get(mem)
        if rep_entry is None or mem_entry is None:
            continue

        rep_seq, rep_md5, rep_len = rep_entry
        mem_seq, mem_md5, mem_len = mem_entry

        rep_annotation = annotations.get(rep_md5)
        if rep_annotation is None:
            continue

        if mem_md5 in annotations:
            continue

        new_entry = clone_rep_annotation(rep_annotation)
        new_entry["locus"] = mem
        new_entry["aa"] = mem_seq
        new_entry["length"] = mem_len

        if args.add_trace:
            new_entry["cluster_rep_locus"] = rep
            new_entry["cluster_rep_aa_hexdigest"] = rep_md5

        annotations[mem_md5] = new_entry
        added_count += 1

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    with args.json_out.open("w") as f:
        json.dump(annotations, f, indent=2)

    print(f"Processed {processed_count} pairs, added {added_count} annotations")


if __name__ == "__main__":
    main()
