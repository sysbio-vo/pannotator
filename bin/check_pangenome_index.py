#!/usr/bin/env python3
"""Check that a pangenome index can be used as --auxiliary_db. Exits with 1 if it cannot."""

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterator

import pangenome_store as pg
import utils as ut

# what merge_annotations_into_pkl.py copies onto a CDS, and the type it needs downstream
EXPECTED_TYPES = {
    "type": (str,),
    "product": (str,),
    "db_xrefs": (list,),
    "gene": (str, type(None)),
    "genes": (list,),
    "pscc": (dict,),
    "pfams": (list,),
    "hypothetical": (bool,),
    "seq_stats": (dict,),
    "ups": (dict,),
    "ips": (dict,),
    "psc": (dict,),
}
REQUIRED = ("type", "product", "db_xrefs", "gene", "genes")
HASH = re.compile(r"[0-9a-f]{32}")


def iter_index(path: Path) -> Iterator[tuple[str, dict]]:
    if path.suffix == ".sqlite":
        conn = pg.connect(path, readonly=True)
        try:
            for aa_hash, payload in conn.execute("SELECT aa_hash, payload FROM protein"):
                yield aa_hash, json.loads(payload)
        finally:
            conn.close()
    else:
        yield from ut.load_pangenome(path).items()


def entry_problems(aa_hash: str, entry: dict, has_sequences: bool) -> Iterator[str]:
    if not HASH.fullmatch(aa_hash):
        yield "key is not an MD5 hex digest"
    if has_sequences and hashlib.md5(entry.get("aa", "").encode("utf-8")).hexdigest() != aa_hash:
        yield "key is not the MD5 of `aa`"
    for field in REQUIRED:
        if field not in entry:
            yield f"`{field}` is missing"
    for field, types in EXPECTED_TYPES.items():
        if field in entry and not isinstance(entry[field], types):
            yield f"`{field}` is {type(entry[field]).__name__}"
    if entry.get("hypothetical") and "seq_stats" not in entry:
        yield "hypothetical protein without `seq_stats`"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", type=Path)
    parser.add_argument("--reference", type=Path, help="Index known to work, used to list fields that differ")
    args = parser.parse_args()

    has_sequences = args.index.suffix != ".sqlite"  # SQLite payloads drop `aa`
    problems: Counter = Counter()
    examples: dict[str, str] = {}
    fields: Counter = Counter()
    count = 0
    for aa_hash, entry in iter_index(args.index):
        count += 1
        fields.update(entry.keys())
        for problem in entry_problems(aa_hash, entry, has_sequences):
            problems[problem] += 1
            examples.setdefault(problem, aa_hash)

    print(f"{args.index}: {count} proteins")
    if args.index.suffix == ".sqlite":
        for key, value in sorted(pg.get_meta(args.index).items()):
            print(f"  {key}: {value}")
    if args.reference:
        reference = Counter()
        for _aa_hash, entry in iter_index(args.reference):
            reference.update(entry.keys())
        print(f"  fields only in {args.reference.name}: {sorted(set(reference) - set(fields)) or 'none'}")
        print(f"  fields only in {args.index.name}: {sorted(set(fields) - set(reference)) or 'none'}")
    if fields["pseudogene"]:
        print(f"  NOTE: {fields['pseudogene']} entries carry `pseudogene`, which holds coordinates of one genome")
    for problem, n in problems.most_common():
        print(f"  PROBLEM: {problem}: {n} (e.g. {examples[problem]})")
    if count == 0:
        print("  PROBLEM: index is empty")
    if problems or count == 0:
        sys.exit(1)
    print("  OK")


if __name__ == "__main__":
    main()
