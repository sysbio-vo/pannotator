#!/usr/bin/env python3
"""Build a pangenome index from existing Bakta annotations (one JSON per genome)."""

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable, NamedTuple

import bakta.constants as bc
import pangenome_store as pg
import utils as ut
from generate_auxiliary_db import (
    MANDATORY_BAKTA_FEATURE_FIELDS,
    sample_feature_to_annotation_entry,
)

INDEXED_TYPES = (bc.FEATURE_CDS, bc.FEATURE_SORF)
JSON_SUFFIXES = (".json", ".json.gz")
OUTPUT_SUFFIXES = (".sqlite", ".json", ".json.gz", ".pkl", ".pickle")
LOW_SHARED_FRACTION = 0.5


class Sample(NamedTuple):
    name: str
    version: tuple[str, str, str]  # Bakta version, database version, database type
    size: int
    gc: float
    entries: dict[str, str]  # aa_hexdigest -> annotation entry serialized as JSON
    skipped: Counter


def skip_reason(feature: dict) -> str | None:
    """Why a protein feature cannot be reused for an identical protein of another genome, if it cannot."""
    if any(field not in feature for field in (*MANDATORY_BAKTA_FEATURE_FIELDS, "aa_hexdigest")):
        return "missing fields"
    # A pseudogene call rests on the flanking DNA and stores coordinates of its own genome. The workflow
    # copies `pseudogene` to every identical protein and never looks up truncated ones.
    if bc.PSEUDOGENE in feature or feature.get("pseudo"):
        return "pseudogene"
    if feature.get("truncated"):
        return "truncated"
    if hashlib.md5(feature["aa"].encode("utf-8")).hexdigest() != feature["aa_hexdigest"]:
        return "hash mismatch"
    return None


def read_sample(path: Path) -> Sample:
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            data = json.load(handle)
        features, version, stats = data["features"], data["version"], data["stats"]
        sample_version = (str(version["bakta"]), str(version["db"]["version"]), str(version["db"].get("type", "")))
    except (OSError, EOFError, ValueError, KeyError, TypeError) as exc:
        raise ValueError(f"{path}: not a readable Bakta annotation JSON ({exc!r})") from exc

    entries: dict[str, str] = {}
    skipped: Counter = Counter()
    for feature in features:
        if feature.get("type") not in INDEXED_TYPES:
            continue
        reason = skip_reason(feature)
        if reason:
            skipped[reason] += 1
            continue
        entry = sample_feature_to_annotation_entry(feature)
        entry["aa_hexdigest"] = feature["aa_hexdigest"]
        entry["length"] = len(feature["aa"])
        # several copies in one genome: the CDS wins over the sORF, otherwise the first copy is kept
        previous = entries.get(feature["aa_hexdigest"])
        if previous is None or (entry["type"] == bc.FEATURE_CDS and json.loads(previous)["type"] != bc.FEATURE_CDS):
            entries[feature["aa_hexdigest"]] = json.dumps(entry, sort_keys=True, separators=(",", ":"))

    return Sample(path.name.split(".")[0], sample_version, stats.get("size", 0), stats.get("gc", 0.0), entries, skipped)


def select_entry(votes: Counter) -> dict:
    """Annotation shared by most genomes. Ties are broken by content, so input order has no effect."""
    candidates = [(-count, text) for text, count in votes.items()]
    decoded = [(json.loads(text)["type"] != bc.FEATURE_CDS, count, text) for count, text in candidates]
    return json.loads(min(decoded)[2])


def build(samples: Iterable[Sample]) -> tuple[dict[str, dict], list[dict], dict]:
    votes: dict[str, Counter] = {}
    genomes: list[tuple[Sample, list[str]]] = []
    versions: set[tuple[str, str, str]] = set()
    skipped: Counter = Counter()
    for sample in samples:
        for aa_hash, text in sample.entries.items():
            votes.setdefault(aa_hash, Counter())[text] += 1
        versions.add(sample.version)
        skipped.update(sample.skipped)
        # keep hashes only, the serialized entries are released
        genomes.append((sample._replace(entries={}), list(sample.entries)))

    genomes_with = {aa_hash: sum(counter.values()) for aa_hash, counter in votes.items()}
    rows = []
    for sample, hashes in sorted(genomes, key=lambda item: item[0].name):
        shared = sum(genomes_with[aa_hash] > 1 for aa_hash in hashes)
        rows.append(
            {
                "sample": sample.name,
                "bakta_version": sample.version[0],
                "db_version": sample.version[1],
                "db_type": sample.version[2],
                "size": sample.size,
                "gc": f"{sample.gc:.4f}",
                "proteins": len(hashes),
                "skipped": sum(sample.skipped.values()),
                "shared_fraction": f"{shared / len(hashes):.4f}" if hashes else "0.0000",
            }
        )
    annotations = {aa_hash: select_entry(votes[aa_hash]) for aa_hash in sorted(votes)}
    conflicts = sum(len(counter) > 1 for counter in votes.values())
    return annotations, rows, {"versions": versions, "skipped": skipped, "conflicts": conflicts}


def expand_inputs(inputs: Iterable[Path]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        if item.is_dir():
            found = sorted(p for p in item.iterdir() if p.is_file() and p.name.endswith(JSON_SUFFIXES))
            if not found:
                raise ValueError(f"No .json or .json.gz files in {item}")
            paths.extend(found)
        elif item.is_file():
            paths.append(item)
        else:
            raise ValueError(f"Input not found: {item}")
    # a genome given twice would vote twice
    names = Counter(p.name.split(".")[0] for p in paths)
    repeated = [name for name, count in names.items() if count > 1]
    if repeated:
        raise ValueError(f"Samples given more than once: {', '.join(repeated[:10])}")
    return paths


def write_index(annotations: dict[str, dict], out: Path, meta: dict) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix == ".sqlite":
        pg.upsert(out, annotations, meta)
    else:
        ut.dump_pangenome(annotations, out)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i",
        "--inputs",
        type=Path,
        nargs="+",
        required=True,
        help="Bakta JSON files, or directories with one .json/.json.gz per genome",
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        nargs="+",
        required=True,
        help=f"Index to write, format by suffix ({', '.join(OUTPUT_SUFFIXES)}). Several may be given.",
    )
    parser.add_argument("--species", default="", help="Species name stored in a SQLite index")
    parser.add_argument("--report", type=Path, help="Per-genome TSV report")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument(
        "--allow_mixed_versions",
        action="store_true",
        help="Accept inputs annotated with different Bakta or Bakta database versions",
    )
    parser.add_argument("--force", action="store_true", help="Replace existing output files")
    args = parser.parse_args()

    for out in args.out:
        if not out.name.endswith(OUTPUT_SUFFIXES):
            parser.error(f"Unsupported output format: {out} (use {', '.join(OUTPUT_SUFFIXES)})")
        if out.exists() and not args.force:
            parser.error(f"{out} exists. Use --force to replace it.")
    try:
        paths = expand_inputs(args.inputs)
        with ProcessPoolExecutor(args.threads) as pool:
            annotations, rows, summary = build(pool.map(read_sample, paths, chunksize=4))
    except ValueError as exc:
        sys.exit(f"ERROR: {exc}")

    versions = sorted(summary["versions"])
    if len(versions) > 1 and not args.allow_mixed_versions:
        listed = "; ".join(f"Bakta {b}, database {v} {t}" for b, v, t in versions)
        sys.exit(f"ERROR: Inputs come from different versions: {listed}. Use --allow_mixed_versions to accept.")

    meta = {
        "bakta_db_version": ",".join(dict.fromkeys(f"{v} {t}".strip() for _, v, t in versions)),
        "bakta_version": ",".join(dict.fromkeys(b for b, _, _ in versions)),
        "source_genomes": len(rows),
        "species": args.species or None,
    }
    for out in args.out:
        if out.exists():
            out.unlink()
        write_index(annotations, out, meta)
        print(f"Saved {len(annotations)} proteins to {out}")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    types = Counter(entry["type"] for entry in annotations.values())
    print(f"{len(rows)} genomes, Bakta {meta['bakta_version']}, database {meta['bakta_db_version']}")
    print(", ".join(f"{count} {name}" for name, count in types.items()))
    print(f"{summary['conflicts']} proteins annotated differently across genomes, most common annotation kept")
    print("Skipped features: " + (", ".join(f"{n} {reason}" for reason, n in summary["skipped"].items()) or "none"))
    if len(rows) > 1:
        for row in rows:
            if float(row["shared_fraction"]) < LOW_SHARED_FRACTION:
                print(
                    f"WARNING: {row['sample']} shares only {float(row['shared_fraction']):.0%} of its proteins with "
                    f"other inputs (GC {row['gc']}). It may be another species.",
                    file=sys.stderr,
                )


if __name__ == "__main__":
    main()
