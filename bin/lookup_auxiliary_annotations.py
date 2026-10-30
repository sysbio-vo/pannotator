#!/usr/bin/env python3

import argparse
from pathlib import Path

import pangenome_store as pg
import utils as ut
from bakta.constants import CDS_NOT_TRUNCATED
from bakta.utils import calc_aa_hash
from Bio import SeqIO


def _is_url(value: str) -> bool:
    return value.startswith("http://") or value.startswith("https://")


def fetch_hits(auxiliary_db: str, queries: list[tuple[str, str]], similarity: bool) -> dict:
    """Return annotations for query hashes. Does not load a remote or SQLite index into a dict."""
    hashes = [aa_hash for aa_hash, _sequence in queries]
    if _is_url(auxiliary_db):
        hits, min_seq_id = pg.remote_lookup(auxiliary_db, hashes)
        if similarity and min_seq_id:
            misses = [(aa_hash, sequence) for aa_hash, sequence in queries if aa_hash not in hits]
            if misses:
                hits.update(pg.remote_search(auxiliary_db, misses))
        return hits

    path = Path(auxiliary_db)
    if path.suffix == ".sqlite":
        hits = pg.lookup(path, hashes)
        if similarity:
            misses = [(aa_hash, sequence) for aa_hash, sequence in queries if aa_hash not in hits]
            hits.update(pg.search_sequences(path, misses))
        return hits

    database = ut.load_pangenome(path)
    return {aa_hash: database[aa_hash] for aa_hash in hashes if aa_hash in database}


def annotate_bulk_proteins_with_auxDB(fasta_path: Path, auxiliary_db: str, similarity: bool = False):
    eligible = []
    remaining_records = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        aa_hexdigest = calc_aa_hash(str(record.seq))[1]
        # TODO: take truncation information into account !!!
        # for now if a protein is not explicitly defined as not-truncated we don't fetch annotation for it
        if CDS_NOT_TRUNCATED not in record.description:
            remaining_records.append(record)
            continue
        eligible.append((aa_hexdigest, str(record.seq), record))

    hits = fetch_hits(auxiliary_db, [(aa_hash, sequence) for aa_hash, sequence, _record in eligible], similarity)

    looked_up_annotations = {}
    for aa_hash, _sequence, record in eligible:
        if aa_hash in hits:
            looked_up_annotations[aa_hash] = hits[aa_hash]
        else:
            remaining_records.append(record)

    return looked_up_annotations, remaining_records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Lookup annotations in an auxiliary database and return a sliced JSON object"
    )
    parser.add_argument(
        "--auxiliary_db",
        type=str,
        required=True,
        help="SQLite pangenome index, legacy JSON/pickle index, or hosted lookup URL "
        "(https://host/v1/<species>/<version>)",
    )

    parser.add_argument("--proteins_fa", type=Path, required=True, help="Input proteins to annotate in FASTA format")

    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output annotation JSON compatible with Pannotator bulk JSON annotation format",
    )

    parser.add_argument(
        "--remaining_proteins_filename",
        type=Path,
        required=True,
        help="Output FASTA for the proteins not found in the auxDB and to be passed to Bakta for annotation",
    )

    parser.add_argument(
        "--similarity",
        action="store_true",
        help="Search representative sequences for exact-hash misses when the index defines a similarity threshold",
    )

    args = parser.parse_args()

    looked_up_annotations, remaining_records = annotate_bulk_proteins_with_auxDB(
        args.proteins_fa, args.auxiliary_db, similarity=args.similarity
    )

    ut.dump_json(looked_up_annotations, args.out)

    with open(args.remaining_proteins_filename, "w") as output_handle:
        count = SeqIO.write(remaining_records, output_handle, "fasta")
    print(f"Saved {count} records to {args.remaining_proteins_filename} for Bakta annotation")
