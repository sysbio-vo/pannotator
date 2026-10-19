#!/usr/bin/env python3
"""SQLite pangenome index: exact hash lookup, JSON import, optional representative search."""

import argparse
import json
import shutil
import sqlite3
import subprocess
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from urllib.parse import quote

import utils as ut

# Legacy SQLite builds cap bound variables at 999. Stay under that.
_BATCH = 500
SCHEMA_VERSION = "1"
DB_NAME = "pangenome.sqlite"
REPS_FASTA = "reps.fasta"
REPS_DB = "reps_db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS protein (
    aa_hash TEXT PRIMARY KEY,
    payload TEXT NOT NULL
);
"""

DEFAULT_META = {
    "schema_version": SCHEMA_VERSION,
    "species": "",
    "bakta_db_version": "",
    "protein_count": "0",
    "similarity_min_seq_id": "",
    "similarity_coverage": "",
}


def default_species(path: Path) -> str:
    if path.name == DB_NAME:
        return path.parent.parent.name or path.parent.name
    return path.stem


def connect(path: Path, readonly: bool = False) -> sqlite3.Connection:
    path = Path(path)
    if readonly:
        uri = "file:" + quote(path.resolve().as_posix(), safe="/") + "?mode=ro"
        return sqlite3.connect(uri, uri=True)
    return sqlite3.connect(path)


def init_db(path: Path, meta: dict | None = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    values = dict(DEFAULT_META)
    values["species"] = default_species(path)
    if meta:
        values.update({k: str(v) for k, v in meta.items() if v is not None})
    conn = connect(path)
    try:
        conn.executescript(SCHEMA)
        conn.executemany(
            "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO NOTHING",
            list(values.items()),
        )
        conn.commit()
    finally:
        conn.close()


def get_meta(path: Path) -> dict[str, str]:
    conn = connect(path, readonly=True)
    try:
        rows = conn.execute("SELECT key, value FROM meta").fetchall()
    finally:
        conn.close()
    return {key: value for key, value in rows}


def _payload(annotation: dict) -> str:
    stored = {key: value for key, value in annotation.items() if key != "aa"}
    return json.dumps(stored, separators=(",", ":"))


def upsert(path: Path, annotations: dict, meta: dict | None = None) -> int:
    path = Path(path)
    init_db(path, meta)
    rows = []
    for aa_hash, annotation in annotations.items():
        if not aa_hash or not isinstance(annotation, dict):
            continue
        rows.append((aa_hash, _payload(annotation)))
    conn = connect(path)
    try:
        conn.execute("BEGIN")
        conn.executemany(
            "INSERT INTO protein(aa_hash, payload) VALUES(?, ?) "
            "ON CONFLICT(aa_hash) DO UPDATE SET payload=excluded.payload",
            rows,
        )
        count = conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
        updates = {"protein_count": str(count)}
        if meta:
            updates.update({k: str(v) for k, v in meta.items() if v is not None})
        conn.executemany(
            "INSERT INTO meta(key, value) VALUES(?, ?) " "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            list(updates.items()),
        )
        conn.commit()
    finally:
        conn.close()
    return count


def lookup(path: Path, hashes: list[str]) -> dict[str, dict]:
    unique = [aa_hash for aa_hash in dict.fromkeys(hashes) if aa_hash]
    if not unique:
        return {}
    hits: dict[str, dict] = {}
    conn = connect(path, readonly=True)
    try:
        for start in range(0, len(unique), _BATCH):
            end = start + _BATCH
            chunk = unique[start:end]
            placeholders = ",".join("?" * len(chunk))
            sql = f"SELECT aa_hash, payload FROM protein WHERE aa_hash IN ({placeholders})"
            for aa_hash, payload in conn.execute(sql, chunk):
                hits[aa_hash] = json.loads(payload)
    finally:
        conn.close()
    return hits


def _copy_sqlite(src: Path, dst: Path) -> None:
    """Copy src onto dst. A symlink is replaced so the link target is not modified."""
    if src.resolve() == dst.resolve():
        if not src.is_symlink():
            return
        real = src.resolve()
        tmp = dst.with_name(f".{dst.name}.copying")
        shutil.copy2(real, tmp)
        dst.unlink()
        tmp.rename(dst)
        return
    shutil.copy2(src, dst)


def prepare_sqlite(src: Path | None, dst: Path) -> None:
    """Create dst, copying an existing sqlite or importing a legacy pangenome file."""
    dst = Path(dst)
    if src is None or not Path(src).exists():
        init_db(dst)
        return
    src = Path(src)
    if src.suffix == ".sqlite":
        _copy_sqlite(src, dst)
        init_db(dst)
        return
    init_db(dst)
    upsert(dst, ut.load_pangenome(src))


def import_json(json_path: Path, sqlite_path: Path, meta: dict | None = None) -> int:
    return upsert(sqlite_path, ut.load_pangenome(json_path), meta)


def hits_from_search_tsv(tsv: Path, min_seq_id: float, coverage: float) -> dict[str, tuple[str, float, float]]:
    """Best query→rep hit. MMseqs pident is 0–100; qcov is 0–1."""
    best: dict[str, tuple[str, float, float]] = {}
    if not tsv.exists():
        return best
    for line in tsv.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        query, rep, pident, qcov = parts[0], parts[1], float(parts[2]), float(parts[3])
        identity = pident / 100.0
        if identity < min_seq_id or qcov < coverage:
            continue
        previous = best.get(query)
        if previous is None or identity > previous[1]:
            best[query] = (rep, identity, qcov)
    return best


def annotations_for_search_hits(sqlite_path: Path, best: dict[str, tuple[str, float, float]]) -> dict[str, dict]:
    if not best:
        return {}
    reps = lookup(sqlite_path, [rep for rep, _identity, _coverage in best.values()])
    hits = {}
    for query, (rep, identity, qcov) in best.items():
        annotation = reps.get(rep)
        if annotation is None:
            continue
        hit = dict(annotation)
        hit["pangenome_rep_aa_hexdigest"] = rep
        hit["pangenome_identity"] = identity
        hit["pangenome_coverage"] = qcov
        hits[query] = hit
    return hits


def _search_target(sqlite_path: Path) -> Path | None:
    parent = Path(sqlite_path).parent
    if (parent / f"{REPS_DB}.dbtype").exists():
        return parent / REPS_DB
    fasta = parent / REPS_FASTA
    if fasta.exists():
        return fasta
    return None


def search_sequences(sqlite_path: Path, sequences: list[tuple[str, str]]) -> dict[str, dict]:
    """Search exact-match misses against representative sequences. Empty if similarity is off."""
    if not sequences:
        return {}
    meta = get_meta(sqlite_path)
    min_raw = meta.get("similarity_min_seq_id", "")
    if not min_raw:
        return {}
    target = _search_target(sqlite_path)
    if target is None:
        return {}
    if shutil.which("mmseqs") is None:
        raise RuntimeError("mmseqs is required for pangenome similarity search")
    min_seq_id = float(min_raw)
    coverage = float(meta.get("similarity_coverage") or 0)
    with tempfile.TemporaryDirectory() as tmp_name:
        tmp = Path(tmp_name)
        query_fasta = tmp / "query.fasta"
        with query_fasta.open("w") as handle:
            for aa_hash, sequence in sequences:
                handle.write(f">{aa_hash}\n{sequence}\n")
        result = tmp / "result.tsv"
        subprocess.run(
            [
                "mmseqs",
                "easy-search",
                str(query_fasta),
                str(target),
                str(result),
                str(tmp / "mmseqs"),
                "--min-seq-id",
                str(min_seq_id),
                "-c",
                str(coverage),
                "--format-output",
                "query,target,pident,qcov",
                "--threads",
                "1",
            ],
            check=True,
        )
        best = hits_from_search_tsv(result, min_seq_id, coverage)
    return annotations_for_search_hits(sqlite_path, best)


def set_similarity(sqlite_path: Path, reps_fasta: Path, min_seq_id: float, coverage: float) -> None:
    sqlite_path = Path(sqlite_path)
    if not sqlite_path.exists():
        raise FileNotFoundError(sqlite_path)
    dest = sqlite_path.parent / REPS_FASTA
    if Path(reps_fasta).resolve() != dest.resolve():
        shutil.copy2(reps_fasta, dest)
    upsert(
        sqlite_path,
        {},
        {"similarity_min_seq_id": str(min_seq_id), "similarity_coverage": str(coverage)},
    )
    mmseqs = shutil.which("mmseqs")
    if mmseqs:
        subprocess.run([mmseqs, "createdb", str(dest), str(sqlite_path.parent / REPS_DB)], check=True)


def post_json(base_url: str, endpoint: str, body: dict) -> dict:
    url = base_url.rstrip("/") + "/" + endpoint
    request = urllib.request.Request(
        url,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(request) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")
        raise RuntimeError(f"Pangenome request failed ({exc.code}) {url}: {detail}") from exc


def remote_lookup(base_url: str, hashes: list[str]) -> tuple[dict[str, dict], str]:
    payload = post_json(base_url, "lookup", {"hashes": list(hashes)})
    return payload.get("hits", {}), payload.get("similarity_min_seq_id", "")


def remote_search(base_url: str, sequences: list[tuple[str, str]]) -> dict[str, dict]:
    body = {"sequences": [{"hash": aa_hash, "aa": sequence} for aa_hash, sequence in sequences]}
    return post_json(base_url, "search", body).get("hits", {})


def main() -> None:
    parser = argparse.ArgumentParser(description="Build and update a SQLite pangenome index")
    sub = parser.add_subparsers(dest="cmd", required=True)

    imp = sub.add_parser("import-json", help="Import a JSON, JSON.gz, or pickle pangenome index")
    imp.add_argument("--json", type=Path, required=True)
    imp.add_argument("--sqlite", type=Path, required=True)
    imp.add_argument("--species", default="")

    sim = sub.add_parser("set-similarity", help="Attach representative FASTA and a similarity threshold")
    sim.add_argument("--sqlite", type=Path, required=True)
    sim.add_argument("--reps", type=Path, required=True, help="FASTA headers must be aa_hexdigest values")
    sim.add_argument("--min-seq-id", type=float, required=True)
    sim.add_argument("--coverage", type=float, required=True)

    args = parser.parse_args()
    if args.cmd == "import-json":
        meta = {"species": args.species} if args.species else None
        count = import_json(args.json, args.sqlite, meta)
        print(f"Imported {count} proteins into {args.sqlite}")
        return
    set_similarity(args.sqlite, args.reps, args.min_seq_id, args.coverage)
    print(f"Similarity index set on {args.sqlite}")


if __name__ == "__main__":
    main()
