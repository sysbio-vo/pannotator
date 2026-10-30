#!/usr/bin/env python3
"""Serve species pangenome indexes. One SQLite file per species version.

Layout:
    root/<species>/<version>/pangenome.sqlite
    root/<species>/<version>/reps.fasta      # optional, for similarity search
"""

import argparse
import json
import re
import subprocess
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pangenome_store as pg

_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
_MAX_BODY = 32 * 1024 * 1024


def version_db(root: Path, species: str, version: str) -> Path | None:
    if not _ID.match(species) or not _ID.match(version):
        return None
    path = root / species / version / pg.DB_NAME
    if not path.is_file():
        return None
    return path


def catalog(root: Path) -> dict:
    entries = []
    if not root.is_dir():
        return {"entries": entries}
    for species_dir in sorted(path for path in root.iterdir() if path.is_dir() and _ID.match(path.name)):
        for version_dir in sorted(path for path in species_dir.iterdir() if path.is_dir() and _ID.match(path.name)):
            db = version_dir / pg.DB_NAME
            if not db.is_file():
                continue
            meta = pg.get_meta(db)
            entries.append(
                {
                    "species": species_dir.name,
                    "version": version_dir.name,
                    "schema_version": meta.get("schema_version", ""),
                    "bakta_db_version": meta.get("bakta_db_version", ""),
                    "protein_count": meta.get("protein_count", "0"),
                    "similarity_min_seq_id": meta.get("similarity_min_seq_id", ""),
                    "similarity_coverage": meta.get("similarity_coverage", ""),
                }
            )
    return {"entries": entries}


def make_server(root: Path, host: str, port: int) -> ThreadingHTTPServer:
    root = root.resolve()

    class Handler(BaseHTTPRequestHandler):
        def _send(self, code: int, obj: dict) -> None:
            body = json.dumps(obj, separators=(",", ":")).encode()
            self.send_response(code)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _read_json(self) -> dict | None:
            length = int(self.headers.get("Content-Length", "0"))
            if length > _MAX_BODY:
                self._send(413, {"error": "body too large"})
                return None
            try:
                return json.loads(self.rfile.read(length).decode())
            except json.JSONDecodeError:
                self._send(400, {"error": "invalid json"})
                return None

        def do_GET(self) -> None:  # noqa: N802
            if self.path.split("?", 1)[0].rstrip("/") == "/v1/catalog":
                self._send(200, catalog(root))
                return
            self._send(404, {"error": "not found"})

        def do_POST(self) -> None:  # noqa: N802
            parts = self.path.split("?", 1)[0].strip("/").split("/")
            if len(parts) != 4 or parts[0] != "v1" or parts[3] not in {"lookup", "search"}:
                self._send(404, {"error": "not found"})
                return
            db = version_db(root, parts[1], parts[2])
            if db is None:
                self._send(404, {"error": "not found"})
                return
            body = self._read_json()
            if body is None:
                return
            try:
                if parts[3] == "lookup":
                    hashes = body.get("hashes") or []
                    if not isinstance(hashes, list) or any(not isinstance(item, str) for item in hashes):
                        self._send(400, {"error": "hashes must be a list of strings"})
                        return
                    hits = pg.lookup(db, hashes)
                    meta = pg.get_meta(db)
                    self._send(
                        200,
                        {"hits": hits, "similarity_min_seq_id": meta.get("similarity_min_seq_id", "")},
                    )
                    return
                sequences = body.get("sequences") or []
                if not isinstance(sequences, list):
                    self._send(400, {"error": "sequences must be a list"})
                    return
                pairs = [
                    (item["hash"], item["aa"])
                    for item in sequences
                    if isinstance(item, dict) and isinstance(item.get("hash"), str) and isinstance(item.get("aa"), str)
                ]
                hits = pg.search_sequences(db, pairs)
            except (OSError, RuntimeError, subprocess.CalledProcessError, json.JSONDecodeError) as exc:
                self._send(500, {"error": str(exc)})
                return
            self._send(200, {"hits": hits})

    return ThreadingHTTPServer((host, port), Handler)


def main() -> None:
    parser = argparse.ArgumentParser(description="Serve SQLite pangenome indexes over HTTP")
    parser.add_argument("--root", type=Path, required=True, help="Directory of <species>/<version>/pangenome.sqlite")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8080)
    args = parser.parse_args()
    server = make_server(args.root, args.host, args.port)
    print(f"Pangenome server on http://{args.host}:{args.port} root {args.root.resolve()}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
