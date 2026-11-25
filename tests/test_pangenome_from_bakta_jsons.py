#!/usr/bin/env python3
"""Pangenome index built from Bakta JSON annotations. Needs the Bakta environment."""

import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import pangenome_from_bakta_jsons as pb  # noqa: E402
import pangenome_store as pg  # noqa: E402


def md5(aa: str) -> str:
    return hashlib.md5(aa.encode()).hexdigest()


def feature(aa: str, **fields) -> dict:
    base = {
        "type": "cds",
        "contig": "c1",
        "start": 1,
        "stop": 9,
        "gene": None,
        "genes": [],
        "product": "kinase",
        "db_xrefs": ["SO:0001217"],
        "aa": aa,
        "aa_hexdigest": md5(aa),
        "locus": "X_1",
    }
    return base | fields


def write_sample(folder: Path, name: str, features: list[dict], db: str = "5.1", suffix: str = ".bakta.json") -> Path:
    path = folder / f"{name}{suffix}"
    data = {
        "features": features,
        "stats": {"size": 100, "gc": 0.5},
        "version": {"bakta": "1.9.4", "db": {"version": db, "type": "full"}},
    }
    opener = gzip.open if suffix.endswith(".gz") else open
    with opener(path, "wt") as handle:
        json.dump(data, handle)
    return path


def run(*args) -> subprocess.CompletedProcess:
    command = [sys.executable, str(BIN / "pangenome_from_bakta_jsons.py"), "-t", "1", *map(str, args)]
    return subprocess.run(command, capture_output=True, text=True)


class BuildTest(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.dir = Path(tmp.name)

    def build(self, *paths):
        return pb.build([pb.read_sample(p) for p in paths])

    def test_entry_has_index_fields_only(self):
        path = write_sample(self.dir, "a", [feature("MKT")])
        annotations, _rows, _summary = self.build(path)
        entry = annotations[md5("MKT")]
        self.assertEqual(entry["length"], 3)
        self.assertEqual(entry["aa_hexdigest"], md5("MKT"))
        self.assertIsNone(entry["gene"])
        for field in ("start", "stop", "contig", "locus"):
            self.assertNotIn(field, entry)

    def test_pseudogene_does_not_reach_an_intact_copy(self):
        pseudo = feature("MKT", truncated="5-prime", pseudo=True, pseudogene={"start": 5}, product="broken")
        one = write_sample(self.dir, "a", [pseudo, feature("MAA", pseudogene={"start": 5})])
        two = write_sample(self.dir, "b", [feature("MKT")])
        annotations, _rows, summary = self.build(one, two)
        self.assertEqual(set(annotations), {md5("MKT")})
        self.assertEqual(annotations[md5("MKT")]["product"], "kinase")
        self.assertNotIn("pseudogene", annotations[md5("MKT")])
        self.assertEqual(summary["skipped"]["pseudogene"], 2)

    def test_truncated_bad_hash_and_rna_are_skipped(self):
        path = write_sample(
            self.dir,
            "a",
            [
                feature("MKT", truncated="3-prime"),
                feature("MAA", aa_hexdigest="0" * 32),
                {"type": "tRNA", "product": "tRNA-Asn"},
                feature("MGG"),
            ],
        )
        annotations, _rows, summary = self.build(path)
        self.assertEqual(set(annotations), {md5("MGG")})
        self.assertEqual(dict(summary["skipped"]), {"truncated": 1, "hash mismatch": 1})

    def test_majority_wins_in_any_order(self):
        paths = [
            write_sample(self.dir, "a", [feature("MKT", gene="lysR")]),
            write_sample(self.dir, "b", [feature("MKT", gene="ybeF")]),
            write_sample(self.dir, "c", [feature("MKT", gene="ybeF")]),
        ]
        for order in (paths, paths[::-1]):
            annotations, _rows, summary = self.build(*order)
            self.assertEqual(annotations[md5("MKT")]["gene"], "ybeF")
            self.assertEqual(summary["conflicts"], 1)
        tied = [self.build(*order)[0] for order in (paths[:2], paths[1::-1])]
        self.assertEqual(tied[0], tied[1])

    def test_cds_wins_over_sorf(self):
        sorf = feature("MKT", type="sorf")
        within = write_sample(self.dir, "a", [sorf, feature("MKT")])
        others = [write_sample(self.dir, name, [sorf]) for name in ("b", "c")]
        self.assertEqual(self.build(within)[0][md5("MKT")]["type"], "cds")
        self.assertEqual(self.build(within, *others)[0][md5("MKT")]["type"], "cds")

    def test_cli_writes_all_formats_and_refuses_overwrite(self):
        write_sample(self.dir, "a", [feature("MKT"), feature("MAA")])
        write_sample(self.dir, "b", [feature("MKT")], suffix=".json.gz")
        out, sqlite = self.dir / "out" / "p.json", self.dir / "out" / "p.sqlite"
        report = self.dir / "out" / "r.tsv"
        args = ("-i", self.dir, "-o", out, sqlite, "--species", "Some_species", "--report", report)
        done = run(*args)
        self.assertEqual(done.returncode, 0, done.stderr)
        index = json.loads(out.read_text())
        self.assertEqual(set(index), {md5("MKT"), md5("MAA")})
        self.assertEqual(index[md5("MKT")]["aa"], "MKT")
        self.assertEqual(set(pg.lookup(sqlite, list(index))), set(index))
        meta = pg.get_meta(sqlite)
        self.assertEqual((meta["species"], meta["bakta_db_version"]), ("Some_species", "5.1 full"))
        self.assertEqual((meta["source_genomes"], meta["protein_count"]), ("2", "2"))
        self.assertEqual(len(report.read_text().splitlines()), 3)
        self.assertEqual(run(*args).returncode, 2)
        self.assertEqual(run(*args, "--force").returncode, 0)
        for path in (out, sqlite):
            checked = subprocess.run([sys.executable, str(BIN / "check_pangenome_index.py"), str(path)])
            self.assertEqual(checked.returncode, 0)

    def test_cli_refuses_mixed_versions_repeated_samples_and_bad_json(self):
        out = self.dir / "p.json"
        one = write_sample(self.dir, "a", [feature("MKT")])
        two = write_sample(self.dir, "b", [feature("MKT")], db="6.0")
        self.assertIn("different versions", run("-i", one, two, "-o", out).stderr)
        self.assertFalse(out.exists())
        self.assertIn("more than once", run("-i", one, one, "-o", out).stderr)
        broken = self.dir / "broken.json"
        broken.write_text('{"features": [')
        self.assertIn("broken.json", run("-i", broken, "-o", out).stderr)
        self.assertFalse(out.exists())
        self.assertEqual(run("-i", one, two, "-o", out, "--allow_mixed_versions").returncode, 0)


if __name__ == "__main__":
    unittest.main()
