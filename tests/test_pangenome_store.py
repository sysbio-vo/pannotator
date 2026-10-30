#!/usr/bin/env python3
"""SQLite pangenome index: partial lookup, JSON import, and HTTP batch query."""

import json
import sys
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

import pangenome_server as srv  # noqa: E402
import pangenome_store as pg  # noqa: E402


class PangenomeStoreTest(unittest.TestCase):
    def test_lookup_returns_hits_only_and_drops_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "ecoli.sqlite"
            pg.upsert(
                db,
                {
                    "aaa": {"aa": "MKT", "product": "kinase", "type": "cds", "db_xrefs": ["COG:1"]},
                    "bbb": {"aa": "AAA", "product": "other", "type": "cds"},
                },
            )
            hits = pg.lookup(db, ["aaa", "missing"])
            self.assertEqual(set(hits), {"aaa"})
            self.assertEqual(hits["aaa"]["product"], "kinase")
            self.assertEqual(hits["aaa"]["db_xrefs"], ["COG:1"])
            self.assertNotIn("aa", hits["aaa"])

            pg.upsert(
                db, {"aaa": {"product": "updated kinase", "type": "cds"}, "ccc": {"product": "new", "type": "cds"}}
            )
            again = pg.lookup(db, ["aaa", "bbb", "ccc"])
            self.assertEqual(again["aaa"]["product"], "updated kinase")
            self.assertEqual(again["bbb"]["product"], "other")
            self.assertEqual(again["ccc"]["product"], "new")
            self.assertEqual(pg.get_meta(db)["protein_count"], "3")

    def test_lookup_spans_batches(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "big.sqlite"
            annotations = {f"{i:032x}": {"product": str(i), "type": "cds"} for i in range(501)}
            pg.upsert(db, annotations)
            first, last = f"{0:032x}", f"{500:032x}"
            hits = pg.lookup(db, [first, last, "missing"])
            self.assertEqual(set(hits), {first, last})

    def test_import_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = Path(tmp) / "index.json"
            src.write_text(json.dumps({"abc": {"aa": "MMM", "product": "imported", "type": "cds"}}))
            db = Path(tmp) / "out.sqlite"
            count = pg.import_json(src, db, {"species": "Escherichia_coli"})
            self.assertEqual(count, 1)
            hit = pg.lookup(db, ["abc"])["abc"]
            self.assertEqual(hit["product"], "imported")
            self.assertNotIn("aa", hit)
            self.assertEqual(pg.get_meta(db)["species"], "Escherichia_coli")

    def test_prepare_replaces_symlink_before_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            real = tmp / "real.sqlite"
            pg.upsert(real, {"aaa": {"product": "kept", "type": "cds"}})
            link = tmp / "staged.sqlite"
            link.symlink_to(real)
            pg.prepare_sqlite(link, link)
            pg.upsert(link, {"bbb": {"product": "added", "type": "cds"}})
            self.assertEqual(set(pg.lookup(real, ["aaa", "bbb"])), {"aaa"})
            self.assertFalse(link.is_symlink())
            self.assertEqual(set(pg.lookup(link, ["aaa", "bbb"])), {"aaa", "bbb"})

    def test_search_tsv_threshold_and_trace(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tsv = tmp / "hits.tsv"
            # pident is percent; qcov is a fraction. Second line loses on identity.
            tsv.write_text("query\trep1\t96\t0.9\nquery\trep2\t99\t0.95\nlow\trep1\t50\t0.99\n")
            best = pg.hits_from_search_tsv(tsv, min_seq_id=0.9, coverage=0.8)
            self.assertEqual(best["query"][0], "rep2")
            self.assertNotIn("low", best)

            db = tmp / "pan.sqlite"
            pg.upsert(db, {"rep2": {"aa": "SHOULD_DROP", "product": "from rep", "type": "cds"}})
            hits = pg.annotations_for_search_hits(db, best)
            self.assertEqual(hits["query"]["product"], "from rep")
            self.assertEqual(hits["query"]["pangenome_rep_aa_hexdigest"], "rep2")
            self.assertAlmostEqual(hits["query"]["pangenome_identity"], 0.99)
            self.assertAlmostEqual(hits["query"]["pangenome_coverage"], 0.95)
            self.assertNotIn("aa", hits["query"])

    def test_similarity_off_without_threshold(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = Path(tmp) / "pan.sqlite"
            pg.upsert(db, {"rep": {"product": "x", "type": "cds"}})
            self.assertEqual(pg.search_sequences(db, [("q", "MKT")]), {})

    def test_http_catalog_and_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            db = root / "ecoli" / "v1" / pg.DB_NAME
            pg.upsert(db, {"aaa": {"aa": "MKT", "product": "kinase", "type": "cds"}})
            server = srv.make_server(root, "127.0.0.1", 0)
            port = server.server_address[1]
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            base = f"http://127.0.0.1:{port}/v1/ecoli/v1"
            try:
                with urllib.request.urlopen(f"http://127.0.0.1:{port}/v1/catalog") as response:
                    catalog = json.load(response)
                self.assertEqual(catalog["entries"][0]["species"], "ecoli")
                self.assertEqual(catalog["entries"][0]["version"], "v1")
                self.assertEqual(catalog["entries"][0]["protein_count"], "1")

                hits, min_seq_id = pg.remote_lookup(base, ["aaa", "missing"])
                self.assertEqual(set(hits), {"aaa"})
                self.assertEqual(min_seq_id, "")
                self.assertNotIn("aa", hits["aaa"])
                self.assertEqual(pg.remote_search(base, [("zzz", "MKT")]), {})

                request = urllib.request.Request(
                    base + "/lookup",
                    data=b"{",
                    headers={"Content-Type": "application/json"},
                    method="POST",
                )
                with self.assertRaises(urllib.error.HTTPError) as caught:
                    urllib.request.urlopen(request)
                self.assertEqual(caught.exception.code, 400)

                bad_hashes = urllib.request.Request(
                    base + "/lookup",
                    data=b'{"hashes":[1]}',
                    headers={"Content-Type": "application/json"},
                    method="POST",
                )
                with self.assertRaises(urllib.error.HTTPError) as caught:
                    urllib.request.urlopen(bad_hashes)
                self.assertEqual(caught.exception.code, 400)
            finally:
                server.shutdown()
                server.server_close()


if __name__ == "__main__":
    unittest.main()
