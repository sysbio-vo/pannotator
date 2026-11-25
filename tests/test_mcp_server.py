"""Tests of the MCP server. Nextflow is replaced by a stub, nothing is annotated.

Run from the repository root, with the `mcp` package installed:
    python -m unittest tests/test_mcp_server.py
"""

import asyncio
import gzip
import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "mcp"))

import annotations as gff  # noqa: E402
import pipeline  # noqa: E402
import server  # noqa: E402

from mcp import Client, StdioServerParameters  # noqa: E402
from mcp.server.mcpserver.exceptions import ToolError  # noqa: E402

STUB = """#!/usr/bin/env bash
if [ "$1" = "-v" ]; then echo "nextflow version 99.0.0.1"; exit; fi
mode=$(cat "$(dirname "$0")/mode")
echo "[aa/111111] Submitted process > FIND_CDSS:FIND_CDS (cds_search)"
echo "[bb/222222] Submitted process > SORF_EXTRA (sample)"
[ "$mode" = sleep ] && sleep 60
outdir=$(python3 -c "import json, sys; print(json.load(open(sys.argv[1]))['outdir'])" "$PWD"/*/params.json)
mkdir -p "$outdir/pipeline_info"
trace="$outdir/pipeline_info/execution_trace_$(date +%Y-%m-%d_%H-%M-%S).txt"
printf 'task_id\\thash\\tname\\tstatus\\texit\\n' > "$trace"
printf '1\\taa/111111\\tFIND_CDSS:FIND_CDS (cds_search)\\tCOMPLETED\\t0\\n' >> "$trace"
if [ "$mode" = fail ]; then
    printf '2\\tbb/222222\\tSORF_EXTRA (sample)\\tFAILED\\t1\\n' >> "$trace"
    echo "ERROR ~ Error executing process > 'SORF_EXTRA (sample)'"
    exit 1
fi
printf '2\\tbb/222222\\tSORF_EXTRA (sample)\\tCOMPLETED\\t0\\n' >> "$trace"
echo "##gff-version 3" > "$outdir/a.gff3"
"""

GFF = """##gff-version 3
##sequence-region contig_1 1 5000
contig_1	Bakta	region	1	5000	.	+	.	ID=contig_1
contig_1	Pyrodigal	gene	1	300	.	+	.	ID=X_1_gene;locus_tag=X_1;gene=recA
contig_1	Pyrodigal	CDS	1	300	.	+	0	ID=X_1;locus_tag=X_1;product=Recombinase%2C RecA;gene=recA;Note=SO:1,UniRef:P0A7G6
contig_1	Pyrodigal	CDS	400	700	.	-	0	ID=X_2;locus_tag=X_2;product=hypothetical protein;pseudo=True
contig_1	tRNAscan-SE	tRNA	800	875	.	-	.	ID=X_3;locus_tag=X_3;product=tRNA-Arg(cct);gene=trnR
##sequence-region contig_2 1 1000
contig_2	Bakta	region	1	1000	.	+	.	ID=contig_2
##FASTA
>contig_1
ACGT
"""


class Base(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.tmp = Path(temporary.name).resolve()

        self.bin = self.tmp / "bin"
        self.bin.mkdir()
        for tool in ("nextflow", "bakta", "mmseqs"):
            (self.bin / tool).write_text(STUB)
            (self.bin / tool).chmod(0o755)
        self.mode("ok")

        self.db = self.tmp / "db-light"
        self.db.mkdir()
        (self.db / "version.json").write_text('{"type": "light"}')

        self.indir = self.tmp / "assemblies"
        self.indir.mkdir()
        (self.indir / "a.fa").write_text(">c\nACGT\n")
        self.out = self.tmp / "out"

        environment = {
            "PATH": f"{self.bin}{os.pathsep}{os.environ['PATH']}",
            "PANNOTATOR_MCP_STATE_DIR": str(self.tmp / "state"),
            "PANNOTATOR_BAKTA_DB": str(self.db),
            "PANNOTATOR_NEXTFLOW": str(self.bin / "nextflow"),
        }
        patcher = mock.patch.dict(os.environ, environment)
        patcher.start()
        self.addCleanup(patcher.stop)

    def mode(self, value):
        (self.bin / "mode").write_text(value)

    def launch(self, **kwargs):
        return server.launch_pipeline(indir=str(self.indir), outdir=str(self.out), **kwargs)

    def wait(self, run_id):
        deadline = time.monotonic() + 20
        while (result := pipeline.status(run_id))["state"] == "running":
            self.assertLess(time.monotonic(), deadline, "run did not finish")
            time.sleep(0.05)
        return result


class InputTest(Base):
    def test_sample_ids_match_main_nf(self):
        for name, expected in (("a.fa.gz", "a"), ("GCF_1.2_x.fna", "GCF_1"), ("a", "a"), ("1_2#3.fa", "1_2#3")):
            self.assertEqual(pipeline.sample_id_from_name(name), expected)

    def test_selection(self):
        (self.indir / "b.fa.gz").write_bytes(gzip.compress(b">c\nACGT\n"))
        (self.indir / ".hidden.fa").write_text("x")
        (self.indir / "directory.fa").mkdir()
        self.assertEqual(pipeline.scan_inputs(str(self.indir))["samples"], ["a", "b"])
        self.assertEqual(pipeline.scan_inputs(str(self.indir), ".fa.gz")["samples"], ["b"])

    def test_problems(self):
        for name, content, expected in (
            ("x;rm -rf y.fa", ">c\nA\n", "metacharacters"),
            ("$(id).fa", ">c\nA\n", "metacharacters"),
            ("notes.txt", "hello", "FASTA"),
            ("empty.fa", "", "FASTA"),
            ("broken.fa.gz", "not gzip", "FASTA"),
            ("a.2.fa", ">c\nA\n", "same sample ID"),
        ):
            with self.subTest(name):
                path = self.indir / name
                path.write_text(content)
                errors = pipeline.scan_inputs(str(self.indir))["errors"]
                self.assertEqual(len(errors), 1, errors)
                self.assertIn(expected, errors[0])
                with self.assertRaises(ToolError):
                    self.launch()
                path.unlink()

    def test_truncated_sample_id_is_a_warning(self):
        (self.indir / "GCF_000005845.2.fna.gz").write_bytes(gzip.compress(b">c\nACGT\n"))
        result = pipeline.scan_inputs(str(self.indir))
        self.assertEqual(result["errors"], [])
        self.assertIn("GCF_000005845.2.fna.gz -> GCF_000005845", result["warnings"][0])
        self.assertNotIn("a.fa", result["warnings"][0])

    def test_relative_and_glob_paths_are_rejected(self):
        for path in ("assemblies", str(self.tmp / "missing")):
            with self.assertRaises(ToolError):
                pipeline.scan_inputs(path)
        globbed = self.tmp / "x[1]"
        globbed.mkdir()
        with self.assertRaisesRegex(ToolError, "glob"):
            pipeline.scan_inputs(str(globbed))


class ParameterTest(Base):
    def params(self, **kwargs):
        return self.launch(dry_run=True, **kwargs)["params"]

    def test_dry_run_writes_nothing(self):
        self.assertFalse(self.launch(dry_run=True)["started"])
        self.assertFalse(self.out.exists())
        self.assertFalse((self.tmp / "state").exists())

    def test_default_clustering_is_recognised_by_main_nf(self):
        expected = self.params(min_seq_id=1, coverage=1)["mmseqs_args"]
        self.assertIn(f"!= '{expected}'", (REPO / "main.nf").read_text())
        self.assertIn(f'mmseqs_args = "{expected}"', (REPO / "nextflow.config").read_text())

    def test_profiles_exist(self):
        configs = (REPO / "nextflow.config").read_text() + (REPO / "configs/pam_generic_base.config").read_text()
        for profile in pipeline.PROFILE_TOOLS:
            self.assertRegex(configs, rf"(?m)^    {profile} \{{")

    def test_missing_database_is_not_downloaded_silently(self):
        with self.assertRaisesRegex(ToolError, "not found"):
            self.params(bakta_db=str(self.tmp / "typo"), download_bakta_db=True)
        with mock.patch.object(pipeline, "PIPELINE_DIR", self.tmp / "nowhere"):
            with self.assertRaisesRegex(ToolError, "33.1 GB"):
                self.params(bakta_db_type="full")
            params = self.params(bakta_db_type="full", download_bakta_db=True)
        self.assertEqual(params["bakta_db"], str(self.tmp / "nowhere/bakta_db/db"))
        self.assertEqual(params["bakta_db_dir"], str(self.tmp / "nowhere/bakta_db"))

    def test_database_checks(self):
        self.assertEqual(self.params()["bakta_db"], str(self.db))
        with self.assertRaisesRegex(ToolError, "'light' database"):
            self.params(bakta_db=str(self.db), bakta_db_type="full")
        with self.assertRaisesRegex(ToolError, "Not a Bakta database"):
            self.params(bakta_db=str(self.indir))

    def test_other_checks(self):
        for kwargs, expected in (
            ({"user_hmms": str(self.tmp / "missing.hmm")}, "user_hmms not found"),
            ({"extend_auxdb": True}, "needs auxiliary_db"),
            ({"auxiliary_db": str(self.tmp / "missing/index.json")}, "does not exist"),
            ({"run_name": "Bad Name"}, "run_name"),
        ):
            with self.subTest(expected), self.assertRaisesRegex(ToolError, expected):
                self.params(**kwargs)
        (self.bin / "mmseqs").unlink()
        with mock.patch.object(pipeline, "find_tool", lambda name: None if name == "mmseqs" else name):
            with self.assertRaisesRegex(ToolError, "mmseqs"):
                self.params()

    def test_environment_hides_the_agent(self):
        variables = {"CLAUDECODE": "1", "CLAUDE_CODE_TOKEN": "x", "AGENT": "1", "AI_AGENT": "x", "NXF_AGENT_MODE": "1"}
        with mock.patch.dict(os.environ, variables):
            environment = pipeline.run_environment()
        self.assertFalse(set(variables) & set(environment))

    def test_singularity_images_are_shared_between_runs(self):
        with mock.patch.dict(os.environ):
            os.environ.pop("NXF_SINGULARITY_CACHEDIR", None)
            self.assertEqual(
                pipeline.run_environment()["NXF_SINGULARITY_CACHEDIR"], str(self.tmp / "state/singularity")
            )
            os.environ["NXF_SINGULARITY_CACHEDIR"] = "/images"
            self.assertEqual(pipeline.run_environment()["NXF_SINGULARITY_CACHEDIR"], "/images")

    def test_environment_hides_the_server_venv(self):
        venv = self.tmp / "venv"
        path = f"{venv}/bin{os.pathsep}/usr/bin"
        with (
            mock.patch.dict(os.environ, {"PATH": path, "VIRTUAL_ENV": str(venv)}),
            mock.patch.object(sys, "prefix", str(venv)),
            mock.patch.object(sys, "executable", f"{venv}/bin/python"),
        ):
            environment = pipeline.run_environment()
        self.assertEqual(environment["PATH"], "/usr/bin")
        self.assertNotIn("VIRTUAL_ENV", environment)
        self.assertEqual(environment["NXF_SYNTAX_PARSER"], "v1")


class RunTest(Base):
    def test_completed(self):
        run = self.launch(min_seq_id=0.9, coverage=0.8)
        result = self.wait(run["run_id"])
        self.assertEqual((result["state"], result["exit_code"]), ("completed", 0))
        self.assertEqual((result["n_genomes"], result["n_genomes_annotated"]), (1, 1))
        self.assertEqual(result["tasks_by_process"], {"FIND_CDS": {"completed": 1}, "SORF_EXTRA": {"completed": 1}})
        params = json.loads((self.out / ".pannotator_mcp" / run["run_id"] / "params.json").read_text())
        self.assertEqual(params["mmseqs_args"], "--min-seq-id 0.9 -c 0.8 --alignment-mode 3")
        self.assertEqual(server.list_runs()["runs"][0]["state"], "completed")
        self.assertIn("Submitted", pipeline.run_log(run["run_id"], "console", 10)["text"])

    def test_failed(self):
        self.mode("fail")
        result = self.wait(self.launch()["run_id"])
        self.assertEqual((result["state"], result["exit_code"]), ("failed", 1))
        self.assertIn("Error executing process", result["error"])
        self.assertEqual(result["failed_tasks"][0]["hash"], "bb/222222")
        self.assertEqual(result["samples_without_annotation"], "a")

    def test_cancel_and_concurrent_launch(self):
        self.mode("sleep")
        run_id = self.launch()["run_id"]
        with self.assertRaisesRegex(ToolError, "still writing"):
            self.launch()
        result = pipeline.cancel(run_id, force=False)
        self.assertEqual(result["state"], "cancelled")
        self.assertIsNotNone(result["exit_code"], "the wrapper should record the exit code of a cancelled run")
        self.assertEqual(pipeline.cancel(run_id, force=False)["message"], "The run was not active")
        self.mode("ok")
        self.assertTrue(self.launch()["started"])

    def test_killed_wrapper_is_reported_as_lost(self):
        self.mode("sleep")
        run_id = self.launch()["run_id"]
        os.killpg(pipeline.load_record(run_id)["pid"], 9)
        self.assertEqual(self.wait(run_id)["state"], "lost")

    def test_reused_pid_is_not_alive(self):
        record = {"pid": os.getpid(), "process_identity": "other"}
        self.assertFalse(pipeline.is_alive(record))
        record["process_identity"] = pipeline.process_identity(os.getpid())
        self.assertTrue(pipeline.is_alive(record))

    def test_task_log(self):
        run_id = self.launch()["run_id"]
        self.wait(run_id)
        task = self.out / ".pannotator_mcp/work/bb/222222abcdef"
        task.mkdir(parents=True)
        (task / ".command.err").write_text("boom")
        (task / ".exitcode").write_text("1")
        result = pipeline.task_log(run_id, "bb/222222", 10)
        self.assertEqual((result["stderr"], result["exit_code"]), ("boom", "1"))
        for task_hash in ("../../etc", "cc/333333"):
            with self.assertRaises(ToolError):
                pipeline.task_log(run_id, task_hash, 10)

    def test_unknown_run(self):
        for run_id in ("nope", "../nope"):
            with self.assertRaises(ToolError):
                pipeline.status(run_id)


class AnnotationTest(Base):
    def setUp(self):
        super().setUp()
        self.out.mkdir()
        (self.out / "a.gff3").write_text(GFF)
        (self.out / "b.gff3").write_text(GFF.replace("recA", "xyzA").replace("RecA", "XyzA"))

    def search(self, query, field="any", exact=False, feature_type=None, samples=None, limit=50):
        return gff.search(self.out, query, field, exact, feature_type, samples, limit)

    def test_summary(self):
        result = server.summarize_annotations(outdir=str(self.out), limit=1)
        self.assertEqual(result["n_genomes"], 2)
        self.assertEqual(
            result["genomes"],
            [
                {
                    "sample": "a",
                    "contigs": 2,
                    "length_bp": 6000,
                    "features": {"CDS": 2, "tRNA": 1},
                    "hypothetical_cds": 1,
                    "pseudogenes": 1,
                }
            ],
        )

    def test_search(self):
        result = self.search("RECA", field="gene")
        self.assertEqual((result["n_hits"], result["genomes_without_hits"]), (1, ["b"]))
        self.assertEqual(result["hits"][0]["product"], "Recombinase, RecA")
        self.assertEqual(self.search("recombinase, r")["n_hits"], 1)  # comma is encoded in the file
        self.assertEqual(self.search("rec", field="gene", exact=True)["n_hits"], 0)
        self.assertEqual(self.search("P0A7G6", field="dbxref", exact=True)["n_genomes_with_hits"], 2)
        self.assertEqual(self.search("X_", feature_type="tRNA", samples=["b"])["n_hits"], 1)
        limited = self.search("X_", limit=1)
        self.assertEqual((len(limited["hits"]), limited["n_hits"], limited["truncated"]), (1, 6, True))

    def test_errors(self):
        for call in (
            lambda: self.search(" "),
            lambda: self.search("x", samples=["missing"]),
            lambda: server.search_annotations("x"),
            lambda: server.search_annotations("x", run_id="a", outdir="/b"),
            lambda: server.summarize_annotations(outdir=str(self.indir)),
        ):
            with self.assertRaises(ToolError):
                call()


class ProtocolTest(Base):
    TOOLS = {
        "check_setup",
        "validate_inputs",
        "launch_pipeline",
        "list_runs",
        "get_run_status",
        "get_run_log",
        "get_task_log",
        "cancel_run",
        "summarize_annotations",
        "search_annotations",
    }

    def test_stdio(self):
        async def session():
            parameters = StdioServerParameters(
                command=sys.executable, args=[str(REPO / "mcp/server.py")], env=dict(os.environ), cwd=str(REPO)
            )
            async with Client(parameters) as client:
                tools = {tool.name: tool for tool in (await client.list_tools()).tools}
                setup = await client.call_tool("check_setup", {})
                started = await client.call_tool("launch_pipeline", {"indir": str(self.indir), "outdir": str(self.out)})
                invalid = await client.call_tool(
                    "launch_pipeline", {"indir": str(self.indir), "outdir": str(self.out), "min_seq_id": 2}
                )
                refused = await client.call_tool("launch_pipeline", {"indir": "relative", "outdir": str(self.out)})
                readme = await client.read_resource("pannotator://readme")
            return tools, setup, started, invalid, refused, readme

        tools, setup, started, invalid, refused, readme = asyncio.run(session())
        self.assertEqual(set(tools), self.TOOLS)
        self.assertEqual(set(tools["launch_pipeline"].input_schema["required"]), {"indir", "outdir"})
        self.assertEqual(setup.structured_content["nextflow"]["version"], "99.0.0.1")
        self.assertTrue(invalid.is_error)
        self.assertTrue(refused.is_error)
        self.assertIn("absolute path", refused.content[0].text)
        self.assertIn("# Pannotator", readme.contents[0].text)
        # the server is gone, the run has to have survived it
        self.assertEqual(self.wait(started.structured_content["run_id"])["state"], "completed")


if __name__ == "__main__":
    unittest.main()
