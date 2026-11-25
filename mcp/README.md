# Pannotator MCP server

Lets an [MCP](https://modelcontextprotocol.io) client, such as Claude Code or Claude Desktop, run Pannotator and query its results. The server runs on the machine where the pipeline runs and communicates over stdio.

## Setup

Needs everything Pannotator itself needs (see the [main README](../README.md)), and either [uv](https://docs.astral.sh/uv/), which installs the dependencies of the server on the fly, or Python >= 3.10 with `pip install "mcp>=2,<3"`.

Claude Code:

```bash
claude mcp add pannotator -- uv run --script /path/to/pannotator/mcp/server.py
```

Claude Desktop and other clients:

```json
{
  "mcpServers": {
    "pannotator": {
      "command": "uv",
      "args": ["run", "--script", "/path/to/pannotator/mcp/server.py"],
      "env": { "PANNOTATOR_BAKTA_DB": "/path/to/bakta_db" }
    }
  }
}
```

Without uv, use `python` as the command and the path of `server.py` as its only argument.

| Environment variable       | Description                                                                                      |
| -------------------------- | ------------------------------------------------------------------------------------------------ |
| `PANNOTATOR_BAKTA_DB`      | Default Bakta database. Otherwise `bakta_db/db-light` or `bakta_db/db` in the Pannotator folder. |
| `PANNOTATOR_NEXTFLOW`      | Nextflow executable, if it is not on the `PATH` the client gives to the server.                  |
| `PANNOTATOR_MCP_STATE_DIR` | Where runs are recorded. Default: `~/.local/state/pannotator-mcp`.                               |
| `NXF_SINGULARITY_CACHEDIR` | Where Singularity images are stored. Default: `singularity` in the state directory.              |

## Tools

| Tool                    | Description                                                                                  |
| ----------------------- | -------------------------------------------------------------------------------------------- |
| `check_setup`           | Usable Nextflow profiles, Bakta databases, versions and resources.                           |
| `validate_inputs`       | Check the assemblies: format, file names, colliding sample IDs.                              |
| `launch_pipeline`       | Start a run in the background. With `dry_run`, only check and show the command.              |
| `list_runs`             | Runs started by the server, also in earlier sessions.                                        |
| `get_run_status`        | State, task counts per process, annotated genomes, and the error of a failed run.            |
| `get_run_log`           | Console or Nextflow log of a run.                                                            |
| `get_task_log`          | Script, output and exit code of one task.                                                    |
| `cancel_run`            | Stop a run.                                                                                  |
| `summarize_annotations` | Contigs, length, features by type, hypothetical proteins and pseudogenes per genome.         |
| `search_annotations`    | Find genes, products, locus tags or accessions, and count the genomes with and without hits. |

The main README is available as the resource `pannotator://readme`. Results of runs started outside the server can be queried too, by `outdir`.

## Behaviour

- Runs are detached from the server and continue when the client quits.
- Logs, parameters, the Nextflow cache and, by default, the work directory of a run are in `<outdir>/.pannotator_mcp/`. Launching again with the same `outdir` resumes. Delete the folder when a collection is finished, the work directory is large.
- Only one run at a time can write to an `outdir`.
- A missing Bakta database is never downloaded unless `download_bakta_db` is set.
- MMseqs2 options are set with `min_seq_id` and `coverage`. Free-form `mmseqs_args` and `bakta_args` are deliberately not available, since the pipeline passes such strings to a shell. Input file names with spaces or shell metacharacters are rejected for the same reason.
- The server runs with the permissions of the user, and there is no network transport.

## Tests

```bash
python -m unittest tests/test_mcp_server.py
```
