"""Preflight checks, launching and monitoring of Pannotator (Nextflow) runs.

Runs are detached from the MCP server: a client usually starts one server per session, and a run has to outlive it.
Everything needed to report on a run is therefore kept on disk:

  <state dir>/runs/<run_id>.json      record of the run
  <outdir>/.pannotator_mcp/           Nextflow launch directory, shared by runs with this outdir so `-resume` works
  <outdir>/.pannotator_mcp/<run_id>/  params.json, console.log, nextflow.log, exit_code
"""

import contextlib
import fcntl
import gzip
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid
from collections import Counter
from datetime import datetime
from pathlib import Path

from mcp.server.mcpserver.exceptions import ToolError

PIPELINE_DIR = Path(__file__).resolve().parent.parent
LAUNCH_DIR_NAME = ".pannotator_mcp"

# Executables that have to be on PATH for a given Nextflow profile
PROFILE_TOOLS = {
    "standard": ["bakta", "mmseqs"],
    "docker": ["docker"],
    "singularity": ["singularity"],
    "conda": ["conda"],
    "sanger_local": ["singularity"],
    "sanger_lsf": ["singularity", "bsub"],
    "conda_lsf": ["conda", "bsub"],
}

# main.nf compares mmseqs_args with this exact string to decide whether annotations are propagated to cluster members
MMSEQS_ARGS_TEMPLATE = "--min-seq-id {min_seq_id} -c {coverage} --alignment-mode 3"

# The pipeline puts file basenames into task scripts unquoted
SAFE_FILENAME = re.compile(r"[A-Za-z0-9_][A-Za-z0-9._#%+=@:,-]*")
SAFE_EXTENSION = re.compile(r"[A-Za-z0-9._-]*")
# Nextflow's own constraint on `-name`
RUN_NAME = re.compile(r"[a-z](?:[a-z\d]|[-_](?=[a-z\d])){0,79}")
TASK_HASH = re.compile(r"[0-9a-f]{2}/[0-9a-f]{6,30}")
CONSOLE_TASK = re.compile(r"\[([0-9a-f]{2}/[0-9a-f]{6})\] (Submitted|Re-submitted|Cached) process > (.+)")
GLOB_CHARS = set("*?[]{}\\")

# Runs Nextflow and records its exit code. The trap keeps bash alive on cancellation, so that the code is still
# written; unlike an ignored signal, it is reset in the child, so Nextflow still receives the signal.
WRAPPER = 'trap : TERM INT; "$@"; code=$?; echo "$code" > "$0.tmp"; mv "$0.tmp" "$0"'


# ---------------------------------------------------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------------------------------------------------


def state_dir() -> Path:
    if os.environ.get("PANNOTATOR_MCP_STATE_DIR"):
        return Path(os.environ["PANNOTATOR_MCP_STATE_DIR"]).expanduser().resolve()
    base = os.environ.get("XDG_STATE_HOME") or "~/.local/state"
    return Path(base).expanduser() / "pannotator-mcp"


def run_environment() -> dict[str, str]:
    """Environment for Nextflow, without what belongs to the MCP server and the agent that started it.

    - Agents leave their own variables, including session tokens, which pipeline tasks must not inherit. Nextflow also
      switches to another console log format if it finds AGENT, CLAUDECODE or NXF_AGENT_MODE, even when the latter is
      false, so they have to be unset.
    - `uv run` prepends the virtual environment of the server to PATH. With the `standard` profile the pipeline's
      Python scripts would then pick up an interpreter that lacks their dependencies.
    """
    agent_variables = {"AGENT", "AI_AGENT", "NXF_AGENT_MODE"}
    env = {
        name: value
        for name, value in os.environ.items()
        if name not in agent_variables and not name.startswith("CLAUDE")
    }
    if sys.prefix != sys.base_prefix:
        env.pop("VIRTUAL_ENV", None)
        server_bin = os.path.dirname(sys.executable)
        env["PATH"] = os.pathsep.join(p for p in env.get("PATH", "").split(os.pathsep) if p != server_bin)
    # by default images are stored in the work directory, i.e. pulled again for every outdir
    env.setdefault("NXF_SINGULARITY_CACHEDIR", str(state_dir() / "singularity"))
    env["NXF_SYNTAX_PARSER"] = "v1"  # same as nxf_bin/pannotator
    env["NXF_ANSI_LOG"] = "false"
    env["NXF_DISABLE_CHECK_LATEST"] = "true"
    return env


def find_tool(name: str) -> str | None:
    return shutil.which(name, path=run_environment().get("PATH"))


def find_nextflow() -> str | None:
    return find_tool(os.environ.get("PANNOTATOR_NEXTFLOW") or "nextflow")


def pipeline_version() -> str:
    match = re.search(r"^\s*version\s*=\s*'([^']+)'", (PIPELINE_DIR / "nextflow.config").read_text(), re.MULTILINE)
    return match.group(1) if match else "unknown"


def nextflow_version(nextflow: str) -> str | None:
    try:
        result = subprocess.run(
            [nextflow, "-v"],
            capture_output=True,
            text=True,
            timeout=60,
            stdin=subprocess.DEVNULL,
            env=run_environment(),
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.search(r"version (\S+)", result.stdout)
    return match.group(1) if match else None


def missing_tools(profile: str) -> list[str]:
    return [tool for tool in PROFILE_TOOLS[profile] if not find_tool(tool)]


def check_setup() -> dict:
    nextflow = find_nextflow()
    missing = {profile: missing_tools(profile) for profile in PROFILE_TOOLS}
    databases = []
    candidates = [default_bakta_db(None), default_bakta_db("light"), default_bakta_db("full")]
    for path in dict.fromkeys(candidates):
        if path.exists():
            databases.append({"path": str(path), "type": read_db_type(path)})
    return {
        "pannotator_version": pipeline_version(),
        "pipeline_dir": str(PIPELINE_DIR),
        "nextflow": {"path": nextflow, "version": nextflow_version(nextflow) if nextflow else None},
        "usable_profiles": [profile for profile, tools in missing.items() if not tools] if nextflow else [],
        "missing_tools_by_profile": {profile: tools for profile, tools in missing.items() if tools},
        "bakta_databases": databases,
        "cpus": os.cpu_count(),
        "total_memory_gb": total_memory_gb(),
    }


def total_memory_gb() -> float | None:
    try:
        return round(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1024**3, 1)
    except (ValueError, OSError):
        return None


# ---------------------------------------------------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------------------------------------------------


def resolve_path(value: str, what: str) -> Path:
    """Absolute paths only: the working directory of an MCP server is up to the client."""
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise ToolError(f"{what} must be an absolute path, got '{value}'")
    return path.resolve()


def sample_id_from_name(name: str) -> str:
    """Same as `sampleIdFromName` in main.nf: everything from the first dot on is dropped."""
    return re.sub(r"(\.[^.]+)+$", "", name, count=1)


def looks_like_fasta(path: Path) -> bool:
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(path, "rb") as handle:
            return handle.read(1024).lstrip().startswith(b">")
    except (OSError, EOFError):
        return False


def scan_inputs(indir: str, infile_extension: str = "") -> dict:
    """List the files the pipeline would pick up, with the problems that would make a run fail or give wrong results."""
    directory = resolve_path(indir, "indir")
    if not directory.is_dir():
        raise ToolError(f"indir is not an existing directory: {directory}")
    if GLOB_CHARS & set(str(directory)):
        raise ToolError(f"indir must not contain glob characters ({' '.join(sorted(GLOB_CHARS))}): {directory}")
    if not SAFE_EXTENSION.fullmatch(infile_extension):
        raise ToolError("infile_extension may only contain letters, digits, '.', '_' and '-'")

    # mirrors Channel.fromPath("${indir}/*${infile_extension}"): no hidden files, no directories
    files = sorted(
        path
        for path in directory.iterdir()
        if path.is_file() and path.name.endswith(infile_extension) and not path.name.startswith(".")
    )
    samples = [sample_id_from_name(path.name) for path in files]

    errors, warnings = [], []
    if not files:
        errors.append(f"No files ending in '{infile_extension}' found in {directory}")

    unsafe = [path.name for path in files if not SAFE_FILENAME.fullmatch(path.name)]
    if unsafe:
        errors.append(f"File names with spaces or shell metacharacters, rename them: {preview(unsafe)}")

    not_fasta = [path.name for path in files if not looks_like_fasta(path)]
    if not_fasta:
        errors.append(
            f"Not readable as FASTA: {preview(not_fasta)}. Restrict the input with infile_extension if these are "
            "not assemblies"
        )

    duplicated = [sample for sample, count in Counter(samples).items() if count > 1]
    if duplicated:
        errors.append(
            "Several files map to the same sample ID (file name up to the first dot), their results would be mixed "
            f"up: {preview(duplicated)}"
        )

    truncated = [
        f"{path.name} -> {sample}"
        for path, sample in zip(files, samples)
        if len(path.name.removesuffix(".gz").removeprefix(sample).split(".")) > 2
    ]
    if truncated:
        warnings.append(f"Sample IDs are cut at the first dot of the file name: {preview(truncated)}")

    return {
        "indir": str(directory),
        "n_genomes": len(files),
        "total_size_mb": round(sum(path.stat().st_size for path in files) / 1e6, 1),
        "samples": samples,
        "errors": errors,
        "warnings": warnings,
    }


def preview(items: list[str], limit: int = 10) -> str:
    text = ", ".join(items[:limit])
    return f"{text}, ... ({len(items)} in total)" if len(items) > limit else text


# ---------------------------------------------------------------------------------------------------------------------
# Bakta database
# ---------------------------------------------------------------------------------------------------------------------


def default_bakta_db(db_type: str | None) -> Path:
    if db_type is None and os.environ.get("PANNOTATOR_BAKTA_DB"):
        return Path(os.environ["PANNOTATOR_BAKTA_DB"]).expanduser().resolve()
    return PIPELINE_DIR / "bakta_db" / ("db" if db_type == "full" else "db-light")


def read_db_type(db: Path) -> str | None:
    try:
        return json.loads((db / "version.json").read_text()).get("type")
    except (OSError, ValueError, AttributeError):
        return None


def resolve_bakta_db(bakta_db: str | None, bakta_db_type: str | None, download: bool) -> dict:
    """Pick the Bakta database and return the corresponding pipeline parameters.

    Left to itself the pipeline downloads a missing database, so a typo in a path would cost up to 33 GB of traffic.
    It also always downloads into bakta_db_dir, whatever bakta_db says. Hence downloads are opt-in and only go to the
    default location.
    """
    if bakta_db:
        db = resolve_path(bakta_db, "bakta_db")
        if not db.exists():
            raise ToolError(f"Bakta database not found: {db}. Omit bakta_db to use or download the default one")
    else:
        db = default_bakta_db(None)
        if bakta_db_type and read_db_type(db) != bakta_db_type:
            db = default_bakta_db(bakta_db_type)

    if not db.exists():
        wanted = bakta_db_type or "light"
        if db != default_bakta_db(wanted):
            raise ToolError(f"PANNOTATOR_BAKTA_DB points to a database that does not exist: {db}")
        if not download:
            size = "33.1 GB" if wanted == "full" else "1.3 GB"
            raise ToolError(
                f"No '{wanted}' Bakta database at {db}. Provide the path to an existing one as bakta_db, or, if the "
                f"user agrees to download {size}, set download_bakta_db"
            )
        return {"bakta_db": str(db), "bakta_db_dir": str(db.parent), "bakta_db_type": wanted}

    db_type = read_db_type(db)
    if db_type is None:
        raise ToolError(f"Not a Bakta database, version.json is missing or has no type: {db}")
    if bakta_db_type and bakta_db_type != db_type:
        raise ToolError(f"bakta_db_type is '{bakta_db_type}', but {db} is a '{db_type}' database")
    return {"bakta_db": str(db), "bakta_db_type": db_type}


# ---------------------------------------------------------------------------------------------------------------------
# Run records
# ---------------------------------------------------------------------------------------------------------------------


def record_path(run_id: str) -> Path:
    if not RUN_NAME.fullmatch(run_id):
        raise ToolError(f"Invalid run ID: '{run_id}'")
    return state_dir() / "runs" / f"{run_id}.json"


def save_record(record: dict) -> None:
    path = record_path(record["run_id"])
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(f".{uuid.uuid4().hex}.tmp")
    temporary.write_text(json.dumps(record, indent=2))
    temporary.replace(path)


def load_record(run_id: str) -> dict:
    try:
        return json.loads(record_path(run_id).read_text())
    except FileNotFoundError:
        raise ToolError(f"Unknown run '{run_id}', see list_runs") from None


def load_records() -> list[dict]:
    records = []
    for path in (state_dir() / "runs").glob("*.json"):
        with contextlib.suppress(OSError, ValueError):
            records.append(json.loads(path.read_text()))
    return sorted(records, key=lambda record: record["started_at"], reverse=True)


def process_identity(pid: int) -> str | None:
    """Start time of a live process, to tell it from a later one with the same PID; None if it is gone or a zombie."""
    stat = Path(f"/proc/{pid}/stat")
    if Path("/proc/self/stat").exists():
        try:
            # the command name, in parentheses, may itself contain spaces and parentheses
            fields = stat.read_text().rpartition(")")[2].split()
        except OSError:
            return None
        return None if fields[0] in "ZX" else fields[19]
    result = subprocess.run(["ps", "-o", "stat=,lstart=", "-p", str(pid)], capture_output=True, text=True)
    status, _, started = result.stdout.strip().partition(" ")
    return started.strip() if status and not status.startswith("Z") else None


def run_files(record: dict) -> Path:
    return Path(record["run_dir"])


def is_alive(record: dict) -> bool:
    return process_identity(record["pid"]) == record["process_identity"]


def run_state(record: dict) -> tuple[str, int | None]:
    # liveness first: the exit code is written before the process goes away
    if is_alive(record):
        return ("cancelling" if record.get("cancel_requested") else "running"), None
    try:
        exit_code = int((run_files(record) / "exit_code").read_text())
    except (OSError, ValueError):
        return ("cancelled" if record.get("cancel_requested") else "lost"), None
    if exit_code == 0:
        return "completed", exit_code
    return ("cancelled" if record.get("cancel_requested") else "failed"), exit_code


# ---------------------------------------------------------------------------------------------------------------------
# Launching and cancelling
# ---------------------------------------------------------------------------------------------------------------------


@contextlib.contextmanager
def launch_lock(launch_dir: Path):
    """Serialises launches into an outdir, also between servers started by different clients."""
    with open(launch_dir / "launch.lock", "w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def prepare(
    indir: str,
    outdir: str,
    infile_extension: str,
    profile: str,
    run_name: str | None,
    resume: bool,
    bakta_db: str | None,
    bakta_db_type: str | None,
    download_bakta_db: bool,
    compliant: bool,
    save_intermediate: bool,
    scratch: bool,
    min_seq_id: float,
    coverage: float,
    mmseqs_command: str,
    user_proteins: str | None,
    user_hmms: str | None,
    auxiliary_db: str | None,
    extend_auxdb: bool,
    work_dir: str | None,
) -> dict:
    """Validate everything and return the record of a run that is ready to start. Nothing is written."""
    nextflow = find_nextflow()
    if not nextflow:
        raise ToolError("nextflow not found on PATH; set PANNOTATOR_NEXTFLOW to the path of the executable")
    missing = missing_tools(profile)
    if missing:
        usable = [name for name in PROFILE_TOOLS if not missing_tools(name)]
        raise ToolError(
            f"Profile '{profile}' needs {', '.join(missing)} on PATH. "
            f"Usable profiles here: {', '.join(usable) or 'none'}"
        )

    inputs = scan_inputs(indir, infile_extension)
    if inputs["errors"]:
        raise ToolError("Input problems:\n- " + "\n- ".join(inputs["errors"]))

    out = resolve_path(outdir, "outdir")
    if out == Path(inputs["indir"]) and not infile_extension:
        raise ToolError("outdir equals indir: without infile_extension the results would become inputs of the next run")

    params = {
        "indir": inputs["indir"],
        "infile_extension": infile_extension,
        "outdir": str(out),
        "compliant": compliant,
        "save_intermediate": save_intermediate,
        "scratch": scratch,
        "mmseqs_command": mmseqs_command,
        "mmseqs_args": MMSEQS_ARGS_TEMPLATE.format(min_seq_id=float(min_seq_id), coverage=float(coverage)),
        "monochrome_logs": True,
        **resolve_bakta_db(bakta_db, bakta_db_type, download_bakta_db),
    }
    for name, value in (("user_proteins", user_proteins), ("user_hmms", user_hmms)):
        if value:
            path = resolve_path(value, name)
            if not path.is_file():
                raise ToolError(f"{name} not found: {path}")
            params[name] = str(path)
    if extend_auxdb and not auxiliary_db:
        raise ToolError("extend_auxdb needs auxiliary_db")
    if auxiliary_db:
        path = resolve_path(auxiliary_db, "auxiliary_db")
        if path.exists() and not path.is_file():
            raise ToolError(f"auxiliary_db is not a file: {path}")
        if not path.parent.is_dir():
            raise ToolError(f"Directory of auxiliary_db does not exist: {path.parent}")
        params["auxiliary_db"] = str(path)
        params["extend_auxdb"] = extend_auxdb

    run_id = run_name or f"pannotator_{datetime.now():%Y%m%d_%H%M%S}_{uuid.uuid4().hex[:4]}"
    if not RUN_NAME.fullmatch(run_id):
        raise ToolError(
            "run_name must start with a lowercase letter and consist of lowercase letters, digits, and single '-' or "
            "'_' in between, 80 characters at most"
        )
    if record_path(run_id).exists():
        raise ToolError(f"Run '{run_id}' already exists, choose another run_name or omit it")

    launch_dir = out / LAUNCH_DIR_NAME
    run_dir = launch_dir / run_id
    work = resolve_path(work_dir, "work_dir") if work_dir else launch_dir / "work"
    command = [nextflow, "-log", str(run_dir / "nextflow.log"), "run", str(PIPELINE_DIR / "main.nf")]
    command += ["-profile", profile, "-name", run_id, "-params-file", str(run_dir / "params.json")]
    command += ["-work-dir", str(work), "-ansi-log", "false"]
    if resume:
        command.append("-resume")

    return {
        "run_id": run_id,
        "profile": profile,
        "outdir": str(out),
        "launch_dir": str(launch_dir),
        "run_dir": str(run_dir),
        "work_dir": str(work),
        "command": command,
        "params": params,
        "samples": inputs["samples"],
        "warnings": inputs["warnings"],
    }


def launch(record: dict) -> dict:
    launch_dir, run_dir = Path(record["launch_dir"]), run_files(record)
    try:
        launch_dir.mkdir(parents=True, exist_ok=True)
    except OSError as error:
        raise ToolError(f"Cannot create {launch_dir}: {error}") from None

    with launch_lock(launch_dir):
        for other in load_records():
            if other["outdir"] == record["outdir"] and is_alive(other):
                raise ToolError(f"Run '{other['run_id']}' is still writing to {record['outdir']}")
        try:
            run_dir.mkdir()
        except FileExistsError:
            raise ToolError(f"{run_dir} already exists, choose another run_name") from None
        (run_dir / "params.json").write_text(json.dumps(record["params"], indent=2))

        record = {**record, "started_at": datetime.now().isoformat(timespec="seconds")}
        with open(run_dir / "console.log", "wb") as console:
            process = subprocess.Popen(
                ["bash", "-c", WRAPPER, str(run_dir / "exit_code"), *record["command"]],
                cwd=launch_dir,
                env=run_environment(),
                stdin=subprocess.DEVNULL,
                stdout=console,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
        # still our unreaped child, so the PID cannot have been reused even if it has already exited
        record |= {"pid": process.pid, "process_identity": process_identity(process.pid) or "exited"}
        save_record(record)
    # reap the process when it ends, it would stay around as a zombie otherwise
    threading.Thread(target=process.wait, daemon=True).start()
    return record


def cancel(run_id: str, force: bool, wait_seconds: float = 15) -> dict:
    record = load_record(run_id)
    if not is_alive(record):
        state, exit_code = run_state(record)
        return {"run_id": run_id, "state": state, "exit_code": exit_code, "message": "The run was not active"}

    record["cancel_requested"] = True
    save_record(record)
    with contextlib.suppress(ProcessLookupError):
        # the wrapper leads its own session, hence its process group holds Nextflow and the local tasks
        os.killpg(record["pid"], signal.SIGKILL if force else signal.SIGTERM)
    deadline = time.monotonic() + wait_seconds
    while is_alive(record) and time.monotonic() < deadline:
        time.sleep(0.2)
    state, exit_code = run_state(record)
    message = "Stopped" if state != "cancelling" else "Still shutting down, check get_run_status"
    if force:
        message += ". Containers or cluster jobs of a killed run may have to be stopped by hand"
    return {"run_id": run_id, "state": state, "exit_code": exit_code, "message": message}


# ---------------------------------------------------------------------------------------------------------------------
# Monitoring
# ---------------------------------------------------------------------------------------------------------------------


def find_trace(record: dict) -> Path | None:
    """Trace files are named by start time; runs sharing an outdir are sequential, so the first one after launch."""
    started = datetime.fromisoformat(record["started_at"]).strftime("%Y-%m-%d_%H-%M-%S")
    later = [other for other in load_records() if other["outdir"] == record["outdir"]]
    next_start = min(
        (other["started_at"] for other in later if other["started_at"] > record["started_at"]),
        default=None,
    )
    end = datetime.fromisoformat(next_start).strftime("%Y-%m-%d_%H-%M-%S") if next_start else "9999"
    for path in sorted((Path(record["outdir"]) / "pipeline_info").glob("execution_trace_*.txt")):
        if started <= path.stem.removeprefix("execution_trace_") < end:
            return path
    return None


def read_trace(record: dict) -> list[dict]:
    path = find_trace(record)
    if path is None:
        return []
    header, *rows = path.read_text(errors="replace").splitlines() or [""]
    columns = header.split("\t")
    return [dict(zip(columns, row.split("\t"))) for row in rows if row.count("\t") == len(columns) - 1]


def short_process_name(task_name: str) -> str:
    return task_name.split(" (")[0].split(":")[-1]


def task_progress(record: dict, active: bool) -> tuple[dict, list[dict]]:
    """Task counts per process, in order of first submission, and the failed tasks."""
    console = read_text(run_files(record) / "console.log")
    submitted = {}  # hash -> task name
    for task_hash, _, name in CONSOLE_TASK.findall(console):
        submitted[task_hash] = name

    progress: dict[str, Counter] = {}
    failed, finished = [], set()
    for name in submitted.values():
        progress.setdefault(short_process_name(name), Counter())
    for row in read_trace(record):
        finished.add(row["hash"])
        status = row["status"].lower()
        progress.setdefault(short_process_name(row["name"]), Counter())[status] += 1
        if status in ("failed", "aborted"):
            failed.append({"task": row["name"], "hash": row["hash"], "status": status, "exit": row["exit"]})
    for task_hash, name in submitted.items():
        if task_hash not in finished:
            progress[short_process_name(name)]["submitted" if active else "unfinished"] += 1
    return {name: dict(counts) for name, counts in progress.items()}, failed


def error_excerpt(console: str, max_lines: int = 60) -> str | None:
    lines = console.splitlines()
    for index, line in enumerate(lines):
        if re.search(r"ERROR ~|^\s*(ERROR|Unknown)\b", line):
            return "\n".join(lines[index:][:max_lines])
    return "\n".join(lines[-20:]) if lines else None


def status(run_id: str) -> dict:
    record = load_record(run_id)
    state, exit_code = run_state(record)
    active = state in ("running", "cancelling")
    progress, failed = task_progress(record, active)
    out = Path(record["outdir"])
    missing = [sample for sample in record["samples"] if not (out / f"{sample}.gff3").is_file()]

    result = {
        "run_id": run_id,
        "state": state,
        "exit_code": exit_code,
        "started_at": record["started_at"],
        "profile": record["profile"],
        "outdir": record["outdir"],
        "tasks_by_process": progress,
        "n_genomes": len(record["samples"]),
        "n_genomes_annotated": len(record["samples"]) - len(missing),
    }
    if failed:
        result["failed_tasks"] = failed[:20]
    if not active:
        with contextlib.suppress(OSError):
            finished = (run_files(record) / "exit_code").stat().st_mtime
            result["finished_at"] = datetime.fromtimestamp(finished).isoformat(timespec="seconds")
        if missing:
            result["samples_without_annotation"] = preview(missing, 20)
    if state == "failed":
        result["error"] = error_excerpt(read_text(run_files(record) / "console.log"))
    if state == "lost":
        result["error"] = "The run ended without an exit code (killed or host restarted). Launch again with resume"
    return result


def summarize_runs(limit: int) -> list[dict]:
    summaries = []
    for record in load_records()[:limit]:
        state, _ = run_state(record)
        summaries.append(
            {
                "run_id": record["run_id"],
                "state": state,
                "started_at": record["started_at"],
                "indir": record["params"]["indir"],
                "outdir": record["outdir"],
                "n_genomes": len(record["samples"]),
            }
        )
    return summaries


def read_text(path: Path) -> str:
    try:
        return path.read_text(errors="replace")
    except OSError:
        return ""


def tail(path: Path, lines: int, max_bytes: int = 200_000) -> str:
    try:
        with open(path, "rb") as handle:
            handle.seek(max(0, path.stat().st_size - max_bytes))
            text = handle.read().decode(errors="replace")
    except OSError:
        return ""
    return "\n".join(text.splitlines()[-lines:])


def run_log(run_id: str, source: str, lines: int) -> dict:
    path = run_files(load_record(run_id)) / f"{source}.log"
    return {"path": str(path), "text": tail(path, lines)}


def task_log(run_id: str, task_hash: str, lines: int) -> dict:
    record = load_record(run_id)
    if not TASK_HASH.fullmatch(task_hash):
        raise ToolError("task_hash must look like 'ab/12cdef', as reported by get_run_status and get_run_log")
    prefix, rest = task_hash.split("/")
    matches = sorted((Path(record["work_dir"]) / prefix).glob(f"{rest}*"))
    if len(matches) != 1:
        raise ToolError(
            f"{'Several' if matches else 'No'} task directories match '{task_hash}' in {record['work_dir']}"
        )
    directory = matches[0]
    return {
        "task_dir": str(directory),
        "exit_code": read_text(directory / ".exitcode").strip() or None,
        "script": tail(directory / ".command.sh", lines),
        "stdout": tail(directory / ".command.out", lines),
        "stderr": tail(directory / ".command.err", lines),
    }
