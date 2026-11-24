# import subprocess
from pathlib import Path

from mcp.server.fastmcp import FastMCP

mcp = FastMCP("pannotator")
PIPELINE_DIR = Path(__file__).resolve().parent.parent


@mcp.tool()
def launch_pipeline(
    indir: str,
    outdir: str = "pannotator_results",
    infile_ext: str = ".fa.gz",
    profile: str = "standard",
    run_name: str = "mcp_run",
    resule: bool = True,
) -> str:
    pass
    """Launch the Pannotator pipeline asynchronously and return the log path."""
    # log_file = Path(outdir) / f"{run_name}.nextflow.log"

    # cmd = ["nextflow", "-log", str(log_file)]
