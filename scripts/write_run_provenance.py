#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _run(cmd: list[str]) -> dict:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        return {
            "cmd": cmd,
            "returncode": proc.returncode,
            "stdout": (proc.stdout or "").strip(),
            "stderr": (proc.stderr or "").strip(),
        }
    except Exception as e:
        return {"cmd": cmd, "error": repr(e)}


def main() -> int:
    ap = argparse.ArgumentParser(description="Write a run-level provenance JSON for RepOrtR runs.")
    ap.add_argument("--out", required=True, help="Output JSON path.")
    ap.add_argument(
        "--workdir",
        default=".",
        help="Repo/workdir path (default: current directory).",
    )
    ap.add_argument(
        "--configfile",
        default="projects/global_config.yaml",
        help="Primary config file path (default: projects/global_config.yaml).",
    )
    ap.add_argument(
        "--snakemake-cmd",
        default=None,
        help="Optional: the exact Snakemake CLI string used for the run (for audit).",
    )
    args = ap.parse_args()

    workdir = Path(args.workdir).resolve()
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    configfile = (workdir / args.configfile).resolve()
    now = datetime.now(timezone.utc).isoformat()

    env = {
        "CONDA_DEFAULT_ENV": os.environ.get("CONDA_DEFAULT_ENV", ""),
        "PATH": os.environ.get("PATH", ""),
    }

    prov = {
        "timestamp_utc": now,
        "workdir": str(workdir),
        "host": {
            "platform": platform.platform(),
            "python": sys.version.replace("\n", " "),
        },
        "config": {
            "configfile": str(configfile),
            "configfile_sha256": _sha256_file(str(configfile)) if configfile.exists() else None,
        },
        "tools": {
            "snakemake_version": _run(["snakemake", "--version"]),
            "conda_info": _run(["conda", "info"]),
            "conda_config_channel_priority": _run(["conda", "config", "--show", "channel_priority"]),
        },
        "git": {
            "rev_parse_HEAD": _run(["git", "rev-parse", "HEAD"]),
            "status_porcelain": _run(["git", "status", "--porcelain"]),
        },
        "invocation": {
            "snakemake_cmd": args.snakemake_cmd,
        },
        "env": env,
    }

    out_path.write_text(json.dumps(prov, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
