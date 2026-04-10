#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument(
        "--no-parse-seqids",
        action="store_true",
        help="Disable makeblastdb -parse_seqids (useful when FASTA IDs are non-unique, e.g. oligo lists).",
    )
    args = ap.parse_args()

    fasta = Path(args.fasta)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    if not fasta.exists():
        raise FileNotFoundError(str(fasta))

    conda_exe = os.environ.get("CONDA_EXE", "conda")
    makeblastdb = ["makeblastdb"]
    if shutil.which("makeblastdb") is None:
        makeblastdb = [conda_exe, "run", "-n", "reportr", "makeblastdb"]

    cmd = makeblastdb + ["-dbtype", "nucl"]
    if not bool(args.no_parse_seqids):
        cmd.append("-parse_seqids")
    cmd += ["-in", str(fasta), "-out", str(out_prefix)]
    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

