#!/usr/bin/env python3
"""
Run RepeatMasker using a user-supplied library (-lib), without Dfam/FamDB.

Intended satMiner-style usage in RepOrtR:
  - "genome" sequences: read-derived FASTA (e.g. prepared_forRE.fasta)
  - "consensus" library: RepeatExplorer/TAREAN consensus or selected contigs FASTA

Produces a stable output layout:
  <outdir>/<basename>.out
  <outdir>/<basename>.align   (if -a enabled)
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--lib-fasta", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--align", action="store_true", help="produce .align (-a)")
    args = ap.parse_args()

    genome = Path(args.genome_fasta)
    lib = Path(args.lib_fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    conda_exe = os.environ.get("CONDA_EXE", "conda")
    have_repeatmasker = shutil.which("RepeatMasker") is not None
    if not genome.exists():
        raise FileNotFoundError(str(genome))
    if not lib.exists():
        raise FileNotFoundError(str(lib))

    genome_abs = genome.resolve()

    cmd: list[str] = []
    if have_repeatmasker:
        cmd.append("RepeatMasker")
    else:
        # Snakemake may run this script inside an ephemeral env that doesn't
        # contain RepeatMasker. Fall back to the user-managed `reportr` env.
        cmd.extend([conda_exe, "run", "-n", "reportr", "RepeatMasker"])

    cmd += [
        "-dir",
        str(outdir.resolve()),
        "-pa",
        str(int(args.threads)),
        "-nolow",
        "-no_is",
        "-engine",
        "rmblast",
        "-lib",
        str(lib.resolve()),
    ]
    if args.align:
        cmd.append("-a")
    cmd.append(str(genome_abs))

    subprocess.run(cmd, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

