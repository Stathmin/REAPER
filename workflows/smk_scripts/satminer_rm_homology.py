#!/usr/bin/env python3
"""
Python3 replacement for satMiner rm_homology.py.

Upstream satMiner uses RepeatMasker+crossmatch in a loop to iteratively remove
homologous sequences from a FASTA. That requires a fully configured RepeatMasker
install (crossmatch engine, libraries).

Here we implement a lightweight, deterministic homology report using BLASTN
all-vs-all on the provided FASTA. This is sufficient for "which sequences are
similar to which" exploration and downstream filtering decisions.
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Ensure repo root is on sys.path when invoked as a script (Snakemake runs `python3 workflows/...`).
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts._rm_homology_core import blast_all_vs_all_tsv


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Input FASTA")
    ap.add_argument("--out", required=True, help="Output TSV (BLAST outfmt 6)")
    ap.add_argument("--min-pident", type=float, default=80.0)
    ap.add_argument("--min-qcovhsp", type=float, default=50.0)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    blast_all_vs_all_tsv(
        fasta=Path(args.fasta),
        out_tsv=Path(args.out),
        min_pident=float(args.min_pident),
        min_qcovhsp=float(args.min_qcovhsp),
        threads=int(args.threads),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

