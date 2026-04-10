#!/usr/bin/env python3
"""
RepeatMasker-based homology report between sequences in a FASTA.

For each sequence, run RepeatMasker with a library composed of the other sequences
and record which library entries hit the query. This mirrors satMiner's intent
but uses rmblast (available via conda) instead of crossmatch.
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Ensure repo root is on sys.path when invoked as a script (Snakemake runs `python3 workflows/...`).
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts._rm_homology_core import repeatmasker_pairwise_edges


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True, help="Output TSV: query\\thit")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    repeatmasker_pairwise_edges(
        fasta=Path(args.fasta),
        out_tsv=Path(args.out),
        threads=int(args.threads),
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

