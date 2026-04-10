#!/usr/bin/env python3
"""
Build a homology edge list from a RepeatMasker .out produced with a custom library.

Expected use:
  RepeatMasker -lib <library.fasta> <genome.fasta>

The .out file reports hits of genome sequences against repeat names (library
sequence IDs). We convert that into a TSV edge list:

  query_id<TAB>hit_id<TAB>score<TAB>pdiv<TAB>qstart<TAB>qend
"""

from __future__ import annotations

import argparse
from pathlib import Path

# Ensure repo root is on sys.path when invoked as a script (Snakemake runs `python3 workflows/...`).
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts._rm_homology_core import rmout_edges


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rm-out", required=True, help="RepeatMasker .out file")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--min-hit-len", type=int, default=0)
    args = ap.parse_args()

    rm_out = Path(args.rm_out)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = rm_out.read_text(errors="replace").splitlines()
    if any("There were no repetitive sequences detected" in ln for ln in lines[:10]):
        out.write_text("query\thit\tscore\tpdiv\tqstart\tqend\thitlen\n")
        return 0

    with out.open("w") as w:
        w.write("query\thit\tscore\tpdiv\tqstart\tqend\thitlen\n")
        for (query, hit, score, pdiv, qstart, qend, hitlen) in rmout_edges(
            rm_out=rm_out, min_hit_len=int(args.min_hit_len)
        ):
            w.write(f"{query}\t{hit}\t{score}\t{pdiv}\t{qstart}\t{qend}\t{hitlen}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

