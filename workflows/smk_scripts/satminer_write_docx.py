#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# Ensure repo root is on sys.path when invoked as a script (Snakemake runs `python3 workflows/...`).
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts._satminer_report_io import write_docx_table


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--table-tsv", required=True)
    ap.add_argument("--out-docx", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--org", default="None")
    ap.add_argument("--genomes", default="None")
    ap.add_argument("--embed-top-images", type=int, default=10)
    args = ap.parse_args()

    df = pd.read_csv(args.table_tsv, sep="\t")
    write_docx_table(
        df=df,
        out_docx=Path(args.out_docx),
        sample=str(args.sample),
        org=str(args.org),
        genomes=str(args.genomes),
        embed_top_images=int(args.embed_top_images),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

