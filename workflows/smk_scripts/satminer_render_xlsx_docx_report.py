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

from workflows.smk_scripts._satminer_report_io import legacy_render_report_via_repeatanalyzer


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--table-tsv", required=True)
    ap.add_argument("--out-xlsx", required=True)
    ap.add_argument("--out-docx", required=True)
    ap.add_argument("--top-tsv", required=True, help="Optional top-N TSV for Word summary")
    args = ap.parse_args()

    project = str(args.project)
    sample = str(args.sample)
    table = Path(args.table_tsv)
    out_xlsx = Path(args.out_xlsx)
    out_docx = Path(args.out_docx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    out_docx.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(table, sep="\t")

    # Try to map BLAST columns into what legacy writer expects.
    # Legacy expects: best_hit, best_evalue, best_pident, coverage, coverage_type, task_type
    # Our BLAST top TSV is outfmt6 with evalue/bitscore etc.
    if "sseqid" in df.columns and "best_hit" not in df.columns:
        df["best_hit"] = df["sseqid"]
    if "evalue" in df.columns and "best_evalue" not in df.columns:
        df["best_evalue"] = df["evalue"]
    if "pident" in df.columns and "best_pident" not in df.columns:
        df["best_pident"] = df["pident"]

    # Coverage proxy: alignment length / qlen if present.
    if "coverage" not in df.columns:
        cov = None
        if "length" in df.columns and "qlen" in df.columns:
            try:
                cov = pd.to_numeric(df["length"], errors="coerce") / pd.to_numeric(df["qlen"], errors="coerce")
            except Exception:
                cov = None
        df["coverage"] = cov if cov is not None else pd.NA
    if "coverage_type" not in df.columns:
        df["coverage_type"] = pd.NA
    if "task_type" not in df.columns:
        df["task_type"] = "blastn-short"

    legacy_render_report_via_repeatanalyzer(
        project=project,
        sample=sample,
        df=df,
        out_xlsx=out_xlsx,
        out_docx=out_docx,
    )

    # Touch top TSV output path if the caller wants it as an output artifact.
    top_path = Path(args.top_tsv)
    if top_path and top_path.exists():
        pass
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

