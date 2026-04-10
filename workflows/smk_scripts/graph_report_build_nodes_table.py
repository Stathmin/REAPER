#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser(description="Build nodes.tsv for repeat interaction graph report.")
    ap.add_argument("--index-tsv", required=True, help="Repeat-to-subject index TSV from queries builder")
    ap.add_argument("--out-nodes-tsv", required=True)
    args = ap.parse_args()

    idx = Path(args.index_tsv)
    out = Path(args.out_nodes_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not idx.exists() or idx.stat().st_size == 0:
        out.write_text("repeat\tsubject_type\tsubject_id\n")
        return 0

    df = pd.read_csv(idx, sep="\t")
    keep = [c for c in ("repeat", "subject_type", "subject_id") if c in df.columns]
    nodes = df[keep].drop_duplicates(subset=["repeat"]).sort_values(["subject_type", "subject_id", "repeat"])
    nodes.to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

