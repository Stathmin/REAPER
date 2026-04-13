#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _read_rich_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


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

    # Base node identity.
    base_keep = [c for c in ("repeat", "subject_type", "subject_id", "source_rich_table") if c in df.columns]
    base = df[base_keep].drop_duplicates(subset=["repeat"]).copy()

    # Optional enrichment from satMiner rich_table.tsv. We join by `repeat`.
    # `source_rich_table` points at .../satminer/report/rich_table.tsv.
    # We keep a stable, small set of fields that are useful for labels/hover.
    enrich_cols = [
        "repeat",
        "smpl",
        "org",
        "genomes",
        "iter",
        "rank",
        "TAREAN_annotation",
        "cons_len",
        "size, %",
        "pdiv_mean",
        "pdiv_median",
        "pdiv_q05",
        "pdiv_q95",
    ]

    if "source_rich_table" in base.columns:
        parts: list[pd.DataFrame] = []
        for src, reps in base.groupby("source_rich_table")["repeat"]:
            src_p = Path(str(src))
            rt = _read_rich_table(src_p)
            if rt.empty or "repeat" not in rt.columns:
                continue
            cols = [c for c in enrich_cols if c in rt.columns]
            if not cols:
                continue
            # Reduce to 1 row per repeat deterministically (keep the first as stored).
            sub = rt[cols].drop_duplicates(subset=["repeat"]).copy()
            parts.append(sub)
        if parts:
            rich = pd.concat(parts, ignore_index=True)
            base = base.merge(rich, on="repeat", how="left")

    # Normalize awkward column names for downstream JS/UI.
    if "size, %" in base.columns:
        base = base.rename(columns={"size, %": "size_pct"})
    if "size_pct" in base.columns:
        base["size_pct"] = pd.to_numeric(base["size_pct"], errors="coerce")

    sort_cols = [c for c in ("subject_type", "subject_id", "repeat") if c in base.columns]
    nodes = base.sort_values(sort_cols) if sort_cols else base

    # Prefer to keep `source_rich_table` for provenance and UI filtering if present.
    nodes.to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

