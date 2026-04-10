#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import matplotlib.pyplot as plt
import pandas as pd

# Ensure repo root is on sys.path when invoked as a script (Snakemake runs `python3 workflows/...`).
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts._repeatmasker_out import iter_hits, parse_repeat_tags, safe_quantile

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rm-out", required=True)
    ap.add_argument("--hits-tsv", required=True)
    ap.add_argument("--divergence-tsv", required=True)
    ap.add_argument("--plots-dir", required=True)
    args = ap.parse_args()

    rm_out = Path(args.rm_out)
    hits_tsv = Path(args.hits_tsv)
    div_tsv = Path(args.divergence_tsv)
    plots_dir = Path(args.plots_dir)
    hits_tsv.parent.mkdir(parents=True, exist_ok=True)
    div_tsv.parent.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    lines = rm_out.read_text(errors="replace").splitlines()
    if any("There were no repetitive sequences detected" in ln for ln in lines[:10]):
        pd.DataFrame(columns=["repeat", "smpl", "iter", "rank", "orig", "hit_count"]).to_csv(
            hits_tsv, sep="\t", index=False
        )
        pd.DataFrame(
            columns=[
                "repeat",
                "smpl",
                "iter",
                "rank",
                "orig",
                "n",
                "pdiv_mean",
                "pdiv_median",
                "pdiv_q05",
                "pdiv_q95",
            ]
        ).to_csv(div_tsv, sep="\t", index=False)
        return 0

    rows = []
    for h in iter_hits(lines):
        tags = parse_repeat_tags(h.repeat)
        rows.append(
            {
                "score": h.score,
                "pdiv": h.pdiv,
                "query": h.query,
                "qbegin": h.qbegin,
                "qend": h.qend,
                "strand": h.strand,
                "repeat": h.repeat,
                **tags,
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        df.to_csv(hits_tsv, sep="\t", index=False)
        df.to_csv(div_tsv, sep="\t", index=False)
        return 0

    # Abundance proxy: number of RM hits per repeat (tagged consensus id)
    hits = (
        df.groupby(["repeat", "smpl", "iter", "rank", "orig"], as_index=False)
        .size()
        .rename(columns={"size": "hit_count"})
        .sort_values(["hit_count"], ascending=False)
    )
    hits.to_csv(hits_tsv, sep="\t", index=False)

    # Divergence summary per repeat
    div_rows = []
    for (rep, smpl, iter_s, rank_s, orig), sub in df.groupby(["repeat", "smpl", "iter", "rank", "orig"]):
        vals = [float(x) for x in sub["pdiv"].dropna().tolist()]
        div_rows.append(
            {
                "repeat": rep,
                "smpl": smpl,
                "iter": iter_s,
                "rank": rank_s,
                "orig": orig,
                "n": len(vals),
                "pdiv_mean": float(sum(vals) / len(vals)) if vals else float("nan"),
                "pdiv_median": safe_quantile(vals, 0.5),
                "pdiv_q05": safe_quantile(vals, 0.05),
                "pdiv_q95": safe_quantile(vals, 0.95),
            }
        )
    div = pd.DataFrame(div_rows).sort_values(["n"], ascending=False)
    div.to_csv(div_tsv, sep="\t", index=False)

    # Plots (PNG + SVG): top repeats by hit_count and pdiv distributions by iter.
    top = hits.head(25).copy()
    if not top.empty:
        plt.figure(figsize=(12, 6))
        plt.barh(range(len(top))[::-1], top["hit_count"].tolist()[::-1])
        plt.yticks(range(len(top))[::-1], top["repeat"].tolist()[::-1], fontsize=6)
        plt.xlabel("RepeatMasker hit count (proxy abundance)")
        plt.title("Top consensus entries by RM hit count")
        plt.tight_layout()
        plt.savefig(plots_dir / "top_hits.png", dpi=200)
        plt.savefig(plots_dir / "top_hits.svg")
        plt.close()

    # pdiv by iter (boxplot)
    df["iter"] = pd.to_numeric(df["iter"], errors="coerce")
    iters = sorted([int(x) for x in df["iter"].dropna().unique().tolist()])
    if iters:
        data = [df.loc[df["iter"] == i, "pdiv"].dropna().tolist() for i in iters]
        if any(len(x) > 0 for x in data):
            plt.figure(figsize=(8, 4))
            plt.boxplot(data, labels=[str(i) for i in iters], showfliers=False)
            plt.xlabel("Iteration depth (ITER)")
            plt.ylabel("perc div. (RepeatMasker pdiv)")
            plt.title("Divergence proxy by iteration")
            plt.tight_layout()
            plt.savefig(plots_dir / "pdiv_by_iter.png", dpi=200)
            plt.savefig(plots_dir / "pdiv_by_iter.svg")
            plt.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

