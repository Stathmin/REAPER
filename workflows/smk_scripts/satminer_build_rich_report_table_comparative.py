#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import SeqIO

# Ensure repo root is on sys.path when invoked as a script (Snakemake runs `python3 workflows/...`).
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts._satminer_table_common import (
    infer_cluster_from_orig,
    normalize_blast_columns,
    parse_tag,
    read_blast_top,
    read_div,
    read_hits,
    resolve_graph_path,
)

def _count_reads_by_prefix(reads_fa: Path, prefixes: List[str]) -> Dict[str, int]:
    counts = {p: 0 for p in prefixes}
    with reads_fa.open() as f:
        for ln in f:
            if not ln.startswith(">"):
                continue
            rid = ln[1:].strip().split()[0]
            for p in prefixes:
                if rid.startswith(p):
                    counts[p] += 1
                    break
    return counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--comparative-analysis", required=True)
    ap.add_argument("--consensus-fasta", required=True)
    ap.add_argument("--reads-fasta", required=True, help="tarean1/comparative_reads.fasta (for per-prefix read counts)")
    ap.add_argument("--hits-tsv", required=True)
    ap.add_argument("--divergence-tsv", required=True)
    ap.add_argument("--blast-best-tsv", required=True)
    ap.add_argument("--legacy-blast-summary-tsv", required=False, default="", help="Optional legacy BLAST per-cluster summary TSV")
    ap.add_argument("--prefixes", required=True, help="Comma-separated sample prefixes (e.g. KA1,KA2)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-top-tsv", required=True)
    ap.add_argument("--top-n", type=int, default=30)
    args = ap.parse_args()

    project = str(args.project)
    comp = str(args.comparative_analysis)
    consensus = Path(args.consensus_fasta)
    reads_fa = Path(args.reads_fasta)
    out_tsv = Path(args.out_tsv)
    out_top = Path(args.out_top_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_top.parent.mkdir(parents=True, exist_ok=True)

    prefixes = [p.strip() for p in str(args.prefixes).split(",") if p.strip()]
    if not prefixes:
        raise ValueError("No prefixes provided")

    hits = read_hits(Path(args.hits_tsv))
    div = read_div(Path(args.divergence_tsv))
    blast_top = read_blast_top(Path(args.blast_best_tsv))
    legacy_summary_path = Path(args.legacy_blast_summary_tsv) if str(args.legacy_blast_summary_tsv or "").strip() else None

    reads_by_prefix = _count_reads_by_prefix(reads_fa, prefixes)
    total_reads = sum(reads_by_prefix.values()) or 1

    rows = []
    for rec in SeqIO.parse(str(consensus), "fasta"):
        rep = str(rec.id)
        tags = parse_tag(rep)
        cluster = infer_cluster_from_orig(tags["orig"])
        graph_path = ""
        if cluster is not None:
            clusters_dir = Path(
                f"projects/{project}/comparative/{comp}/tarean{tags['iter']}/tarean/seqclust/clustering/clusters"
            )
            graph_path = resolve_graph_path(clusters_dir, cluster)

        seq = str(rec.seq)
        seq_wrapped = seq if len(seq) <= 1000 else (seq[:1000] + "...")
        rows.append(
            {
                "repeat": rep,
                "smpl": tags["smpl"],
                "org": tags["org"],
                "genomes": tags["genomes"],
                "iter": int(tags["iter"]),
                "rank": int(tags["rank"]),
                "orig": tags["orig"],
                "cluster_num": int(cluster) if cluster is not None else pd.NA,
                "cluster": f"CL{cluster}" if cluster is not None else "",
                "iter_cluster": (f"tarean{tags['iter']}:CL{cluster}" if cluster is not None else ""),
                "cons_len": len(seq),
                "seq": seq_wrapped,
                "pic_path": graph_path,
                "TAREAN_annotation": (
                    "Putative satellites (high confidence)" if int(tags["rank"]) == 1 else "Putative satellites (low confidence)"
                ),
            }
        )

    base = pd.DataFrame(rows)

    # Join abundance: includes hit_count + hit_count_<prefix> columns.
    if not hits.empty:
        base = base.merge(hits, on=["repeat", "smpl", "iter", "rank", "orig"], how="left")
    else:
        base["hit_count"] = pd.NA
    base["hit_count"] = pd.to_numeric(base.get("hit_count"), errors="coerce")

    # Compute overall size,% from total reads.
    base["size, %"] = (base["hit_count"] / float(total_reads) * 100.0).round(2)

    # Per-prefix size,% columns if per-prefix hit counts exist.
    for p in prefixes:
        hc_col = f"hit_count_{p}"
        if hc_col in base.columns:
            denom = max(int(reads_by_prefix.get(p, 0)), 1)
            base[f"size, %_{p}"] = (pd.to_numeric(base[hc_col], errors="coerce") / float(denom) * 100.0).round(2)

    # Join divergence on repeat tag
    if not div.empty:
        keep_cols = ["repeat", "n", "pdiv_mean", "pdiv_median", "pdiv_q05", "pdiv_q95"]
        have = [c for c in keep_cols if c in div.columns]
        base = base.merge(div[have], on="repeat", how="left")
    else:
        for c in ["n", "pdiv_mean", "pdiv_median", "pdiv_q05", "pdiv_q95"]:
            base[c] = pd.NA

    # Join BLAST best by qseqid (=repeat tag).
    if not blast_top.empty and "qseqid" in blast_top.columns:
        b = blast_top.copy().rename(columns={"qseqid": "repeat"})
        base = base.merge(b, on="repeat", how="left", suffixes=("", "_blast"))

    # Normalize BLAST columns into stable names for report writers.
    base = normalize_blast_columns(base)

    # Optional: merge legacy BLAST evidence by inferred cluster key (COMP_CL####).
    if legacy_summary_path is not None and legacy_summary_path.exists():
        legacy = pd.read_csv(legacy_summary_path, sep="\t")
        if not legacy.empty and "legacy_cluster_key" in legacy.columns and "cluster_num" in base.columns:
            def _legacy_key(cluster_num):
                try:
                    if pd.isna(cluster_num):
                        return ""
                    return f"{comp}_CL{int(cluster_num):04d}"
                except Exception:
                    return ""

            b2 = base.copy()
            b2["legacy_cluster_key"] = b2["cluster_num"].apply(_legacy_key)
            base = b2.merge(legacy, on="legacy_cluster_key", how="left", suffixes=("", "_legacy"))

    # Sort by total abundance.
    base = base.sort_values(
        by=["hit_count", "iter", "rank", "repeat"],
        ascending=[False, True, True, True],
        na_position="last",
    )

    base.to_csv(out_tsv, sep="\t", index=False)
    top_n = int(args.top_n)
    top_df = base.head(top_n).copy() if top_n > 0 else base.iloc[0:0].copy()
    top_df.to_csv(out_top, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

