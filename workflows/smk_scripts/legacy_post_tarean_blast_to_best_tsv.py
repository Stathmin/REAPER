#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _fmt_hit(row: pd.Series) -> str:
    # Compact multi-hit evidence (safe for TSV cell).
    return (
        f"sseqid={row.get('sseqid','')}"
        f",pident={row.get('pident','')}"
        f",length={row.get('length','')}"
        f",evalue={row.get('evalue','')}"
        f",task={row.get('task','')}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-blast-results-csv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--evalue-sig", type=float, default=1e-3)
    ap.add_argument("--top-n", type=int, default=8)
    args = ap.parse_args()

    inp = Path(args.in_blast_results_csv)
    outp = Path(args.out_tsv)
    outp.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(inp)
    if df.empty:
        pd.DataFrame(
            columns=[
                "legacy_cluster_key",
                "legacy_best_sseqid",
                "legacy_best_evalue",
                "legacy_best_pident",
                "legacy_best_length",
                "legacy_best_coverage",
                "legacy_coverage_type",
                "legacy_task_type",
                "legacy_hits_topN",
            ]
        ).to_csv(outp, sep="\t", index=False)
        return 0

    # Minimal required columns.
    required = {"qseqid", "sseqid", "evalue", "pident", "length", "qlength", "task"}
    missing = [c for c in sorted(required) if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns in legacy blast_results.csv: {missing}")

    df = df.copy()
    # Normalize numeric columns.
    for c in ["evalue", "pident", "length", "qlength"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    rows = []
    for qseqid, g in df.groupby("qseqid", dropna=False):
        if not isinstance(qseqid, str) or qseqid.strip() == "":
            continue
        g = g.dropna(subset=["evalue"]).sort_values("evalue", ascending=True)
        sig = g[g["evalue"] <= float(args.evalue_sig)]
        best = (sig.iloc[0] if not sig.empty else g.iloc[0])

        cov = None
        try:
            cov = float(best["length"]) / float(best["qlength"]) if float(best["qlength"]) > 0 else None
        except Exception:
            cov = None

        if cov is None:
            cov_type = "unknown"
        elif cov >= 0.95:
            cov_type = "near-full"
        elif cov >= 0.80:
            cov_type = "composite"
        elif cov >= 0.60:
            cov_type = "partial"
        else:
            cov_type = "weak"

        task_desc = {
            "megablast": "high homology",
            "dc-megablast": "moderate homology",
            "blastn": "general homology",
        }.get(str(best.get("task") or ""), "homology")

        topN = ";".join(_fmt_hit(r) for _, r in g.head(int(args.top_n)).iterrows())

        rows.append(
            {
                "legacy_cluster_key": qseqid,
                "legacy_best_sseqid": str(best.get("sseqid") or ""),
                "legacy_best_evalue": float(best.get("evalue")) if pd.notna(best.get("evalue")) else pd.NA,
                "legacy_best_pident": float(best.get("pident")) if pd.notna(best.get("pident")) else pd.NA,
                "legacy_best_length": float(best.get("length")) if pd.notna(best.get("length")) else pd.NA,
                "legacy_best_coverage": cov if cov is not None else pd.NA,
                "legacy_coverage_type": cov_type,
                "legacy_task_type": task_desc,
                "legacy_hits_topN": topN,
            }
        )

    pd.DataFrame(rows).to_csv(outp, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

