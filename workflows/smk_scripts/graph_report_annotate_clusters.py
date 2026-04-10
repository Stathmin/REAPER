#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Ensure repo root is on sys.path when invoked as a script.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import yaml  # type: ignore[import-not-found]

from workflows.smk_scripts.reportr_config import resolve_value as resolve_cfg_value


def _load_config() -> dict:
    cfg_path = Path("projects/global_config.yaml")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize cluster consensus annotations (NCBI + oligo best hits).")
    ap.add_argument("--project", required=True)
    ap.add_argument("--clusters-tsv", required=True)
    ap.add_argument("--ncbi-best-tsv", required=True)
    ap.add_argument("--oligo-best-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    project = str(args.project)
    cfg = _load_config()
    enabled = bool(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "enabled"]))
    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not enabled:
        out.write_text("")
        return 0

    clusters = pd.read_csv(Path(args.clusters_tsv), sep="\t") if Path(args.clusters_tsv).stat().st_size else pd.DataFrame()
    ncbi = pd.read_csv(Path(args.ncbi_best_tsv), sep="\t") if Path(args.ncbi_best_tsv).stat().st_size else pd.DataFrame()
    oligo = pd.read_csv(Path(args.oligo_best_tsv), sep="\t") if Path(args.oligo_best_tsv).stat().st_size else pd.DataFrame()

    # ncbi/oligo best tables are per-query (qseqid). Here qseqid == cluster_id (from align script).
    keep_ncbi = [c for c in ncbi.columns if c.startswith("best_")] + ["qseqid"]
    keep_oligo = [c for c in oligo.columns if c.startswith("oligo_") or c in ("qseqid", "oligo_fitting")]
    ncbi2 = ncbi[keep_ncbi].copy() if not ncbi.empty else pd.DataFrame(columns=["qseqid"])
    oligo2 = oligo[keep_oligo].copy() if not oligo.empty else pd.DataFrame(columns=["qseqid"])

    # Cluster sizes
    if not clusters.empty:
        sizes = clusters.groupby("cluster_id", as_index=False).size().rename(columns={"size": "n_repeats"})
    else:
        sizes = pd.DataFrame(columns=["cluster_id", "n_repeats"])

    ann = sizes.rename(columns={"cluster_id": "qseqid"}).merge(ncbi2, on="qseqid", how="left").merge(oligo2, on="qseqid", how="left")
    ann = ann.rename(columns={"qseqid": "cluster_id"})
    ann.to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

