#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# Ensure repo root is on sys.path when invoked as a script.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import yaml  # type: ignore[import-not-found]

from workflows.smk_scripts._blast_grouping import group_hsps_to_pairs
from workflows.smk_scripts.reportr_config import resolve_value as resolve_cfg_value


def _load_config() -> dict:
    cfg_path = Path("projects/global_config.yaml")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def _blastn_cmd() -> list[str]:
    if shutil.which("blastn"):
        return ["blastn"]
    return ["conda", "run", "-n", "reportr", "blastn"]


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description="Build LocalDB similarity edges for repeat interaction report.")
    ap.add_argument("--project", required=True)
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--localdb-prefix", required=True, help="BLAST db prefix (projects/<project>/blast_db/multifasta.fasta)")
    ap.add_argument("--out-hsps-tsv", required=True)
    ap.add_argument("--out-edges-full-tsv", required=True)
    ap.add_argument("--out-edges-filtered-tsv", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    project = str(args.project)
    query_fa = Path(args.query_fasta)
    localdb = str(args.localdb_prefix)
    out_hsps = Path(args.out_hsps_tsv)
    out_full = Path(args.out_edges_full_tsv)
    out_filt = Path(args.out_edges_filtered_tsv)
    for p in (out_hsps, out_full, out_filt):
        p.parent.mkdir(parents=True, exist_ok=True)

    cfg = _load_config()
    enabled = bool(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "enabled"]))
    if not enabled:
        out_hsps.write_text("")
        out_full.write_text("")
        out_filt.write_text("")
        return 0

    tasks_keep = list(
        resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "tasks_keep"])
    )
    cov_keep = set(
        resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "coverage_types_keep"])
    )
    min_pident = float(
        resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "min_pident"])
    )
    max_evalue = float(
        resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "max_evalue"])
    )
    max_qlen = int(
        resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "max_query_len_nt"])
    )

    cols = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
    ]
    outfmt = "6 " + " ".join(cols)
    blastn = _blastn_cmd()

    dfs: list[pd.DataFrame] = []
    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        for task in tasks_keep:
            task_out = td_p / f"localdb.{task}.hsps.tsv"
            cmd = (
                blastn
                + [
                    "-query",
                    str(query_fa),
                    "-db",
                    localdb,
                    "-task",
                    str(task),
                    "-num_threads",
                    str(int(args.threads)),
                    "-outfmt",
                    outfmt,
                    "-out",
                    str(task_out),
                ]
            )
            _run(cmd)
            if task_out.exists() and task_out.stat().st_size:
                df = pd.read_csv(task_out, sep="\t", names=cols)
                df["task"] = task
                df["db"] = "local"
                dfs.append(df)

    if not dfs:
        out_hsps.write_text("")
        out_full.write_text("")
        out_filt.write_text("")
        return 0

    hsps = pd.concat(dfs, ignore_index=True)

    # Basic sanitation / keep size reasonable
    hsps["qlen"] = pd.to_numeric(hsps["qlen"], errors="coerce")
    hsps = hsps[hsps["qlen"].fillna(0) <= max_qlen].copy()
    hsps = hsps[hsps["qseqid"].astype(str) != hsps["sseqid"].astype(str)].copy()

    out_hsps.write_text("")  # always materialize for reproducibility
    hsps.to_csv(out_hsps, sep="\t", index=False)

    edges_full = group_hsps_to_pairs(hsps, group_cols=("qseqid", "sseqid", "task", "db"))
    if edges_full.empty:
        out_full.write_text("")
        out_filt.write_text("")
        return 0

    edges_full.to_csv(out_full, sep="\t", index=False)

    edges = edges_full.copy()
    edges = edges[edges["coverage_type"].isin(cov_keep)].copy()
    edges = edges[pd.to_numeric(edges["pident"], errors="coerce").fillna(0) >= min_pident].copy()
    edges = edges[pd.to_numeric(edges["evalue"], errors="coerce").fillna(1.0) <= max_evalue].copy()

    edges.to_csv(out_filt, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

