#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
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


def _sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def _read_fasta_to_dict(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    parts: list[str] = []
    for ln in path.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            if cur_id is not None:
                seqs[cur_id] = "".join(parts).strip()
            cur_id = ln[1:].strip().split()[0]
            parts = []
        else:
            parts.append(ln.strip())
    if cur_id is not None:
        seqs[cur_id] = "".join(parts).strip()
    return seqs


def main() -> int:
    ap = argparse.ArgumentParser(description="Build LocalDB similarity edges with per-source caching.")
    ap.add_argument("--project", required=True)
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--query-index-tsv", required=True)
    ap.add_argument("--localdb-prefix", required=True, help="BLAST db prefix (projects/<project>/blast_db/multifasta.fasta)")
    ap.add_argument("--out-hsps-tsv", required=True)
    ap.add_argument("--out-edges-full-tsv", required=True)
    ap.add_argument("--out-edges-filtered-tsv", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    project = str(args.project)
    query_fa = Path(args.query_fasta)
    index_tsv = Path(args.query_index_tsv)
    localdb_prefix = str(args.localdb_prefix)
    localdb_fa = Path(localdb_prefix)
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

    tasks_keep = list(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "tasks_keep"]))
    cov_keep = set(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "coverage_types_keep"]))
    min_pident = float(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "min_pident"]))
    max_evalue = float(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "max_evalue"]))
    max_qlen = int(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "max_query_len_nt"]))

    if not localdb_fa.exists():
        raise FileNotFoundError(str(localdb_fa))
    if not query_fa.exists():
        raise FileNotFoundError(str(query_fa))
    if not index_tsv.exists():
        raise FileNotFoundError(str(index_tsv))

    # Cache keys
    localdb_fp = _sha256_file(localdb_fa)
    params_fp = _sha256_text(
        json.dumps(
            {
                "tasks_keep": tasks_keep,
                "cov_keep": sorted(cov_keep),
                "min_pident": min_pident,
                "max_evalue": max_evalue,
                "max_qlen": max_qlen,
            },
            sort_keys=True,
        )
    )
    cache_dir = Path(f"projects/{project}/reports/repeat_interaction_graph/cache/local_edges")
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Partition repeats by source rich table.
    idx = pd.read_csv(index_tsv, sep="\t")
    if idx.empty:
        out_hsps.write_text("")
        out_full.write_text("")
        out_filt.write_text("")
        return 0
    if "repeat" not in idx.columns or "source_rich_table" not in idx.columns:
        raise KeyError(f"{index_tsv} missing required columns: repeat, source_rich_table")

    by_source: dict[str, list[str]] = defaultdict(list)
    for _, r in idx.iterrows():
        rep = str(r.get("repeat", "")).strip()
        src = str(r.get("source_rich_table", "")).strip()
        if not rep or not src:
            continue
        by_source[src].append(rep)

    seqs = _read_fasta_to_dict(query_fa)

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
        "task",
        "db",
    ]
    outfmt_cols = cols[:-2]  # task/db not emitted by blast
    outfmt = "6 " + " ".join(outfmt_cols)
    blastn = _blastn_cmd()

    hsps_parts: list[pd.DataFrame] = []

    for src, reps in sorted(by_source.items(), key=lambda kv: kv[0]):
        src_p = Path(src)
        if not src_p.exists():
            raise FileNotFoundError(f"source_rich_table missing: {src_p}")
        src_fp = _sha256_file(src_p)
        src_key = _sha256_text(str(src_p.resolve()))
        key_dir = cache_dir / f"ldb_{localdb_fp[:12]}_p_{params_fp[:12]}"
        key_dir.mkdir(parents=True, exist_ok=True)
        meta_path = key_dir / f"{src_key}.meta.json"
        data_path = key_dir / f"{src_key}.hsps.tsv"

        cache_hit = False
        if meta_path.exists() and data_path.exists() and data_path.stat().st_size:
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                cache_hit = (
                    meta.get("localdb_fp") == localdb_fp
                    and meta.get("params_fp") == params_fp
                    and meta.get("source_fp") == src_fp
                    and meta.get("source_path") == str(src_p)
                )
            except Exception:
                cache_hit = False

        if cache_hit:
            df = pd.read_csv(data_path, sep="\t")
            if not df.empty:
                hsps_parts.append(df)
            continue

        # Cache miss: run BLAST for this source's queries.
        reps = [r for r in reps if r in seqs]
        if not reps:
            meta_path.write_text(
                json.dumps(
                    {"localdb_fp": localdb_fp, "params_fp": params_fp, "source_fp": src_fp, "source_path": str(src_p), "empty": True},
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            data_path.write_text("", encoding="utf-8")
            continue

        with tempfile.TemporaryDirectory() as td:
            td_p = Path(td)
            q_subset = td_p / "query_subset.fasta"
            with q_subset.open("w", encoding="utf-8") as out:
                for rep in reps:
                    seq = seqs.get(rep, "")
                    if not seq:
                        continue
                    out.write(f">{rep}\n{seq}\n")

            dfs: list[pd.DataFrame] = []
            for task in tasks_keep:
                task_out = td_p / f"localdb.{task}.hsps.tsv"
                cmd = (
                    blastn
                    + [
                        "-query",
                        str(q_subset),
                        "-db",
                        localdb_prefix,
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
                    df = pd.read_csv(task_out, sep="\t", names=outfmt_cols)
                    df["task"] = task
                    df["db"] = "local"
                    dfs.append(df)

        if not dfs:
            meta_path.write_text(
                json.dumps(
                    {"localdb_fp": localdb_fp, "params_fp": params_fp, "source_fp": src_fp, "source_path": str(src_p), "empty": True},
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            data_path.write_text("", encoding="utf-8")
            continue

        df_all = pd.concat(dfs, ignore_index=True)
        df_all.to_csv(data_path, sep="\t", index=False)
        meta_path.write_text(
            json.dumps(
                {"localdb_fp": localdb_fp, "params_fp": params_fp, "source_fp": src_fp, "source_path": str(src_p), "rows": int(df_all.shape[0])},
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        hsps_parts.append(df_all)

    if not hsps_parts:
        out_hsps.write_text("")
        out_full.write_text("")
        out_filt.write_text("")
        return 0

    hsps = pd.concat(hsps_parts, ignore_index=True)

    hsps["qlen"] = pd.to_numeric(hsps["qlen"], errors="coerce")
    hsps = hsps[hsps["qlen"].fillna(0) <= max_qlen].copy()
    hsps = hsps[hsps["qseqid"].astype(str) != hsps["sseqid"].astype(str)].copy()

    out_hsps.write_text("")
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

