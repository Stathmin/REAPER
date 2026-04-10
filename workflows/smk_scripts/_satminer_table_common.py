from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence

import pandas as pd


TAG_RE = re.compile(
    r"^SMPL=(?P<smpl>[^|]+)\|ORG=(?P<org>[^|]+)\|GENOMES=(?P<genomes>[^|]+)\|"
    r"ITER=(?P<iter>[0-9]+)\|RANK=(?P<rank>[0-9]+)\|ORIG=(?P<orig>.+)$"
)
ORIG_CLUSTER_RE = re.compile(r"^CL(?P<cluster>[0-9]+)_")


def parse_tag(tag: str) -> Dict[str, str]:
    m = TAG_RE.match(tag)
    if not m:
        raise ValueError(f"Unexpected tagged consensus id format: {tag}")
    d = m.groupdict()
    return {
        "smpl": d["smpl"],
        "org": d["org"],
        "genomes": d["genomes"],
        "iter": d["iter"],
        "rank": d["rank"],
        "orig": d["orig"],
    }


def infer_cluster_from_orig(orig: str) -> Optional[str]:
    m = ORIG_CLUSTER_RE.match(orig)
    if not m:
        return None
    return m.group("cluster")


def resolve_graph_path(clusters_dir: Path, cluster_num: str) -> str:
    if not clusters_dir.exists():
        return ""

    # Prefer exact padded directory if present.
    for pad in (4, 3, 2, 1):
        d = clusters_dir / f"dir_CL{int(cluster_num):0{pad}d}"
        for name in ("graph_layout.png", "graph_layout.jpg", "graph_layout.jpeg", "graph_layout_tmb.png"):
            p = d / name
            if p.exists():
                return str(p)

    # Fallback: scan a small set of candidates (avoid walking whole tree).
    candidates = sorted(clusters_dir.glob("dir_CL*/graph_layout.png"))
    want = f"dir_CL{int(cluster_num)}"
    for p in candidates:
        if p.parent.name == want:
            return str(p)
        if re.fullmatch(r"dir_CL0*%s" % re.escape(str(int(cluster_num))), p.parent.name):
            return str(p)
    return ""


def _empty_df(columns: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def read_hits(hits_tsv: Path, required: Iterable[str] = ("repeat",)) -> pd.DataFrame:
    if not hits_tsv.exists() or hits_tsv.stat().st_size == 0:
        return _empty_df(["repeat", "hit_count"])
    df = pd.read_csv(hits_tsv, sep="\t")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"hits TSV missing columns {missing}: {hits_tsv}")
    return df


def read_div(div_tsv: Path, required: Iterable[str] = ("repeat",)) -> pd.DataFrame:
    if not div_tsv.exists() or div_tsv.stat().st_size == 0:
        return _empty_df(["repeat", "n", "pdiv_mean", "pdiv_median", "pdiv_q05", "pdiv_q95"])
    df = pd.read_csv(div_tsv, sep="\t")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"divergence TSV missing columns {missing}: {div_tsv}")
    return df


def read_blast_top(top_tsv: Path) -> pd.DataFrame:
    if not top_tsv.exists() or top_tsv.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(top_tsv, sep="\t")


def normalize_blast_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Normalize BLAST columns into stable names for report writers.
    # satminer_blast_x3.py emits best_* columns, but we also tolerate older schemas.
    if "best_sseqid" in df.columns and "best_hit" not in df.columns:
        df["best_hit"] = df.get("best_accession", pd.NA).fillna(df["best_sseqid"])
    if "best_evalue" not in df.columns and "evalue" in df.columns:
        df["best_evalue"] = df["evalue"]
    if "best_pident" not in df.columns and "pident" in df.columns:
        df["best_pident"] = df["pident"]
    if "task_used" not in df.columns and "task" in df.columns:
        df["task_used"] = df["task"]
    return df

