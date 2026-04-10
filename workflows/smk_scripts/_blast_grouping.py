from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping

import pandas as pd


def coverage_type(cov: float) -> str:
    if cov >= 0.95:
        return "near-full"
    if cov >= 0.80:
        return "composite"
    if cov >= 0.60:
        return "partial"
    if cov > 0:
        return "weak"
    return "none"


def _norm_interval(a: float, b: float) -> tuple[int, int]:
    x = int(a)
    y = int(b)
    return (x, y) if x <= y else (y, x)


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda t: t[0])
    merged: list[tuple[int, int]] = [intervals[0]]
    for cur in intervals[1:]:
        prev = merged[-1]
        if cur[0] <= prev[1] + 1:
            merged[-1] = (prev[0], max(prev[1], cur[1]))
        else:
            merged.append(cur)
    return merged


def union_len(intervals: list[tuple[int, int]]) -> int:
    return sum((b - a + 1) for a, b in merge_intervals(intervals))


@dataclass(frozen=True)
class BestPickPolicy:
    coverage_order: Mapping[str, int] = None  # type: ignore[assignment]
    task_order: Mapping[str, int] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.coverage_order is None:
            object.__setattr__(
                self,
                "coverage_order",
                {"near-full": 0, "composite": 1, "partial": 2, "weak": 3, "none": 9},
            )
        if self.task_order is None:
            object.__setattr__(
                self,
                "task_order",
                {"megablast": 0, "dc-megablast": 1, "blastn": 2, "blastn-short": 3},
            )


def group_hsps_to_pairs(
    hsps: pd.DataFrame,
    *,
    group_cols: Iterable[str] = ("qseqid", "sseqid", "task", "db"),
) -> pd.DataFrame:
    """
    Convert per-HSP outfmt6 rows into per-(qseqid,sseqid,task,db) grouped rows using
    union-of-HSP coverage on both query and subject.

    Required HSP columns:
      qseqid,sseqid,qstart,qend,sstart,send,qlen,slen,length,evalue,bitscore,pident
    Optional:
      db, task
    """
    if hsps.empty:
        return pd.DataFrame()

    need = {"qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qlen", "slen", "length", "evalue", "bitscore", "pident"}
    missing = [c for c in sorted(need) if c not in hsps.columns]
    if missing:
        raise ValueError(f"HSP table missing required columns: {missing}")

    df = hsps.copy()
    for col in ("qlen", "slen", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "pident"):
        df[col] = pd.to_numeric(df[col], errors="coerce")

    rows: list[dict] = []
    for key, g in df.groupby(list(group_cols), dropna=False):
        g2 = g.dropna(subset=["qstart", "qend", "qlen", "sstart", "send", "slen", "length", "evalue", "bitscore", "pident"])
        if g2.empty:
            continue
        q_intervals = [_norm_interval(a, b) for a, b in zip(g2["qstart"], g2["qend"])]
        s_intervals = [_norm_interval(a, b) for a, b in zip(g2["sstart"], g2["send"])]
        q_union = union_len(q_intervals)
        s_union = union_len(s_intervals)
        qlen = float(g2["qlen"].iloc[0])
        slen = float(g2["slen"].iloc[0])
        qcov_union = (q_union / qlen) if qlen > 0 else float("nan")
        scov_union = (s_union / slen) if slen > 0 else float("nan")
        cov = max(qcov_union, scov_union)
        cov_type = coverage_type(float(cov)) if pd.notna(cov) else "none"

        rep = g2.sort_values(["length", "bitscore"], ascending=[False, False]).iloc[0]
        row = {
            "evalue": float(rep["evalue"]),
            "bitscore": float(rep["bitscore"]),
            "pident": float(rep["pident"]),
            "length": float(rep["length"]),
            "qlen": float(rep["qlen"]),
            "slen": float(rep["slen"]),
            "q_union": int(q_union),
            "s_union": int(s_union),
            "qcov": float(qcov_union),
            "scov": float(scov_union),
            "coverage": float(cov),
            "coverage_type": cov_type,
        }

        if isinstance(key, tuple):
            for kcol, kval in zip(list(group_cols), key):
                row[kcol] = kval
        else:
            row[list(group_cols)[0]] = key

        rows.append(row)
    return pd.DataFrame(rows)


def pick_best_per_query(
    pairs: pd.DataFrame,
    *,
    policy: BestPickPolicy | None = None,
    query_col: str = "qseqid",
    x3_detector_col: str = "db",
) -> pd.DataFrame:
    if pairs.empty:
        return pairs.copy()
    pol = policy or BestPickPolicy()
    df = pairs.copy()
    df["_covp"] = df["coverage_type"].map(pol.coverage_order).fillna(9)
    df["_taskp"] = df["task"].map(pol.task_order).fillna(9) if "task" in df.columns else 9
    if x3_detector_col in df.columns:
        df["_x3p"] = df[x3_detector_col].astype(str).apply(lambda s: 1 if "_x3" in s else 0)
    else:
        df["_x3p"] = 0
    df = df.sort_values(
        [query_col, "_covp", "_taskp", "_x3p", "evalue", "bitscore"],
        ascending=[True, True, True, True, True, False],
    )
    best = df.groupby(query_col, as_index=False).head(1).copy()
    return best.drop(columns=["_covp", "_taskp", "_x3p"], errors="ignore")

