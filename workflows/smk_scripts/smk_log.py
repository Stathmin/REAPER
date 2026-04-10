#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path


def _now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def log_header(*, rule: str, category: str, **kv: str) -> None:
    parts = [f"ts={_now_iso()}", f"RULE={rule}", f"CATEGORY={category}"]
    for k, v in kv.items():
        if v is None:
            continue
        s = str(v)
        if s == "":
            continue
        parts.append(f"{k}={s}")
    print(" ".join(parts), flush=True)


def log_footer(*, ok: bool, elapsed_s: float, **kv: str) -> None:
    parts = [f"ts={_now_iso()}", "STATUS=ok" if ok else "STATUS=fail", f"elapsed_s={elapsed_s:.2f}"]
    for k, v in kv.items():
        if v is None:
            continue
        s = str(v)
        if s == "":
            continue
        parts.append(f"{k}={s}")
    print(" ".join(parts), flush=True)


@dataclass(frozen=True)
class FileSummary:
    path: str
    exists: bool
    size_bytes: int
    n_records: int | None


def summarize_fasta(path: str) -> FileSummary:
    p = Path(path)
    if not p.exists():
        return FileSummary(path=str(p), exists=False, size_bytes=0, n_records=None)
    n = 0
    with p.open("r", errors="replace") as f:
        for ln in f:
            if ln.startswith(">"):
                n += 1
    return FileSummary(path=str(p), exists=True, size_bytes=p.stat().st_size, n_records=n)


def summarize_fastq(path: str) -> FileSummary:
    p = Path(path)
    if not p.exists():
        return FileSummary(path=str(p), exists=False, size_bytes=0, n_records=None)
    # Exact: count reads = lines/4 (fast, streaming).
    lines = 0
    with p.open("r", errors="replace") as f:
        for _ in f:
            lines += 1
    n = lines // 4
    return FileSummary(path=str(p), exists=True, size_bytes=p.stat().st_size, n_records=n)


def main() -> int:
    ap = argparse.ArgumentParser(description="Tiny logging helpers for Snakemake shell rules.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_h = sub.add_parser("header")
    p_h.add_argument("--rule", required=True)
    p_h.add_argument("--category", required=True)
    p_h.add_argument("--kv", action="append", default=[], help="k=v pairs (repeatable)")

    p_f = sub.add_parser("footer")
    p_f.add_argument("--ok", choices=("0", "1"), required=True)
    p_f.add_argument("--elapsed-s", required=True)
    p_f.add_argument("--kv", action="append", default=[], help="k=v pairs (repeatable)")

    p_sf = sub.add_parser("summarize-fasta")
    p_sf.add_argument("--path", required=True)

    p_sq = sub.add_parser("summarize-fastq")
    p_sq.add_argument("--path", required=True)

    args = ap.parse_args()

    def parse_kv(pairs: list[str]) -> dict[str, str]:
        out: dict[str, str] = {}
        for item in pairs:
            if "=" not in item:
                raise ValueError(f"Expected k=v, got: {item!r}")
            k, v = item.split("=", 1)
            out[k] = v
        return out

    if args.cmd == "header":
        log_header(rule=args.rule, category=args.category, **parse_kv(args.kv))
        return 0
    if args.cmd == "footer":
        log_footer(ok=args.ok == "1", elapsed_s=float(args.elapsed_s), **parse_kv(args.kv))
        return 0
    if args.cmd == "summarize-fasta":
        s = summarize_fasta(args.path)
        print(f"path={s.path} exists={int(s.exists)} size_bytes={s.size_bytes} n_records={s.n_records}", flush=True)
        return 0
    if args.cmd == "summarize-fastq":
        s = summarize_fastq(args.path)
        print(f"path={s.path} exists={int(s.exists)} size_bytes={s.size_bytes} n_records={s.n_records}", flush=True)
        return 0
    raise AssertionError(args.cmd)


if __name__ == "__main__":
    raise SystemExit(main())

