#!/usr/bin/env python3
"""Summarize Snakemake benchmark TSVs: rank jobs by wall time and I/O (io_out).

Reads files under benchmarks/ (or --benchmark-dir), one data row per file
(Snakemake benchmark format: header + single stats row).

Example:
  python3 scripts/summarize_benchmarks.py --top 25 --sort s
  python3 scripts/summarize_benchmarks.py --sort io_out --by-rule
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any
import yaml


def _load_project_ids(config_path: Path) -> list[str]:
    """Return configured project ids (longest first)."""
    if not config_path.is_file():
        return []
    with config_path.open("r", encoding="utf-8", errors="replace") as f:
        cfg = yaml.safe_load(f) or {}
    projects = (cfg.get("projects", {}) or {})
    if not isinstance(projects, dict):
        return []
    ids = [str(k) for k in projects.keys()]
    ids.sort(key=len, reverse=True)
    return ids


def _parse_stem(stem: str, project_ids: list[str]) -> tuple[str, str]:
    """Return (rule_prefix, tail) by detecting any configured project id in the stem."""
    for pid in project_ids:
        needle = f"_{pid}"
        if needle in stem:
            before, after = stem.split(needle, 1)
            if after.startswith("_"):
                after = after[1:]
            return before, after
    return stem, ""


def _read_benchmark_tsv(path: Path, *, project_ids: list[str]) -> dict[str, Any] | None:
    try:
        with path.open(newline="", encoding="utf-8", errors="replace") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    except OSError:
        return None
    if not rows:
        return None
    row = rows[0]
    out: dict[str, Any] = {"file": str(path), "stem": path.stem}
    for key in ("s", "max_rss", "io_in", "io_out", "cpu_time", "mean_load"):
        raw = row.get(key)
        if raw is None or raw == "":
            out[key] = None
        else:
            try:
                out[key] = float(raw)
            except ValueError:
                out[key] = None
    rule, tail = _parse_stem(path.stem, project_ids)
    out["rule"] = rule
    out["tail"] = tail
    # Sample id if present: KA1, KA2, ...
    sm = re.search(r"\b(KA\d+)\b", path.stem)
    out["sample"] = sm.group(1) if sm else ""
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--benchmark-dir",
        type=Path,
        default=Path("benchmarks"),
        help="Directory containing *.tsv benchmark files (default: benchmarks)",
    )
    ap.add_argument(
        "--sort",
        choices=("s", "io_out", "max_rss", "cpu_time"),
        default="s",
        help="Primary sort key (descending)",
    )
    ap.add_argument("--top", type=int, default=50, help="Max rows to print (default: 50)")
    ap.add_argument(
        "--by-rule",
        action="store_true",
        help="Aggregate by rule prefix: max wall time and sum io_out per rule",
    )
    ap.add_argument("--csv", action="store_true", help="CSV output to stdout")
    args = ap.parse_args()

    project_ids = _load_project_ids(Path("projects/global_config.yaml"))

    root = args.benchmark_dir
    if not root.is_dir():
        print(f"Not a directory: {root}", file=sys.stderr)
        return 1

    rows: list[dict[str, Any]] = []
    for path in sorted(root.glob("*.tsv")):
        rec = _read_benchmark_tsv(path, project_ids=project_ids)
        if rec is None or rec.get("s") is None:
            continue
        rows.append(rec)

    if not rows:
        print("No benchmark TSVs with numeric 's' found.", file=sys.stderr)
        return 1

    key = args.sort

    if args.by_rule:
        agg: dict[str, dict[str, Any]] = {}
        for r in rows:
            rule = str(r.get("rule") or "")
            if rule not in agg:
                agg[rule] = {
                    "rule": rule,
                    "n": 0,
                    "s_max": 0.0,
                    "io_out_sum": 0.0,
                    "max_rss_max": 0.0,
                }
            a = agg[rule]
            a["n"] += 1
            s = float(r["s"])
            a["s_max"] = max(a["s_max"], s)
            if r.get("io_out") is not None:
                a["io_out_sum"] += float(r["io_out"])
            if r.get("max_rss") is not None:
                a["max_rss_max"] = max(a["max_rss_max"], float(r["max_rss"]))
        agg_rows = list(agg.values())
        agg_rows.sort(key=lambda x: x["s_max"], reverse=True)
        agg_rows = agg_rows[: args.top]

        if args.csv:
            w = csv.DictWriter(
                sys.stdout,
                fieldnames=["rule", "n", "s_max", "io_out_sum_mib", "max_rss_max_mib"],
            )
            w.writeheader()
            for a in agg_rows:
                w.writerow(
                    {
                        "rule": a["rule"],
                        "n": a["n"],
                        "s_max": f"{a['s_max']:.2f}",
                        "io_out_sum_mib": f"{a['io_out_sum']:.2f}",
                        "max_rss_max_mib": f"{a['max_rss_max']:.2f}",
                    }
                )
            return 0

        print(
            f"Top {len(agg_rows)} rules by max wall time (s) under {root} "
            f"(n = benchmark files per rule)"
        )
        print(
            f"{'rule':<42} {'n':>5} {'s_max':>10} {'io_out_sum':>14} {'max_rss_max':>12}"
        )
        print("-" * 90)
        for a in agg_rows:
            print(
                f"{a['rule']:<42} {a['n']:>5} {a['s_max']:>10.2f} {a['io_out_sum']:>14.2f} {a['max_rss_max']:>12.2f}"
            )
        print()
        print("Columns: s_max = max seconds per file; io_out_sum = sum of io_out (MiB); max_rss_max = max RSS (MiB).")
        return 0

    rows.sort(key=lambda r: (r.get(key) is not None, r.get(key) or 0.0), reverse=True)
    rows = [r for r in rows if r.get(key) is not None][: args.top]

    if args.csv:
        w = csv.DictWriter(
            sys.stdout,
            fieldnames=["stem", "rule", "sample", "s", "io_out", "max_rss", "cpu_time", "file"],
        )
        w.writeheader()
        for r in rows:
            w.writerow(
                {
                    "stem": r["stem"],
                    "rule": r["rule"],
                    "sample": r.get("sample") or "",
                    "s": f"{r['s']:.2f}",
                    "io_out": f"{r['io_out']:.2f}" if r.get("io_out") is not None else "",
                    "max_rss": f"{r['max_rss']:.2f}" if r.get("max_rss") is not None else "",
                    "cpu_time": f"{r['cpu_time']:.2f}" if r.get("cpu_time") is not None else "",
                    "file": r["file"],
                }
            )
        return 0

    print(f"Top {len(rows)} benchmarks by {key} (desc) under {root}")
    print(
        f"{'s':>10} {'io_out':>12} {'max_rss':>10} {'rule':<36} {'tail':<26}"
    )
    print("-" * 105)
    for r in rows:
        io_out = r.get("io_out")
        rss = r.get("max_rss")
        print(
            f"{r['s']:>10.2f} {io_out if io_out is not None else float('nan'):>12.2f} "
            f"{rss if rss is not None else float('nan'):>10.2f} {r['rule']:<36} {r.get('tail', ''):<26}"
        )
    print()
    print("Units: s = wall seconds; io_out / max_rss = MiB (Snakemake benchmark).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
