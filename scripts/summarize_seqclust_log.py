#!/usr/bin/env python3
r"""
Summarize RepeatExplorer2 / seqclust driver logs (logs/seqclust_*.log).

seqclust logs do **not** record per-step RSS/CPU; they mix:
  - Wrapper lines (ctime) from RepOrtR `run_tarean_step.py`
  - Python logging with ISO timestamps (when present)
  - Plain stdout from tools (often **without** timestamps)

This script finds **approximate phase durations** by pairing milestone strings with the
most recent ISO timestamp seen before each milestone line. For long stretches of
unstamped BLAST output, the "blast" segment is only bounded when the next stamped
section appears — use **Snakemake benchmark** `benchmarks/run_tarean_*.tsv` for
**whole-rule** wall time, max_rss (MiB), io_out, cpu_time.

Usage:
  python3 scripts/summarize_seqclust_log.py logs/seqclust_triticeae_F21FTSEUHT1241_KA1_iter3.log
  python3 scripts/summarize_seqclust_log.py logs/seqclust_*.log --csv
  python3 scripts/summarize_seqclust_log.py logs/seqclust_...log --benchmark benchmarks/run_tarean_...tsv
"""

from __future__ import annotations

import argparse
import csv
import glob
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


ISO_RE = re.compile(
    r"^(?P<ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}),\d{3}"
)
# e.g. "Building a new DB, current time: 04/02/2026 11:33:22"
LEGACY_TIME_RE = re.compile(r"current time:\s*(?P<m>\d{1,2})/(?P<d>\d{1,2})/(?P<y>\d{4}) (?P<H>\d{2}):(?P<M>\d{2}):(?P<S>\d{2})")
@dataclass(frozen=True)
class Milestone:
    needle: str
    label: str
    occurrence: int = 1  # 1 = first match, 2 = second, ...
    line_match: bool = True


# Order matters for reporting: first match wins per (label, occurrence).
MILESTONES: tuple[Milestone, ...] = (
    Milestone("Starting seqclust process", "seqclust_start"),
    Milestone("Building a new DB", "blast_db_build", 1),
    Milestone("Trying to start Rserve", "rserve_start"),
    Milestone("running in parallel using", "parallel_blast_launch"),
    Milestone("all to all blast finished", "all_to_all_blast_done", 1),
    Milestone("running louvain clustering", "louvain_clustering", 1),
    Milestone("all to all blast finished", "all_to_all_blast_done_2", 2),
    Milestone("running louvain clustering", "louvain_clustering_2", 2),
    Milestone("hitsort with", "hitsort_loaded"),
    Milestone("assembling..", "assembling"),
    Milestone("detecting LTR in assembly", "ltr_detection"),
    Milestone("Creating report for superclusters", "report_superclusters"),
    Milestone("Creating report for individual clusters", "report_clusters"),
    Milestone("Creating main html report", "report_main_html"),
    Milestone("SEQCLUST_EXECUTION_SUMMARY", "summary_appended"),
    Milestone("TAREAN_SUCCESS", "tarean_success"),
)


def _parse_iso_line(line: str) -> datetime | None:
    m = ISO_RE.match(line)
    if not m:
        return None
    return datetime.strptime(m.group("ts"), "%Y-%m-%d %H:%M:%S")


def _parse_legacy_time(line: str) -> datetime | None:
    m = LEGACY_TIME_RE.search(line)
    if not m:
        return None
    return datetime(
        int(m.group("y")),
        int(m.group("m")),
        int(m.group("d")),
        int(m.group("H")),
        int(m.group("M")),
        int(m.group("S")),
    )


def _collect_milestones(path: Path) -> tuple[list[tuple[datetime | None, str, str]], datetime | None, datetime | None]:
    """Returns (events as (time, key, line_snip)), file_start_hint, file_end_hint."""
    text = path.read_text(errors="replace")
    lines = text.splitlines()
    last_iso: datetime | None = None
    # Per-needle occurrence index (lines can repeat the same substring, e.g. two blast rounds).
    needle_line_hits: dict[str, int] = {}

    events: list[tuple[datetime | None, str, str]] = []

    for line in lines:
        iso = _parse_iso_line(line)
        if iso:
            last_iso = iso
        leg = _parse_legacy_time(line)
        if leg and "Building a new DB" in line:
            last_iso = leg

        needles_here = {ms.needle for ms in MILESTONES if ms.needle in line}
        for n in needles_here:
            needle_line_hits[n] = needle_line_hits.get(n, 0) + 1
        for ms in MILESTONES:
            if ms.needle not in line:
                continue
            if needle_line_hits[ms.needle] != ms.occurrence:
                continue
            snip = line.strip()[:140]
            events.append((last_iso, ms.label, snip))
            break

    file_start: datetime | None = None
    file_end: datetime | None = None
    for line in lines[:30]:
        if "Starting TAREAN analysis" in line:
            # "Thu Apr  2 11:33:09 2026: Starting..."
            try:
                part = line.split(":", 1)[0]
                file_start = datetime.strptime(part.strip(), "%a %b %d %H:%M:%S %Y")
            except ValueError:
                pass
    for line in lines[-40:]:
        iso = _parse_iso_line(line)
        if iso:
            file_end = iso
    if file_end is None:
        for line in reversed(lines[-20:]):
            if "Analysis completed for" in line or "TAREAN analysis completed" in line:
                try:
                    part = line.split(":", 1)[0]
                    file_end = datetime.strptime(part.strip(), "%a %b %d %H:%M:%S %Y")
                except ValueError:
                    pass
                break

    return events, file_start, file_end


def _fmt_delta(a: datetime | None, b: datetime | None) -> str:
    if a is None or b is None:
        return "n/a"
    sec = (b - a).total_seconds()
    if sec < 0:
        return "n/a"
    if sec >= 3600:
        return f"{sec / 3600:.2f}h ({sec:.0f}s)"
    if sec >= 60:
        return f"{sec / 60:.1f}m ({sec:.0f}s)"
    return f"{sec:.1f}s"


def _read_benchmark_row(tsv: Path) -> dict[str, str] | None:
    try:
        with tsv.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    except OSError:
        return None
    if not rows:
        return None
    return {k: str(v) for k, v in rows[0].items()}


def _report_one(
    path: Path,
    benchmark: dict[str, str] | None,
    as_csv: bool,
    writer: csv.DictWriter | None,
) -> None:
    events, t0, t1 = _collect_milestones(path)
    rows_out: list[dict[str, str]] = []

    if not as_csv:
        print(f"=== {path} ===")
        if t0 and t1:
            print(f"Wrapper span (ctime first/last lines): {_fmt_delta(t0, t1)}")
        print(
            "Note: phases below use the last ISO timestamp *before* each milestone; "
            "long unstamped BLAST sections are only crudely bounded.\n"
        )

    prev_t: datetime | None = t0
    prev_label = "file_start"
    for t, label, snip in events:
        delta = _fmt_delta(prev_t, t) if prev_t and t else ("n/a" if not t else "0s")
        row = {
            "log": str(path),
            "phase_to": label,
            "delta_from_prev": delta,
            "timestamp": t.isoformat(sep=" ") if t else "",
            "prev": prev_label,
        }
        rows_out.append(row)
        if not as_csv:
            ts = t.strftime("%Y-%m-%d %H:%M:%S") if t else "?"
            print(f"  {prev_label:22} -> {label:22}  {delta:16}  (t≈{ts})")
            if snip and len(snip) < 120:
                print(f"      {snip[:119]}")
        if t:
            prev_t = t
        prev_label = label

    if benchmark:
        if not as_csv:
            print("\nSnakemake benchmark (whole run_tarean job, not per seqclust phase):")
            for k in ("s", "max_rss", "io_out", "cpu_time", "mean_load"):
                if k in benchmark:
                    print(f"  {k}: {benchmark[k]}")
        for row in rows_out:
            row["benchmark_s"] = benchmark.get("s", "")
            row["benchmark_max_rss_mib"] = benchmark.get("max_rss", "")

    if writer and rows_out:
        for row in rows_out:
            writer.writerow(row)
    elif as_csv and not rows_out and writer:
        writer.writerow(
            {
                "log": str(path),
                "phase_to": "",
                "delta_from_prev": "",
                "timestamp": "",
                "prev": "",
                "benchmark_s": benchmark.get("s", "") if benchmark else "",
                "benchmark_max_rss_mib": benchmark.get("max_rss", "") if benchmark else "",
            }
        )

    if not as_csv:
        print()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logs", nargs="+", help="seqclust log path(s) or glob")
    ap.add_argument(
        "--benchmark",
        type=Path,
        help="Optional matching run_tarean benchmark TSV (same stem as log if omitted)",
    )
    ap.add_argument("--csv", action="store_true", help="CSV rows to stdout")
    args = ap.parse_args()

    paths: list[Path] = []
    for pat in args.logs:
        paths.extend(Path(p) for p in sorted(glob.glob(pat)))
    paths = [p for p in paths if p.is_file()]
    if not paths:
        print("No log files matched.", file=sys.stderr)
        return 1

    writer: csv.DictWriter | None = None
    if args.csv:
        fieldnames = [
            "log",
            "phase_to",
            "delta_from_prev",
            "timestamp",
            "prev",
            "benchmark_s",
            "benchmark_max_rss_mib",
        ]
        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
        writer.writeheader()

    for path in paths:
        bench_path = args.benchmark
        if bench_path is None:
            # logs/seqclust_PROJ_KA1_iter3.log -> benchmarks/run_tarean_PROJ_KA1_iter3.tsv
            stem = path.stem
            if stem.startswith("seqclust_"):
                cand = Path("benchmarks") / f"run_tarean_{stem[len('seqclust_'):]}.tsv"
                if cand.is_file():
                    bench_path = cand
        bench: dict[str, str] | None = _read_benchmark_row(bench_path) if bench_path and bench_path.is_file() else None

        _report_one(path, bench, args.csv, writer)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
