#!/usr/bin/env python3
"""Aggregate Snakemake benchmark TSVs into manuscript summary tables.

Reads benchmarks/*.tsv (Snakemake benchmark format) and writes
docs/tables/table3_benchmarks_triticeae.tsv. Run from repository root:

  python3 docs/tables/build_tables.py

If no benchmark files are present, keep the committed summary TSV as the
published artifact and regenerate after a full Triticeae run.
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / "benchmarks"
OUT = Path(__file__).resolve().parent / "table3_benchmarks_triticeae.tsv"

RULES = ("clean_reads", "prepare_reads", "run_tarean")


def _median_iqr(values: list[float]) -> tuple[float, float, float]:
    values = sorted(values)
    med = statistics.median(values)
    n = len(values)
    if n < 2:
        return med, med, med
    q1 = statistics.median(values[: n // 2]) if n % 2 else statistics.median(values[: n // 2])
    q3 = statistics.median(values[(n + 1) // 2 :])
    return med, q1, q3


def main() -> None:
    rows_out: list[tuple[str, float, float, float, float, float, float]] = []
    for rule in RULES:
        files = sorted(BENCH.glob(f"{rule}*.tsv"))
        if not files:
            print(f"No benchmark files matched {rule}; leaving {OUT} unchanged.")
            return
        times: list[float] = []
        rss_bytes: list[float] = []
        for f in files:
            with f.open(newline="") as fp:
                r = csv.DictReader(fp, delimiter="\t")
                for row in r:
                    try:
                        times.append(float(row["s"]))
                        rss_bytes.append(float(row["max_rss"]))
                    except (KeyError, ValueError):
                        continue
        if not times:
            continue
        t_med, t_lo, t_hi = _median_iqr(times)
        r_med, r_lo, r_hi = _median_iqr(rss_bytes)
        gib = 1024.0**3
        rows_out.append(
            (
                rule,
                t_med,
                t_lo,
                t_hi,
                r_med / gib,
                r_lo / gib,
                r_hi / gib,
            )
        )
    if not rows_out:
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fp:
        w = csv.writer(fp, delimiter="\t")
        w.writerow(
            [
                "rule",
                "median_s",
                "IQR_low_s",
                "IQR_high_s",
                "median_peak_rss_gib",
                "IQR_low_rss_gib",
                "IQR_high_rss_gib",
            ]
        )
        for row in rows_out:
            w.writerow([f"{x:.2f}" if isinstance(x, float) else x for x in row])
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
