#!/usr/bin/env python
"""
Aggregate RepOrtR benchmark and validation outputs into tabular summaries
for the manuscript tables and figures described in docs/REPORTR_ARTICLE_DRAFT.md.

Usage examples (run from repo root, in the reportr environment):

  # Aggregate rule-level benchmarks for a given project
  python scripts/aggregate_benchmarks.py benchmarks \
      --project triticeae_F21FTSEUHT1241 \
      --output benchmarks/summary_triticeae.tsv

  # Summarise datasets (Table 1-style)
  python scripts/aggregate_benchmarks.py datasets \
      --output benchmarks/datasets.tsv

  # Summarise validation outcomes (Table 5-style, coarse)
  python scripts/aggregate_benchmarks.py validation \
      --output benchmarks/validation_summary.tsv
"""

import argparse
import csv
import glob
import os
from pathlib import Path
from typing import Dict, List, Tuple

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_benchmark_file(path: Path) -> Dict[str, str]:
    """Parse a Snakemake benchmark TSV (one header row + one data row)."""
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        row = next(reader, None)
        return row or {}


def aggregate_benchmarks(project: str, output: Path) -> None:
    """Aggregate clean_reads / prepare_reads / run_tarean benchmarks for a project."""
    patterns = [
        f"benchmarks/clean_reads_{project}_*.tsv",
        f"benchmarks/prepare_reads_{project}_*.tsv",
        f"benchmarks/run_tarean_{project}_*.tsv",
    ]

    records: List[Dict[str, str]] = []

    for pattern in patterns:
        for path_str in glob.glob(str(REPO_ROOT / pattern)):
            path = Path(path_str)
            basename = path.name  # e.g. clean_reads_triticeae_F21FTSEUHT1241_KA3.tsv
            parts = basename.replace(".tsv", "").split("_")
            rule_name = parts[0]  # clean_reads / prepare_reads / run_tarean
            sample_id = parts[-1]

            data = parse_benchmark_file(path)
            if not data:
                continue

            record = {
                "project": project,
                "sample": sample_id,
                "rule": rule_name,
            }
            record.update(data)
            records.append(record)

    if not records:
        raise SystemExit(f"No benchmark files found for project {project}")

    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(records[0].keys())
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


def summarise_datasets(output: Path) -> None:
    """Summarise projects/samples from projects/global_config.yaml (Table 1-style)."""
    cfg_path = REPO_ROOT / "projects" / "global_config.yaml"
    with cfg_path.open() as fh:
        cfg = yaml.safe_load(fh)

    rows: List[Tuple[str, str, int, float]] = []

    for project_id, project_cfg in cfg.get("projects", {}).items():
        if project_id == "probename_project":
            # Example project in the draft; include for completeness.
            taxonomy = project_cfg.get("taxonomy", "")
        else:
            taxonomy = project_cfg.get("taxonomy", "")

        samples_cfg = project_cfg.get("samples", {})
        sample_count = len(samples_cfg)

        # Genome sizes are optional; use 0.0 if missing.
        total_genome_size = 0.0
        for sample_id, scfg in samples_cfg.items():
            total_genome_size += float(scfg.get("genome_size", 0.0))

        rows.append((project_id, taxonomy, sample_count, total_genome_size))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["project_id", "taxonomy", "sample_count", "sum_genome_size"])
        for row in rows:
            writer.writerow(row)


def summarise_validation(output: Path) -> None:
    """Summarise configuration/structure validation outcomes from logs (Table 5-style, coarse)."""
    logs = {
        "config_validation": REPO_ROOT / "logs" / "config_validation.txt",
        "hardcoded_check": REPO_ROOT / "logs" / "hardcoded_check.txt",
        "tool_validation": REPO_ROOT / "logs" / "tool_validation.txt",
        "legacy_files_check": REPO_ROOT / "logs" / "legacy_files_check.txt",
    }

    rows: List[Tuple[str, str]] = []
    for name, path in logs.items():
        status = "missing"
        if path.exists():
            text = path.read_text()
            if "ERROR" in text or "Error" in text:
                status = "errors_present"
            else:
                status = "ok_or_informational"
        rows.append((name, status))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["check", "status"])
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Aggregate RepOrtR benchmarks and validation logs.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_bench = subparsers.add_parser("benchmarks", help="Aggregate rule-level benchmark TSVs.")
    p_bench.add_argument("--project", required=True, help="Project ID (e.g. triticeae_F21FTSEUHT1241).")
    p_bench.add_argument("--output", required=True, type=Path, help="Output TSV path.")

    p_data = subparsers.add_parser("datasets", help="Summarise projects/samples from global_config.yaml.")
    p_data.add_argument("--output", required=True, type=Path, help="Output TSV path.")

    p_valid = subparsers.add_parser("validation", help="Summarise validation log status.")
    p_valid.add_argument("--output", required=True, type=Path, help="Output TSV path.")

    args = parser.parse_args()

    if args.command == "benchmarks":
        aggregate_benchmarks(args.project, args.output)
    elif args.command == "datasets":
        summarise_datasets(args.output)
    elif args.command == "validation":
        summarise_validation(args.output)
    else:
        parser.error(f"Unknown command {args.command}")


if __name__ == "__main__":
    main()

