#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO

# Ensure repo root is on sys.path when invoked as a script.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import yaml  # type: ignore[import-not-found]

from workflows.smk_scripts.reportr_config import resolve_value as resolve_cfg_value


def _load_config() -> dict:
    cfg_path = Path("projects/global_config.yaml")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def _discover_rich_tables(project: str) -> List[Path]:
    base = Path("projects") / project
    paths: List[Path] = []
    paths.extend(sorted(base.glob("samples/*/post_tarean/satminer/report/rich_table.tsv")))
    paths.extend(sorted(base.glob("comparative/*/post_tarean/satminer/report/rich_table.tsv")))
    return [p for p in paths if p.exists()]


def _infer_subject(rich_table: Path) -> Tuple[str, str]:
    parts = list(rich_table.parts)
    if "samples" in parts:
        i = parts.index("samples")
        return "sample", parts[i + 1]
    if "comparative" in parts:
        i = parts.index("comparative")
        return "comparative", parts[i + 1]
    return "unknown", rich_table.stem


def _consensus_fasta_for(rich_table: Path) -> Path:
    # .../satminer/report/rich_table.tsv → .../satminer/consensus_all_iters.tagged.fasta
    satminer_dir = rich_table.parent.parent
    cand = satminer_dir / "consensus_all_iters.tagged.fasta"
    if not cand.exists():
        raise FileNotFoundError(str(cand))
    return cand


def _all_repeat_ids(rich_table: Path) -> List[str]:
    df = pd.read_csv(rich_table, sep="\t")
    if df.empty:
        return []
    if "repeat" not in df.columns:
        raise KeyError(f"{rich_table} missing required column 'repeat'")
    vals = [str(x).strip() for x in df["repeat"].tolist()]
    return [v for v in vals if v and v.lower() != "nan"]


def _seqs_by_id(fa: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for rec in SeqIO.parse(str(fa), "fasta"):
        out[str(rec.id)] = str(rec.seq).replace("\n", "").strip()
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build combined query FASTA for repeat interaction report (all repeats present in satMiner rich_table.tsv)."
    )
    ap.add_argument("--project", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-index-tsv", required=True, help="Repeat-to-subject index TSV")
    ap.add_argument(
        "--rich-table",
        nargs="*",
        default=[],
        help="Optional explicit rich_table.tsv paths. When provided, these are used instead of filesystem discovery.",
    )
    args = ap.parse_args()

    project = str(args.project)
    out_fa = Path(args.out_fasta)
    out_idx = Path(args.out_index_tsv)
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    out_idx.parent.mkdir(parents=True, exist_ok=True)

    cfg = _load_config()
    enabled = bool(
        resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "enabled"])
    )
    if not enabled:
        out_fa.write_text("")
        out_idx.write_text("repeat\tsubject_type\tsubject_id\tsource_rich_table\n")
        return 0

    rich_tables = [Path(p) for p in (args.rich_table or [])]
    if rich_tables:
        # Preserve CLI order but make deterministic and ensure files exist.
        rich_tables = sorted(rich_tables, key=lambda p: str(p))
        missing = [str(p) for p in rich_tables if not p.exists()]
        if missing:
            raise FileNotFoundError("Missing rich_table.tsv inputs:\n" + "\n".join(missing))
    else:
        rich_tables = _discover_rich_tables(project)
    if not rich_tables:
        # No "ready" subjects yet. Emit empty outputs so the report rule can
        # short-circuit cleanly (and remain re-runnable when data appears).
        out_fa.write_text("")
        out_idx.write_text("repeat\tsubject_type\tsubject_id\tsource_rich_table\n")
        return 0

    # Deterministic: iterate sources in sorted path order.
    repeat_to_seq: Dict[str, str] = {}
    index_rows: List[dict] = []

    for rich_table in rich_tables:
        subject_type, subject_id = _infer_subject(rich_table)
        consensus = _consensus_fasta_for(rich_table)
        wanted = _all_repeat_ids(rich_table)
        if not wanted:
            continue
        seqs = _seqs_by_id(consensus)
        for rep in wanted:
            seq = seqs.get(rep, "")
            if not seq:
                continue
            # De-duplicate by repeat ID; if the same repeat ID appears twice (shouldn't),
            # keep first encountered for determinism and record all sources in index.
            if rep not in repeat_to_seq:
                repeat_to_seq[rep] = seq
            index_rows.append(
                {
                    "repeat": rep,
                    "subject_type": subject_type,
                    "subject_id": subject_id,
                    "source_rich_table": str(rich_table),
                }
            )

    if not repeat_to_seq:
        out_fa.write_text("")
        pd.DataFrame(index_rows).to_csv(out_idx, sep="\t", index=False)
        return 0

    with out_fa.open("w") as out:
        for rep in sorted(repeat_to_seq.keys()):
            out.write(f">{rep}\n{repeat_to_seq[rep]}\n")

    pd.DataFrame(index_rows).to_csv(out_idx, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

