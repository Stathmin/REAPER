#!/usr/bin/env python3

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _blastn_cmd() -> list[str]:
    # Prefer active environment, else fall back to reportr env.
    if subprocess.run(["bash", "-lc", "command -v blastn >/dev/null 2>&1"]).returncode == 0:
        return ["blastn"]
    conda_exe = os.environ.get("CONDA_EXE", "conda")
    return [conda_exe, "run", "-n", "reportr", "blastn"]


def _verify_db_prefix(prefix: str) -> None:
    required = [f"{prefix}.nsq", f"{prefix}.nin", f"{prefix}.nhr"]
    missing = [p for p in required if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            "BLAST DB index files missing.\n"
            f"- expected: {', '.join(required)}\n"
            f"- missing: {', '.join(missing)}\n"
        )


def main() -> int:
    ap = argparse.ArgumentParser(description="BLAST tagged consensus vs an explicit BLAST db prefix (full + best TSV).")
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--db-prefix", required=True, help="Prefix suitable for blastn -db (expects .nsq/.nin/.nhr)")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-top-tsv", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--task", default="megablast")
    ap.add_argument("--evalue", type=float, default=10.0)
    ap.add_argument("--max-target-seqs", type=int, default=25)
    args = ap.parse_args()

    query = Path(args.query_fasta)
    db_prefix = str(args.db_prefix)
    out_tsv = Path(args.out_tsv)
    out_top = Path(args.out_top_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_top.parent.mkdir(parents=True, exist_ok=True)

    _verify_db_prefix(db_prefix)

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
    cmd = (
        blastn
        + [
            "-query",
            str(query),
            "-db",
            db_prefix,
            "-outfmt",
            outfmt,
            "-task",
            str(args.task),
            "-evalue",
            str(float(args.evalue)),
            "-num_threads",
            str(int(args.threads)),
            "-max_target_seqs",
            str(int(args.max_target_seqs)),
            "-out",
            str(out_tsv),
        ]
    )
    _run(cmd)

    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        pd.DataFrame(columns=cols).to_csv(out_top, sep="\t", index=False)
        return 0

    df = pd.read_csv(out_tsv, sep="\t", names=cols)
    df_sorted = df.sort_values(["qseqid", "evalue", "bitscore"], ascending=[True, True, False])
    top = df_sorted.groupby("qseqid", as_index=False).head(1)
    top.to_csv(out_top, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

