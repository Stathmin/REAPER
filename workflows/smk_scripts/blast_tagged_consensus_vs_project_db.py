#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class DbChoice:
    prefix: str
    nsq_path: str | None


def _verify_blastdb_prefix(prefix: str) -> None:
    # BLAST db v5 commonly includes these; we require the core trio.
    required = [f"{prefix}.nsq", f"{prefix}.nin", f"{prefix}.nhr"]
    missing = [p for p in required if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(
            "BLAST DB index files missing.\n"
            f"- expected: {', '.join(required)}\n"
            f"- missing: {', '.join(missing)}\n"
            "Hint: (re)build with makeblastdb -dbtype nucl -in <fasta> -out <prefix>"
        )


def _pick_newest_db_prefix(blast_db_dir: Path) -> DbChoice:
    # Prefer actual BLAST DB index files.
    nsq = sorted(blast_db_dir.glob("*.nsq"), key=lambda p: p.stat().st_mtime, reverse=True)
    if nsq:
        nsq_path = nsq[0]
        prefix = str(nsq_path).removesuffix(".nsq")
        _verify_blastdb_prefix(prefix)
        return DbChoice(prefix=prefix, nsq_path=str(nsq_path))

    # Fallback: parse *.njs (BLAST DB v5 metadata) to infer expected files.
    njs = sorted(blast_db_dir.glob("*.njs"), key=lambda p: p.stat().st_mtime, reverse=True)
    if njs:
        meta = json.loads(njs[0].read_text())
        dbname = meta.get("dbname")
        if not dbname:
            raise FileNotFoundError(f"{njs[0]} missing 'dbname' field")
        prefix = str((blast_db_dir / dbname).resolve())
        _verify_blastdb_prefix(prefix)
        return DbChoice(prefix=prefix, nsq_path=None)

    raise FileNotFoundError(
        f"No BLAST database found in {blast_db_dir}.\n"
        "Expected either '*.nsq' (plus matching .nin/.nhr) or a '*.njs' metadata file."
    )


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True, help="Tagged consensus FASTA")
    ap.add_argument("--blast-db-dir", required=True, help="Project blast_db directory")
    ap.add_argument("--out-tsv", required=True, help="Full outfmt6 TSV")
    ap.add_argument("--out-top-tsv", required=True, help="Top-hit per query TSV")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--max-target-seqs", type=int, default=25)
    ap.add_argument(
        "--task",
        default="blastn-short",
        help="BLAST task (default: blastn-short for short consensus sequences)",
    )
    ap.add_argument(
        "--word-size",
        type=int,
        default=7,
        help="BLAST word size (default: 7 for sensitivity on short queries)",
    )
    ap.add_argument(
        "--evalue",
        type=float,
        default=10.0,
        help="BLAST e-value threshold (default: 10.0 for sensitivity)",
    )
    ap.add_argument(
        "--dust",
        default="no",
        choices=["yes", "no"],
        help="Enable DUST filtering (default: no)",
    )
    args = ap.parse_args()

    query = Path(args.query)
    blast_db_dir = Path(args.blast_db_dir)
    out_tsv = Path(args.out_tsv)
    out_top = Path(args.out_top_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_top.parent.mkdir(parents=True, exist_ok=True)

    choice = _pick_newest_db_prefix(blast_db_dir)

    conda_exe = os.environ.get("CONDA_EXE", "conda")
    blastn = ["blastn"]
    # If blastn isn't on PATH in the active env, fall back to reportr env.
    if subprocess.run(["bash", "-lc", "command -v blastn >/dev/null 2>&1"]).returncode != 0:
        blastn = [conda_exe, "run", "-n", "reportr", "blastn"]

    outfmt = (
        "6 qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore qlen slen"
    )
    cmd = (
        blastn
        + [
            "-query",
            str(query),
            "-db",
            choice.prefix,
            "-outfmt",
            outfmt,
            "-task",
            str(args.task),
            "-word_size",
            str(int(args.word_size)),
            "-evalue",
            str(float(args.evalue)),
            "-dust",
            "yes" if str(args.dust).lower() == "yes" else "no",
            "-num_threads",
            str(int(args.threads)),
            "-max_target_seqs",
            str(int(args.max_target_seqs)),
            "-out",
            str(out_tsv),
        ]
    )
    _run(cmd)

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
    if out_tsv.stat().st_size == 0:
        pd.DataFrame(columns=cols).to_csv(out_top, sep="\t", index=False)
        return 0

    df = pd.read_csv(out_tsv, sep="\t", names=cols)
    # Pick best by evalue then bitscore (desc)
    df_sorted = df.sort_values(["qseqid", "evalue", "bitscore"], ascending=[True, True, False])
    top = df_sorted.groupby("qseqid", as_index=False).head(1)
    top.to_csv(out_top, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

