#!/usr/bin/env python3
"""
BLASTn comparison of TAREAN consensus across iterations (e.g. iter2/3 vs iter1).

Typical layout:
  projects/<pid>/samples/<sample>/tarean1/TAREAN_consensus_rank_*.fasta
  projects/<pid>/samples/<sample>/tarean2/...

Also supports a legacy single `tarean/` directory (treated as iteration 1).

Requires `makeblastdb` and `blastn` on PATH (NCBI BLAST+).

Example:
  python3 scripts/compare_iterative_consensus_blast.py \\
    --sample-dir projects/triticeae_F21FTSEUHT1241/samples/KA1

  python3 scripts/compare_iterative_consensus_blast.py \\
    --iter-dirs path/to/tarean1 path/to/tarean2 \\
    --out-tsv redundancy_iter_vs_iter1.tsv
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple

_ITER_DIR = re.compile(r"^tarean(\d+)$")


def _discover_iteration_dirs(sample_dir: Path) -> List[Tuple[int, Path]]:
    out: List[Tuple[int, Path]] = []
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"Not a directory: {sample_dir}")
    for p in sorted(sample_dir.iterdir()):
        if not p.is_dir():
            continue
        m = _ITER_DIR.match(p.name)
        if m:
            out.append((int(m.group(1)), p))
        elif p.name == "tarean":
            out.append((1, p))
    return sorted(out, key=lambda x: x[0])


def _concat_consensus_rank_fastas(tarean_dir: Path, out_fasta: Path) -> int:
    """Concatenate TAREAN_consensus_rank_*.fasta; return number of sequences."""
    parts = sorted(tarean_dir.glob("TAREAN_consensus_rank_*.fasta"))
    if not parts:
        raise FileNotFoundError(f"No TAREAN_consensus_rank_*.fasta under {tarean_dir}")
    nseq = 0
    with out_fasta.open("w") as w:
        for p in parts:
            text = p.read_text()
            w.write(text)
            nseq += sum(1 for line in text.splitlines() if line.startswith(">"))
    return nseq


def _run_blastn(
    subject_fasta: Path,
    query_fasta: Path,
    *,
    threads: int,
) -> List[str]:
    with tempfile.TemporaryDirectory(prefix="iterblast_") as td:
        tdp = Path(td)
        db = tdp / "subj"
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                str(subject_fasta),
                "-dbtype",
                "nucl",
                "-out",
                str(db),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        proc = subprocess.run(
            [
                "blastn",
                "-task",
                "megablast",
                "-query",
                str(query_fasta),
                "-db",
                str(db),
                "-outfmt",
                "6 qseqid sseqid pident length qlen slen",
                "-max_hsps",
                "1",
                "-max_target_seqs",
                "5",
                "-num_threads",
                str(threads),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    return [ln for ln in proc.stdout.strip().splitlines() if ln.strip()]


def _best_hit_per_query(lines: List[str]) -> List[Tuple[str, str, float, int, int, int]]:
    """First row per qseqid (blastn output sorted by query order; take best = first seen)."""
    seen: set[str] = set()
    rows: List[Tuple[str, str, float, int, int, int]] = []
    for ln in lines:
        parts = ln.split("\t")
        if len(parts) < 6:
            continue
        q, s, pident, length, qlen, slen = parts[0], parts[1], float(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
        if q in seen:
            continue
        seen.add(q)
        rows.append((q, s, pident, length, qlen, slen))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description="BLASTn iter2+ consensus vs iter1 reference.")
    ap.add_argument(
        "--sample-dir",
        help="Sample directory containing tarean1/, tarean2/, ... (or legacy tarean/)",
    )
    ap.add_argument(
        "--iter-dirs",
        nargs="*",
        help="Explicit tarean directories in iteration order (overrides --sample-dir discovery)",
    )
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument(
        "--min-pident",
        type=float,
        default=95.0,
        help="Report fraction of queries with best hit pident >= this",
    )
    ap.add_argument("--out-tsv", type=Path, default=None, help="Write per-query best hits")
    args = ap.parse_args()

    if not shutil.which("makeblastdb") or not shutil.which("blastn"):
        print("ERROR: makeblastdb and blastn must be on PATH.", file=sys.stderr)
        return 2

    if args.iter_dirs:
        ordered = [(i + 1, Path(p)) for i, p in enumerate(args.iter_dirs)]
    elif args.sample_dir:
        ordered = _discover_iteration_dirs(Path(args.sample_dir))
    else:
        ap.error("Provide --sample-dir or --iter-dirs")

    if len(ordered) < 2:
        print(
            "Need at least two iteration directories (e.g. tarean1 and tarean2). "
            f"Found: {[p for _, p in ordered]}",
            file=sys.stderr,
        )
        return 2

    ref_iter, ref_dir = ordered[0]
    print(f"Reference iteration: {ref_iter} ({ref_dir})")

    with tempfile.TemporaryDirectory(prefix="itercons_") as tmp:
        tmp = Path(tmp)
        ref_cat = tmp / "ref_concat.fasta"
        nref = _concat_consensus_rank_fastas(ref_dir, ref_cat)
        print(f"  Reference sequences: {nref}")

        first_tsv = True
        for q_iter, q_dir in ordered[1:]:
            q_cat = tmp / f"q_iter{q_iter}.fasta"
            nq = _concat_consensus_rank_fastas(q_dir, q_cat)
            lines = _run_blastn(ref_cat, q_cat, threads=args.threads)
            hits = _best_hit_per_query(lines)
            if not hits:
                print(f"Query iter {q_iter}: no BLAST hits (nq={nq})")
                continue
            high = sum(1 for h in hits if h[2] >= args.min_pident)
            mean_p = sum(h[2] for h in hits) / len(hits)
            frac = high / len(hits) if hits else 0.0
            msg = (
                f"iter{q_iter} vs iter{ref_iter}: queries={len(hits)} "
                f"mean_best_pident={mean_p:.2f} frac_pident>={args.min_pident:g}={frac:.3f}"
            )
            print(msg)
            if args.out_tsv:
                mode = "w" if first_tsv else "a"
                first_tsv = False
                with args.out_tsv.open(mode) as w:
                    w.write(f"# {msg}\n")
                    w.write("qseqid\tsseqid\tpident\tlength\tqlen\tslen\n")
                    for q, s, pi, ln, ql, sl in hits:
                        w.write(f"{q}\t{s}\t{pi:.2f}\t{ln}\t{ql}\t{sl}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
