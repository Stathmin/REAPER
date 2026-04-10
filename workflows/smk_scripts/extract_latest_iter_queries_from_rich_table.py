#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd
from Bio import SeqIO


def _load_seqs_by_id(fa: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for rec in SeqIO.parse(str(fa), "fasta"):
        out[str(rec.id)] = str(rec.seq).replace("\n", "").strip()
    return out


def _read_latest_repeat_ids(rich_top: Path) -> List[str]:
    df = pd.read_csv(rich_top, sep="\t")
    if df.empty:
        return []
    if "iter" in df.columns:
        it = pd.to_numeric(df["iter"], errors="coerce")
        it_max = it.max()
        if pd.notna(it_max):
            df = df[it == it_max].copy()

    col = "repeat" if "repeat" in df.columns else ("qseqid" if "qseqid" in df.columns else "")
    if not col:
        raise KeyError(f"{rich_top} missing expected column 'repeat' (or 'qseqid')")

    vals = [str(x).strip() for x in df[col].tolist()]
    return [v for v in vals if v and v.lower() != "nan"]


def main() -> int:
    ap = argparse.ArgumentParser(description="Extract latest-iter consensus sequences into a query FASTA from a satMiner rich_table.*.top.tsv")
    ap.add_argument("--rich-top-tsv", required=True)
    ap.add_argument("--consensus-fasta", required=True, help="consensus_all_iters.tagged.fasta")
    ap.add_argument("--out-fasta", required=True)
    args = ap.parse_args()

    rich_top = Path(args.rich_top_tsv)
    consensus = Path(args.consensus_fasta)
    out_fa = Path(args.out_fasta)
    out_fa.parent.mkdir(parents=True, exist_ok=True)

    wanted = _read_latest_repeat_ids(rich_top)
    if not wanted:
        out_fa.write_text("")
        return 0

    seqs = _load_seqs_by_id(consensus)
    written = 0
    with out_fa.open("w") as out:
        for rid in wanted:
            seq = seqs.get(rid, "")
            if not seq:
                continue
            out.write(f">{rid}\n{seq}\n")
            written += 1

    if written == 0:
        out_fa.write_text("")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

