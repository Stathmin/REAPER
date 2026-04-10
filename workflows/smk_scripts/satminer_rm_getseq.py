#!/usr/bin/env python3
"""
Python3 adapter for satMiner rm_getseq.py:
Extract RepeatMasker hits from a FASTA according to a RepeatMasker .out file.
Writes <rm_out>.fas like upstream.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--rm-out", required=True, help="RepeatMasker .out file")
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--out", required=True, help="Output FASTA of extracted sequences")
    args = ap.parse_args()

    fafile = Path(args.fasta)
    rmfile = Path(args.rm_out)
    outpath = Path(args.out)
    outpath.parent.mkdir(parents=True, exist_ok=True)

    seqs = {str(r.id): str(r.seq) for r in SeqIO.parse(str(fafile), "fasta")}
    rmout = rmfile.read_text().splitlines()

    with outpath.open("w") as out:
        if any("There were no repetitive sequences detected" in ln for ln in rmout[:10]):
            return 0
        for line in rmout[3:]:
            if not line.strip():
                continue
            line = line.replace("(", "").replace(")", "")
            info = line.split()
            if len(info) < 14:
                continue
            name = info[4]
            begin_q = int(info[5])
            end_q = int(info[6])
            sense = info[8]
            begin_r = int(info[11])
            end_r = int(info[12])
            left_r = int(info[13])
            double = info[15] if len(info) > 15 else ""

            if double:
                continue

            if name not in seqs:
                continue

            if sense == "+":
                len_rep = end_r - begin_r + 1
                if len_rep < args.min_len:
                    continue
                secu = seqs[name][begin_q - 1 : end_q]
                out.write(f">{name}\n{secu}\n")
            elif sense == "C":
                len_rep = end_r - left_r + 1
                if len_rep < args.min_len:
                    continue
                secu = seqs[name][begin_q - 1 : end_q]
                secu_inv = Seq(secu).reverse_complement()
                out.write(f">{name}\n{secu_inv}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

