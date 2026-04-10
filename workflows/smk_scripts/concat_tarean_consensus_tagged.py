#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
from pathlib import Path

from Bio import SeqIO


_RANK_RE = re.compile(r"TAREAN_consensus_rank_(?P<rank>[0-9]+)\.fasta$")
_ITER_RE = re.compile(r"/tarean(?P<iter>[0-9]+)/")


def _norm_tag_value(v: str | None) -> str:
    s = ("" if v is None else str(v)).strip()
    if not s:
        return "None"
    # Keep headers parseable and avoid breaking our pipe-separated tags.
    s = s.replace("|", "_")
    s = re.sub(r"\s+", "_", s)
    return s


def _infer_rank(path: Path) -> int:
    m = _RANK_RE.search(path.name)
    if not m:
        raise ValueError(f"Could not infer rank from filename: {path.name}")
    return int(m.group("rank"))


def _infer_iter(path: Path) -> int:
    m = _ITER_RE.search(path.as_posix())
    if not m:
        raise ValueError(f"Could not infer iter from path: {path}")
    return int(m.group("iter"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--organism", default="None")
    ap.add_argument("--genomes", default="None")
    ap.add_argument("--inputs", nargs="+", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    sample = str(args.sample)
    organism = _norm_tag_value(args.organism)
    genomes = _norm_tag_value(args.genomes)
    inputs = [Path(p) for p in args.inputs]

    records = []
    for p in inputs:
        if not p.exists() or p.stat().st_size == 0:
            continue
        iter_n = _infer_iter(p)
        rank = _infer_rank(p)
        for r in SeqIO.parse(str(p), "fasta"):
            orig_id = str(r.id)
            new_id = (
                f"SMPL={sample}|ORG={organism}|GENOMES={genomes}|ITER={iter_n}|RANK={rank}|ORIG={orig_id}"
            )
            r.id = new_id
            r.name = new_id
            r.description = ""
            records.append(r)

    SeqIO.write(records, str(out), "fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

