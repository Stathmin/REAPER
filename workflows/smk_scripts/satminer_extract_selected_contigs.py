#!/usr/bin/env python3
"""
RepOrtR adapter for satMiner step 1c extraction:
Given a FASTA (contig-info concatenation) and a list of contig IDs (e.g. CL12Contig7),
extract matching sequences into a new FASTA.

Upstream reference: third_party/satminer/extract_seq.py
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterable, Set

from Bio import SeqIO


def _load_ids(path: Path) -> Set[str]:
    ids: Set[str] = set()
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        ids.add(s)
    return ids


def _matches_any_id(description: str, wanted: Set[str]) -> bool:
    # Most robust: check if the synthetic id appears in the description.
    for wid in wanted:
        if wid in description:
            return True

    # Fallback: parse CL/Contig and reconstruct, then compare.
    m = re.search(r"\bCL\s*0*([0-9]+)\b", description, re.IGNORECASE)
    m2 = re.search(r"\bContig\s*([A-Za-z0-9]+)\b", description, re.IGNORECASE)
    if m and m2:
        synth = f"CL{int(m.group(1))}Contig{m2.group(1)}"
        return synth in wanted
    return False


def iter_selected_records(fasta_path: Path, wanted: Set[str]):
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        if _matches_any_id(rec.description, wanted):
            yield rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="Input FASTA containing contigs")
    ap.add_argument("--list", required=True, help="List of selected contig IDs (out.list)")
    ap.add_argument("--out", required=True, help="Output FASTA with selected contigs")
    args = ap.parse_args()

    fasta_path = Path(args.fasta)
    list_path = Path(args.list)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    wanted = _load_ids(list_path)
    if not wanted:
        raise ValueError(f"Empty contig list: {list_path}")

    selected = list(iter_selected_records(fasta_path, wanted))
    if not selected:
        raise RuntimeError(
            f"No contigs matched list entries (n={len(wanted)}) in {fasta_path}. "
            "This usually indicates an unexpected header format."
        )

    with out_path.open("w") as w:
        SeqIO.write(selected, w, "fasta")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

