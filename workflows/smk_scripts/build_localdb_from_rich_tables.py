#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO


def _infer_origin_and_subject(rich_table: Path) -> Tuple[str, str]:
    parts = rich_table.parts
    if "samples" in parts:
        i = parts.index("samples")
        return "sample", parts[i + 1]
    if "comparative" in parts:
        i = parts.index("comparative")
        return "comparative", parts[i + 1]
    return "unknown", rich_table.stem


def _infer_consensus_fasta(rich_table: Path) -> Optional[Path]:
    # Expected layout:
    #  .../post_tarean/satminer/report/rich_table.tsv
    #  .../post_tarean/satminer/consensus_all_iters.tagged.fasta
    satminer_dir = rich_table.parent.parent  # report/ -> satminer/
    cand = satminer_dir / "consensus_all_iters.tagged.fasta"
    return cand if cand.exists() else None


def _read_repeat_ids(rich_table: Path) -> List[str]:
    df = pd.read_csv(rich_table, sep="\t")
    if df.empty:
        return []
    for col in ("repeat", "qseqid"):
        if col in df.columns:
            vals = [str(x).strip() for x in df[col].tolist()]
            return [v for v in vals if v and v.lower() != "nan"]
    raise KeyError(f"{rich_table} missing expected repeat id column ('repeat' or 'qseqid')")


def _load_sequences_by_id(consensus_fasta: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for rec in SeqIO.parse(str(consensus_fasta), "fasta"):
        out[str(rec.id)] = str(rec.seq).replace("\n", "").strip()
    return out


def _sha256(seq: str) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()


def _wrap_fasta(header: str, seq: str, width: int = 80) -> str:
    seq = seq.replace("\n", "").strip()
    lines = [seq[i : i + width] for i in range(0, len(seq), width)] if seq else [""]
    return f">{header}\n" + "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description="Build canonical LocalDB multifasta.fasta from satMiner rich_table.tsv files")
    ap.add_argument("--project", required=True)
    ap.add_argument(
        "--rich-tsv",
        required=True,
        nargs="+",
        help="One or more rich_table.tsv paths (sample + comparative)",
    )
    ap.add_argument("--out-fasta", required=True, help="Output FASTA path (projects/<project>/blast_db/multifasta.fasta)")
    ap.add_argument("--out-manifest", required=False, default="", help="Optional TSV manifest mapping headers to sources/hashes")
    args = ap.parse_args()

    project = str(args.project)
    rich_paths = [Path(p) for p in args.rich_tsv]
    out_fasta = Path(args.out_fasta)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.out_manifest) if str(args.out_manifest).strip() else None

    # Deterministic ordering: sort by path string.
    rich_paths = sorted(rich_paths, key=lambda p: str(p))

    seen_hashes: Dict[str, str] = {}  # hash -> header
    records: List[str] = []
    manifest_rows: List[Dict[str, str]] = []

    for rich_table in rich_paths:
        if not rich_table.exists():
            raise FileNotFoundError(str(rich_table))

        origin, subject = _infer_origin_and_subject(rich_table)
        consensus = _infer_consensus_fasta(rich_table)
        if consensus is None:
            raise FileNotFoundError(f"Could not locate consensus_all_iters.tagged.fasta adjacent to {rich_table}")

        wanted = _read_repeat_ids(rich_table)
        if not wanted:
            continue

        seqs = _load_sequences_by_id(consensus)
        for rep in wanted:
            seq = seqs.get(rep, "")
            if not seq:
                continue
            h = _sha256(seq)
            header = f"{project}|{origin}|{subject}|{rep}"
            if h in seen_hashes:
                # Keep DB compact: dedupe by exact sequence hash.
                manifest_rows.append(
                    {
                        "sequence_hash": h,
                        "kept_header": seen_hashes[h],
                        "duplicate_header": header,
                        "source_rich_table_tsv": str(rich_table),
                        "source_consensus_fasta": str(consensus),
                    }
                )
                continue
            seen_hashes[h] = header
            records.append(_wrap_fasta(header, seq))
            manifest_rows.append(
                {
                    "sequence_hash": h,
                    "kept_header": header,
                    "duplicate_header": "",
                    "source_rich_table_tsv": str(rich_table),
                    "source_consensus_fasta": str(consensus),
                }
            )

    if not records:
        raise RuntimeError(
            "No sequences collected for LocalDB. Ensure satMiner rich tables exist and include 'repeat' entries."
        )

    # Deterministic output: sort by kept_header (stable across reruns).
    # Note: `records` already in deterministic traversal order, but headers include repeat IDs
    # so we normalize ordering explicitly for reproducibility.
    parsed: List[Tuple[str, str]] = []
    for rec in records:
        header = rec.split("\n", 1)[0][1:]
        parsed.append((header, rec))
    parsed.sort(key=lambda t: t[0])

    with out_fasta.open("w") as f:
        for _, rec in parsed:
            f.write(rec)

    if manifest_path is not None:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(manifest_rows).to_csv(manifest_path, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

