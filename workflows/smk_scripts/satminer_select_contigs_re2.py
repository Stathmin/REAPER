#!/usr/bin/env python3
"""
RepOrtR adapter for satMiner step 1c (RepeatExplorer2):
Select contigs per cluster until reaching half of total coverage.

Upstream reference (Python2): third_party/satminer/rexp_get_contigs_re2.py

Inputs:
  - clusters_dir containing dir_CL####/contigs.info.fasta (or any *.fasta with contig info)

Outputs:
  - out_fasta: concatenated contig-info fasta (equivalent to upstream out.txt)
  - out_list: selected contig IDs (one per line), like "CL12Contig7"
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from Bio import SeqIO


def _find_contig_info_fastas(clusters_dir: Path) -> List[Path]:
    # Prefer RE2 contig info files; fall back to any fasta under dir_CL*/.
    candidates = sorted(clusters_dir.glob("dir_CL*/contigs.info.fasta"))
    if candidates:
        return candidates
    return sorted(clusters_dir.glob("dir_CL*/*.fasta"))


_COV_PATTERNS: List[re.Pattern[str]] = [
    # RepeatExplorer2 contigs.info.fasta commonly uses: "CL1Contig1 (169-6.0-1016)"
    # Where the middle field is coverage.
    re.compile(r"\((?P<len>[0-9]+)-(?P<cov>[0-9]+(?:\.[0-9]+)?)-(?P<n>[0-9]+)\)"),
    # Common-ish: "... cov 12.34" or "... cov:12.34"
    re.compile(r"\bcov[:= ]+(?P<cov>[0-9]+(?:\.[0-9]+)?)\b", re.IGNORECASE),
    # Fallback: last floating number in description (allow trailing ')', etc.)
    re.compile(r"(?P<cov>[0-9]+(?:\.[0-9]+)?)(?:\s*[\)\]]\s*)?$"),
]


def _parse_cluster_contig_cov(description: str) -> Tuple[int, str, float]:
    """
    Try to extract:
      - cluster number (int)
      - contig id (string, digits or alnum)
      - coverage (float)

    Upstream splits on: "CL", "Contig", "-", " " and expects:
      CL = info[1], contig = info[2], cov = float(info[4])
    """
    # Accept both "CL1Contig1 ..." and "CL 0001 Contig 1 ..." style headers.
    m = re.search(r"\bCL\s*0*([0-9]+)", description, re.IGNORECASE)
    if not m:
        raise ValueError(f"Could not parse cluster from description: {description!r}")
    cl = int(m.group(1))

    m2 = re.search(r"\bContig\s*([A-Za-z0-9]+)\b", description, re.IGNORECASE)
    if not m2:
        # Some RE2 headers may not include 'Contig' literal; fall back to first token after CL.
        m2 = re.search(r"\bCL\s*0*[0-9]+\s*([A-Za-z0-9]+)\b", description, re.IGNORECASE)
    if not m2:
        raise ValueError(f"Could not parse contig from description: {description!r}")
    contig = m2.group(1)
    # Some header patterns (depending on the fallback match) can yield "Contig76".
    # Normalize so downstream IDs are always like "CL1Contig76".
    if contig.lower().startswith("contig"):
        contig = contig[len("contig") :]

    cov: float | None = None
    for pat in _COV_PATTERNS:
        m3 = pat.search(description)
        if m3:
            try:
                cov_str = m3.groupdict().get("cov") or m3.group(1)
                cov = float(cov_str)
                break
            except ValueError:
                continue
    if cov is None:
        raise ValueError(f"Could not parse coverage from description: {description!r}")

    return cl, contig, cov


def _iter_records(paths: Iterable[Path]):
    for p in paths:
        yield from SeqIO.parse(str(p), "fasta")


def select_contigs_by_coverage_fraction(
    records: Iterable,  # SeqRecord
    coverage_fraction: float,
) -> Dict[int, List[str]]:
    """
    Select contigs per cluster until cumulative coverage reaches
    ``coverage_fraction`` of total cluster coverage (satMiner upstream: 0.5 = half).
    Use 1.0 to include all contigs in the reference (stricter iterative masking).
    """
    if not 0.0 < coverage_fraction <= 1.0:
        raise ValueError(f"coverage_fraction must be in (0, 1], got {coverage_fraction!r}")
    # Build cluster -> contig -> coverage
    covs: Dict[int, Dict[str, float]] = {}
    for rec in records:
        cl, contig, cov = _parse_cluster_contig_cov(rec.description)
        covs.setdefault(cl, {})[contig] = max(covs.setdefault(cl, {}).get(contig, 0.0), cov)

    selected: Dict[int, List[str]] = {}
    for cl, contigs in covs.items():
        sorted_contigs = sorted(contigs.items(), key=lambda kv: kv[1], reverse=True)
        total_cov = sum(c for _, c in sorted_contigs)
        threshold = total_cov * float(coverage_fraction)
        acc = 0.0
        pick: List[str] = []
        for contig, cov in sorted_contigs:
            if acc < threshold:
                pick.append(contig)
                acc += cov
            else:
                break
        selected[cl] = pick

    return dict(sorted(selected.items(), key=lambda kv: kv[0]))


def select_contigs_half_coverage(
    records: Iterable,  # SeqRecord
) -> Dict[int, List[str]]:
    return select_contigs_by_coverage_fraction(records, 0.5)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters-dir", required=True, help="Path to seqclust/clustering/clusters")
    ap.add_argument("--out-fasta", required=True, help="Concatenated contig-info fasta (out.txt equivalent)")
    ap.add_argument("--out-list", required=True, help="Selected contig IDs list (out.list equivalent)")
    ap.add_argument(
        "--coverage-fraction",
        type=float,
        default=0.5,
        help="Target fraction of per-cluster cumulative coverage (0.5=half, 1.0=all contigs)",
    )
    args = ap.parse_args()

    clusters_dir = Path(args.clusters_dir)
    out_fasta = Path(args.out_fasta)
    out_list = Path(args.out_list)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_list.parent.mkdir(parents=True, exist_ok=True)

    fastas = _find_contig_info_fastas(clusters_dir)
    if not fastas:
        raise FileNotFoundError(f"No contig-info FASTA files found under {clusters_dir}")

    # Write concatenated fasta (preserve upstream behavior for downstream extraction).
    with out_fasta.open("w") as w:
        for rec in _iter_records(fastas):
            SeqIO.write(rec, w, "fasta")

    # Re-parse from the concatenated output (mirrors upstream that parses out.txt).
    selected = select_contigs_by_coverage_fraction(
        SeqIO.parse(str(out_fasta), "fasta"),
        args.coverage_fraction,
    )

    with out_list.open("w") as w:
        for cl, contigs in selected.items():
            for contig in contigs:
                w.write(f"CL{cl}Contig{contig}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

