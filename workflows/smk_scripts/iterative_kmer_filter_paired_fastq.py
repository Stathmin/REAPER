"""
Paired-end k-mer filtering to unmapped reads (BBTools).

Mirrors the shell pipeline in workflows/iterative_tarean_rules.smk rule
`deconseq_filter_reads_iter` (repair -> interleaved reformat -> bbduk -> split).
Used by workflows/smk_scripts/prepare_comparative_reads.py for iterative comparatives so
filter semantics stay aligned with single-sample iterative runs.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path


def _run(cmd: list[str], *, cwd: Path | None = None) -> None:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed (exit {proc.returncode}): {' '.join(cmd)}\n{proc.stderr}"
        )


def filter_paired_reads_to_unmapped(
    r1: Path,
    r2: Path,
    ref_fasta: Path,
    out_r1: Path,
    out_r2: Path,
    job_tmp: Path,
    *,
    kmer_k: int,
    threads: int,
    minkmerhits: int = 1,
    removeifeitherbad: bool = True,
) -> None:
    """
    K-mer filter paired reads to an unmapped interleaved stream (bbduk), then split to PE.

    ``removeifeitherbad`` matches BBDuk: if True, discard the pair when either mate
    matches the reference (stricter); if False, discard only when both mates match.
    """
    job_tmp.mkdir(parents=True, exist_ok=True)
    paired_r1 = job_tmp / "paired_R1.fq"
    paired_r2 = job_tmp / "paired_R2.fq"
    singletons = job_tmp / "singletons.fq"
    paired_int = job_tmp / "paired_int.fq"
    unmapped_int = job_tmp / "unmapped_int.fq"

    _run(
        [
            "repair.sh",
            f"in={r1}",
            f"in2={r2}",
            f"out1={paired_r1}",
            f"out2={paired_r2}",
            f"outs={singletons}",
            "repair=t",
            "trimreaddescription=t",
            "ain=t",
            "overwrite=t",
        ]
    )
    _run(
        [
            "reformat.sh",
            f"in1={paired_r1}",
            f"in2={paired_r2}",
            f"out={paired_int}",
            "ow=t",
        ]
    )
    rib = "t" if removeifeitherbad else "f"
    _run(
        [
            "bbduk.sh",
            f"in={paired_int}",
            "interleaved=t",
            f"ref={ref_fasta}",
            f"k={int(kmer_k)}",
            f"threads={int(threads)}",
            "rcomp=t",
            f"minkmerhits={int(minkmerhits)}",
            f"removeifeitherbad={rib}",
            f"out={unmapped_int}",
            "overwrite=t",
        ]
    )
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    _run(
        [
            "reformat.sh",
            f"in={unmapped_int}",
            "interleaved=t",
            f"out1={out_r1}",
            f"out2={out_r2}",
            "ow=t",
        ]
    )
