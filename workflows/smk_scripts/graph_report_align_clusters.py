#!/usr/bin/env python3

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Ensure repo root is on sys.path when invoked as a script.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import yaml  # type: ignore[import-not-found]

from workflows.smk_scripts.reportr_config import resolve_value as resolve_cfg_value


def _load_config() -> dict:
    cfg_path = Path("projects/global_config.yaml")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def _run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def _mafft_cmd() -> list[str]:
    if shutil.which("mafft"):
        return ["mafft"]
    # mafft is expected in repeatexplorer env typically.
    return ["conda", "run", "-n", "repeatexplorer", "mafft"]


def _trim_alignment(aln_path: Path, *, min_non_gap_frac: float = 0.5) -> None:
    aln = AlignIO.read(str(aln_path), "fasta")
    nseq = len(aln)
    if nseq == 0:
        return
    ncol = aln.get_alignment_length()

    def non_gap_frac(col_idx: int) -> float:
        non_gap = 0
        for i in range(nseq):
            c = aln[i, col_idx]
            if c not in ("-", "."):
                non_gap += 1
        return non_gap / nseq

    left = 0
    while left < ncol and non_gap_frac(left) < min_non_gap_frac:
        left += 1
    right = ncol - 1
    while right >= left and non_gap_frac(right) < min_non_gap_frac:
        right -= 1

    if left == 0 and right == ncol - 1:
        return

    trimmed = aln[:, left : right + 1]
    AlignIO.write(trimmed, str(aln_path), "fasta")


def _consensus_from_alignment(aln_path: Path, *, out_fa: Path, consensus_id: str) -> None:
    aln = AlignIO.read(str(aln_path), "fasta")
    if len(aln) == 0:
        out_fa.write_text("")
        return
    ncol = aln.get_alignment_length()
    consensus_chars: list[str] = []
    for j in range(ncol):
        counts: dict[str, int] = {}
        for i in range(len(aln)):
            c = aln[i, j].upper()
            if c in ("-", "."):
                continue
            counts[c] = counts.get(c, 0) + 1
        if not counts:
            continue
        best = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
        consensus_chars.append(best)
    seq = "".join(consensus_chars).replace("U", "T")
    rec = SeqRecord(Seq(seq), id=consensus_id, description="")
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write([rec], str(out_fa), "fasta")


def main() -> int:
    ap = argparse.ArgumentParser(description="Align per-cluster sequences with MAFFT and write consensuses.")
    ap.add_argument("--project", required=True)
    ap.add_argument("--queries-fasta", required=True)
    ap.add_argument("--clusters-tsv", required=True)
    ap.add_argument("--out-dir", required=True, help="Directory to write cluster FASTAs/alignments/consensuses")
    ap.add_argument("--out-consensus-fasta", required=True, help="Combined consensus FASTA (one per cluster)")
    args = ap.parse_args()

    project = str(args.project)
    cfg = _load_config()
    enabled = bool(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "enabled"]))
    if not enabled:
        Path(args.out_consensus_fasta).write_text("")
        return 0

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_cons = Path(args.out_consensus_fasta)
    out_cons.parent.mkdir(parents=True, exist_ok=True)

    mafft_threads = int(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "mafft_threads"]))
    mafft_args = str(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "mafft_args"])).split()
    trim_enabled = bool(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "trim_enabled"]))

    # Load sequences
    seqs = {rec.id: rec for rec in SeqIO.parse(str(Path(args.queries_fasta)), "fasta")}
    cl = Path(args.clusters_tsv)
    df = pd.read_csv(cl, sep="\t")
    if df.empty:
        out_cons.write_text("")
        return 0

    mafft = _mafft_cmd()
    consensus_recs: list[SeqRecord] = []

    for cluster_id, g in df.groupby("cluster_id", dropna=False):
        members = [str(x) for x in g["repeat"].tolist()]
        fa_path = out_dir / f"{cluster_id}.fasta"
        aln_path = out_dir / f"{cluster_id}.aln.fasta"
        cons_path = out_dir / f"{cluster_id}_consensus.fasta"

        recs = [seqs[m] for m in members if m in seqs]
        if not recs:
            continue
        SeqIO.write(recs, str(fa_path), "fasta")

        # MAFFT → alignment
        with aln_path.open("w") as out_f:
            cmd = mafft + mafft_args + ["--thread", str(mafft_threads), str(fa_path)]
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.DEVNULL)

        if trim_enabled:
            _trim_alignment(aln_path, min_non_gap_frac=0.5)

        _consensus_from_alignment(aln_path, out_fa=cons_path, consensus_id=str(cluster_id))
        # Add to combined
        for rec in SeqIO.parse(str(cons_path), "fasta"):
            consensus_recs.append(rec)

    SeqIO.write(consensus_recs, str(out_cons), "fasta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

