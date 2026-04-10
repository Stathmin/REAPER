from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Iterator, Optional, Set

from Bio import SeqIO


def _run(cmd: list[str], *, cwd: Optional[Path] = None) -> None:
    subprocess.run(cmd, check=True, cwd=(str(cwd) if cwd is not None else None))


def require_tools(tools: Iterable[str]) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise RuntimeError(f"Missing required tool(s) in PATH: {', '.join(missing)}")


def parse_rm_out_repeat_hits(rm_out_path: Path) -> Set[str]:
    lines = rm_out_path.read_text(errors="replace").splitlines()
    if any("There were no repetitive sequences detected" in ln for ln in lines[:10]):
        return set()
    hits: set[str] = set()
    for ln in lines[3:]:
        if not ln.strip():
            continue
        parts = ln.replace("(", "").replace(")", "").split()
        if len(parts) < 10:
            continue
        hits.add(parts[9])
    return hits


def rmout_edges(
    *,
    rm_out: Path,
    min_hit_len: int = 0,
) -> Iterator[tuple[str, str, str, str, str, str, int]]:
    lines = rm_out.read_text(errors="replace").splitlines()
    if any("There were no repetitive sequences detected" in ln for ln in lines[:10]):
        return iter(())

    def _iter() -> Iterator[tuple[str, str, str, str, str, str, int]]:
        for ln in lines[3:]:
            if not ln.strip():
                continue
            parts = ln.replace("(", "").replace(")", "").split()
            if len(parts) < 14:
                continue
            score = parts[0]
            pdiv = parts[1]
            query = parts[4]
            qstart = parts[5]
            qend = parts[6]
            hit = parts[9]
            try:
                hitlen = abs(int(qend) - int(qstart)) + 1
            except Exception:
                hitlen = 0
            if hitlen < int(min_hit_len):
                continue
            if query == hit:
                continue
            yield (query, hit, score, pdiv, qstart, qend, hitlen)

    return _iter()


def repeatmasker_pairwise_edges(
    *,
    fasta: Path,
    out_tsv: Path,
    threads: int = 8,
) -> None:
    recs = list(SeqIO.parse(str(fasta), "fasta"))
    seq_by_id = {r.id: str(r.seq) for r in recs}
    ids = list(seq_by_id.keys())

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w") as w:
        w.write("query\thit\n")
        if len(ids) <= 1:
            return

        with tempfile.TemporaryDirectory() as td_s:
            td = Path(td_s)
            qfa = td / "query.fasta"
            libfa = td / "lib.fasta"
            for qid in ids:
                with qfa.open("w") as qh:
                    qh.write(f">{qid}\n{seq_by_id[qid]}\n")
                with libfa.open("w") as lh:
                    for sid in ids:
                        if sid == qid:
                            continue
                        lh.write(f">{sid}\n{seq_by_id[sid]}\n")

                _run(
                    [
                        "RepeatMasker",
                        "-pa",
                        str(int(threads)),
                        "-nolow",
                        "-no_is",
                        "-engine",
                        "rmblast",
                        "-lib",
                        str(libfa),
                        str(qfa),
                    ],
                    cwd=td,
                )
                hits = parse_rm_out_repeat_hits(td / "query.fasta.out")
                for h in sorted(hits):
                    w.write(f"{qid}\t{h}\n")


def blast_all_vs_all_tsv(
    *,
    fasta: Path,
    out_tsv: Path,
    min_pident: float = 80.0,
    min_qcovhsp: float = 50.0,
    threads: int = 8,
) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    if not fasta.exists() or fasta.stat().st_size == 0:
        out_tsv.write_text(
            "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcovhsp\n"
        )
        return

    require_tools(["makeblastdb", "blastn"])

    db_prefix = out_tsv.with_suffix("")
    _run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_prefix)])

    tmp_out = out_tsv.with_suffix(out_tsv.suffix + ".tmp")
    _run(
        [
            "blastn",
            "-query",
            str(fasta),
            "-db",
            str(db_prefix),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovhsp",
            "-num_threads",
            str(int(threads)),
            "-perc_identity",
            str(float(min_pident)),
            "-qcov_hsp_perc",
            str(float(min_qcovhsp)),
            "-max_target_seqs",
            "1000000",
            "-out",
            str(tmp_out),
        ]
    )

    with tmp_out.open("r") as r, out_tsv.open("w") as w:
        w.write(
            "\t".join(
                [
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
                    "qcovhsp",
                ]
            )
            + "\n"
        )
        for line in r:
            if not line.strip():
                continue
            qseqid = line.split("\t", 1)[0]
            sseqid = line.split("\t", 2)[1]
            if qseqid == sseqid:
                continue
            w.write(line)

    tmp_out.unlink(missing_ok=True)

