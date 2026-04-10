#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import shutil
import tempfile
from pathlib import Path
from typing import List

import pandas as pd

def _ensure_repo_on_path() -> None:
    # Allow importing `workflows.*` when invoked by snakemake / subprocess.
    import sys

    repo_root = Path(__file__).resolve().parents[2]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))


_ensure_repo_on_path()

from workflows.smk_scripts._blast_grouping import group_hsps_to_pairs, pick_best_per_query


def _write_x3_fasta(in_fa: Path, out_fa: Path) -> None:
    # Simple FASTA tripling; preserve headers exactly (qseqid == repeat tag).
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with in_fa.open() as f, out_fa.open("w") as out:
        header = None
        seq_parts: List[str] = []
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).strip()
                    out.write(header + "\n")
                    out.write((seq * 3) + "\n")
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if header is not None:
            seq = "".join(seq_parts).strip()
            out.write(header + "\n")
            out.write((seq * 3) + "\n")


def _run(cmd: List[str]) -> None:
    import subprocess

    subprocess.run(cmd, check=True)


def _blastn_cmd() -> List[str]:
    if shutil.which("blastn"):
        return ["blastn"]
    # Fallback: reach BLAST+ from the reportr env.
    return ["conda", "run", "-n", "reportr", "blastn"]

def _blastdbcmd_cmd() -> List[str]:
    if shutil.which("blastdbcmd"):
        return ["blastdbcmd"]
    return ["conda", "run", "-n", "reportr", "blastdbcmd"]


_SPECIES_RE = re.compile(r"\[(?P<species>[^\]]+)\]\s*$")


def _parse_title(title: str) -> dict[str, str]:
    t = (title or "").strip()
    m = _SPECIES_RE.search(t)
    species = m.group("species") if m else ""
    name = t[: m.start()].strip() if m else t
    return {"best_title": t, "best_species": species, "best_name": name}

def _oligo_fitting(*, length: float, pident: float, gapopen: float, scov: float, evalue: float) -> str:
    # Primer/oligo hits are short, so query coverage is not informative. Use subject coverage + strict quality.
    if length >= 16 and pident >= 80 and gapopen == 0 and scov >= 0.80:
        if evalue <= 1e-3:
            return "oligo"
        if evalue <= 50:
            return "oligo_weak"
    return "weak"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", required=True)
    ap.add_argument("--query-fasta", required=True)
    ap.add_argument("--blast-db-prefix", required=True, help="Prefix suitable for blastn -db")
    ap.add_argument("--oligo-db-prefix", default="", help="Optional BLAST db prefix for oligo/primer db")
    ap.add_argument("--out-full-tsv", required=True)
    ap.add_argument("--out-best-tsv", required=True)
    ap.add_argument("--raw-hsps-dir", default="", help="If set, persist raw per-HSP outfmt6 TSVs here")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    query = Path(args.query_fasta)
    db_prefix = str(args.blast_db_prefix)
    oligo_db_prefix = str(args.oligo_db_prefix or "").strip()
    out_full = Path(args.out_full_tsv)
    out_best = Path(args.out_best_tsv)
    out_full.parent.mkdir(parents=True, exist_ok=True)
    out_best.parent.mkdir(parents=True, exist_ok=True)

    blastn = _blastn_cmd()
    tasks = ["megablast", "dc-megablast", "blastn"]

    raw_dir = Path(str(args.raw_hsps_dir)).resolve() if str(args.raw_hsps_dir).strip() else None
    if raw_dir is not None:
        raw_dir.mkdir(parents=True, exist_ok=True)

    # Build an x3 query fasta and run task ladder for x1 + x3 queries.
    with tempfile.TemporaryDirectory() as td:
        td_p = Path(td)
        q_x3 = td_p / "queries_x3.fasta"
        _write_x3_fasta(query, q_x3)

        cols = [
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
        ]
        outfmt = "6 " + " ".join(cols)

        dfs: List[pd.DataFrame] = []
        for db_name, qfa in [("ncbi", query.resolve()), ("ncbi_x3", q_x3.resolve())]:
            for task in tasks:
                out_tsv = (raw_dir / f"{db_name}.{task}.hsps.tsv") if raw_dir is not None else (td_p / f"{db_name}.{task}.tsv")
                cmd = (
                    blastn
                    + [
                        "-query",
                        str(qfa),
                        "-db",
                        db_prefix,
                        "-task",
                        task,
                        *([] if task != "blastn" else ["-word_size", "6"]),
                        "-num_threads",
                        str(int(args.threads)),
                        "-outfmt",
                        outfmt,
                        "-out",
                        str(out_tsv),
                    ]
                )
                _run(cmd)
                # If we are persisting raw HSPs, always materialize the file
                # (even if empty) so users can debug "no hits" cases.
                if raw_dir is not None and not out_tsv.exists():
                    out_tsv.write_text("")

                if out_tsv.exists() and out_tsv.stat().st_size:
                    df = pd.read_csv(out_tsv, sep="\t", names=cols)
                    df["task"] = task
                    df["db"] = db_name
                    dfs.append(df)

        # Optional: oligo/primer db (blastn-short, no masking).
        if oligo_db_prefix:
            task = "blastn-short"
            for db_name, qfa in [("oligo", query.resolve())]:
                out_tsv = (raw_dir / f"{db_name}.{task}.hsps.tsv") if raw_dir is not None else (td_p / f"{db_name}.{task}.tsv")
                cmd = (
                    blastn
                    + [
                        "-query",
                        str(qfa),
                        "-db",
                        oligo_db_prefix,
                        "-task",
                        task,
                        "-dust",
                        "no",
                        "-soft_masking",
                        "false",
                        "-num_threads",
                        str(int(args.threads)),
                        "-outfmt",
                        outfmt,
                        "-out",
                        str(out_tsv),
                    ]
                )
                _run(cmd)
                if raw_dir is not None and not out_tsv.exists():
                    out_tsv.write_text("")
                if out_tsv.exists() and out_tsv.stat().st_size:
                    df = pd.read_csv(out_tsv, sep="\t", names=cols)
                    df["task"] = task
                    df["db"] = db_name
                    dfs.append(df)

        if not dfs:
            out_full.write_text("")
            # header-only best table
            pd.DataFrame(
                columns=[
                    "qseqid",
                    "sseqid",
                    "task",
                    "db",
                    "evalue",
                    "bitscore",
                    "pident",
                    "length",
                    "qlen",
                    "slen",
                    "qcov",
                    "scov",
                    "coverage",
                    "coverage_type",
                    "task_used",
                    "oligo_best_sseqid",
                    "oligo_best_pident",
                    "oligo_best_length",
                    "oligo_best_evalue",
                    "oligo_best_bitscore",
                    "oligo_best_scov",
                    "oligo_best_sum_scov",
                    "oligo_fitting",
                    "oligo_hits_topN",
                ]
            ).to_csv(out_best, sep="\t", index=False)
            return 0

        full = pd.concat(dfs, ignore_index=True)

        full["qlen"] = pd.to_numeric(full["qlen"], errors="coerce")
        full["slen"] = pd.to_numeric(full["slen"], errors="coerce")
        full["length"] = pd.to_numeric(full["length"], errors="coerce")
        full["qstart"] = pd.to_numeric(full["qstart"], errors="coerce")
        full["qend"] = pd.to_numeric(full["qend"], errors="coerce")
        full["sstart"] = pd.to_numeric(full["sstart"], errors="coerce")
        full["send"] = pd.to_numeric(full["send"], errors="coerce")
        full["evalue"] = pd.to_numeric(full["evalue"], errors="coerce")
        full["bitscore"] = pd.to_numeric(full["bitscore"], errors="coerce")

        # ------------------------------------------------------------------
        # Legacy-style union-of-HSP coverage per (qseqid,sseqid,task,db)
        # ------------------------------------------------------------------
        per_group = group_hsps_to_pairs(full, group_cols=("qseqid", "sseqid", "task", "db"))
        if not per_group.empty:
            per_group["task_used"] = per_group["task"]
        if per_group.empty:
            out_full.write_text("")
            per_group.to_csv(out_best, sep="\t", index=False)
            return 0

        # Best hit per qseqid after legacy sorting.
        best = pick_best_per_query(per_group, query_col="qseqid", x3_detector_col="db")

        # Enrich best hits with accession/species/name/title via blastdbcmd.
        blastdbcmd = _blastdbcmd_cmd()
        # blastdbcmd expects DB entry IDs; our sseqid values may include extra pipe-suffixes.
        sids = [str(x) for x in best["sseqid"].dropna().unique()]
        sid_to_entry = {sid: sid.split("|")[0].strip() for sid in sids}
        entries = sorted(set(sid_to_entry.values()))
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as f_ids:
            for ent in entries:
                f_ids.write(ent + "\n")
            ids_path = f_ids.name
        titles_tsv = td_p / "blastdbcmd_titles.tsv"
        try:
            # Some BLAST DBs lack accession fields; %i (entry id) + %t (title) is robust.
            import subprocess

            proc = subprocess.run(
                blastdbcmd
                + [
                    "-db",
                    db_prefix,
                    "-entry_batch",
                    ids_path,
                    "-outfmt",
                    "%i\t%t",
                    "-out",
                    str(titles_tsv),
                ],
                capture_output=True,
                text=True,
            )
            # blastdbcmd returns non-zero if *any* entries are missing ("Skipped ..."),
            # but it can still write a useful partial titles TSV. Treat missing entries
            # as non-fatal enrichment loss.
            if proc.returncode != 0 and not (titles_tsv.exists() and titles_tsv.stat().st_size):
                titles_tsv.write_text("")
        except Exception:
            # If blastdbcmd isn't usable, keep enrichment empty but do not fail BLAST step.
            titles_tsv.write_text("")

        entry_to_meta: dict[str, dict[str, str]] = {}
        if titles_tsv.exists() and titles_tsv.stat().st_size:
            for ln in titles_tsv.read_text().splitlines():
                parts = ln.split("\t", 1)
                if len(parts) != 2:
                    continue
                ent, title = parts
                # Use the DB entry id as accession-like identifier.
                meta = {"best_accession": ent.strip()}
                meta.update(_parse_title(title))
                entry_to_meta[ent.strip()] = meta

        def _enrich_row(sid: str) -> dict[str, str]:
            s = str(sid)
            ent = sid_to_entry.get(s, s.split("|")[0].strip())
            if ent in entry_to_meta:
                return entry_to_meta[ent]
            # Fallback: derive accession from sseqid.
            species = ""
            if "|" in s:
                # Common pattern in our DB: ACC|SpeciesOrGenus
                species = s.split("|", 1)[1].strip()
            return {"best_accession": ent, "best_title": "", "best_species": species, "best_name": ""}

        enrich = best["sseqid"].map(_enrich_row)
        best = pd.concat([best, enrich.apply(pd.Series)], axis=1)
        best = best.rename(columns={"sseqid": "best_sseqid"})
        best["best_hit"] = best["best_accession"].fillna(best["best_sseqid"])
        best = best.rename(
            columns={
                "evalue": "best_evalue",
                "pident": "best_pident",
                "bitscore": "best_bitscore",
            }
        )

        # ------------------------------------------------------------
        # Oligo/primer evidence (multi-hit + fitting), joined into best
        # ------------------------------------------------------------
        if oligo_db_prefix and (per_group["db"] == "oligo").any():
            og = per_group[per_group["db"] == "oligo"].copy()
            # Use subject coverage and strict quality for oligo fitting.
            og["oligo_fitting"] = og.apply(
                lambda r: _oligo_fitting(
                    length=float(r.get("length", 0.0)),
                    pident=float(r.get("pident", 0.0)),
                    gapopen=float(r.get("gapopen", 0.0)),
                    scov=float(r.get("scov", 0.0)),
                    evalue=float(r.get("evalue", 1.0)),
                ),
                axis=1,
            )
            fit_rank = {"oligo": 0, "oligo_weak": 1, "weak": 9}
            og["_fitp"] = og["oligo_fitting"].map(fit_rank).fillna(9)
            og = og.sort_values(
                ["qseqid", "_fitp", "evalue", "bitscore", "length"],
                ascending=[True, True, True, False, False],
            )
            # Best oligo hit per qseqid.
            obest = og.groupby("qseqid", as_index=False).head(1).copy()
            obest = obest.rename(
                columns={
                    "sseqid": "oligo_best_sseqid",
                    "pident": "oligo_best_pident",
                    "length": "oligo_best_length",
                    "evalue": "oligo_best_evalue",
                    "bitscore": "oligo_best_bitscore",
                    "scov": "oligo_best_scov",
                }
            )
            # Multi-hit string (top N distinct oligos per qseqid).
            top_n = 10
            top_rows = og.groupby("qseqid", as_index=False).head(top_n).copy()

            def _fmt_hit(r) -> str:
                return (
                    f"{r['sseqid']}|fit={r['oligo_fitting']}|pident={float(r['pident']):.1f}|len={int(float(r['length']))}"
                    f"|scov={float(r.get('scov', 0.0)):.3f}|e={float(r.get('evalue', 1.0)):.2g}"
                )

            hits_top = (
                top_rows.groupby("qseqid")  # type: ignore[call-arg]
                .apply(lambda g: ";".join(_fmt_hit(r) for _, r in g.iterrows()))
                .rename("oligo_hits_topN")
                .reset_index()
            )
            # Join into best
            keep_cols = [
                "qseqid",
                "oligo_best_sseqid",
                "oligo_best_pident",
                "oligo_best_length",
                "oligo_best_evalue",
                "oligo_best_bitscore",
                "oligo_best_scov",
                "oligo_fitting",
            ]
            best = best.merge(obest[keep_cols], on="qseqid", how="left")
            best = best.merge(hits_top, on="qseqid", how="left")

        # Write outputs
        per_group.drop(columns=["_covp", "_taskp", "_x3p"], errors="ignore").to_csv(out_full, sep="\t", index=False)
        best.to_csv(out_best, sep="\t", index=False)
        return 0


if __name__ == "__main__":
    raise SystemExit(main())

