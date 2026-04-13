#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

_RE_REF_PIPE = re.compile(r"ref\|([^|]+)\|", re.IGNORECASE)
_RE_ACCESSION = re.compile(r"^[A-Za-z0-9._-]+$")

# Ensure repo root is on sys.path when invoked as a script.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import yaml  # type: ignore[import-not-found]

from workflows.smk_scripts.reportr_config import resolve_value as resolve_cfg_value


def _load_config() -> dict:
    cfg_path = Path("projects/global_config.yaml")
    return yaml.safe_load(cfg_path.read_text(encoding="utf-8")) or {}


def _pick_newest_ncbi_metadata_csv(project: str) -> Path | None:
    meta = Path("projects") / project / "metadata"
    if not meta.is_dir():
        return None
    cands = sorted(meta.glob("ncbi_metadata_*.csv"), key=lambda p: p.stat().st_mtime, reverse=True)
    return cands[0] if cands else None


def _load_ncbi_metadata_map(path: Path) -> dict[str, tuple[str, str]]:
    df = pd.read_csv(path)
    need = {"accession", "organism", "description"}
    missing = need - set(df.columns)
    if missing:
        return {}
    out: dict[str, tuple[str, str]] = {}
    for _, row in df.iterrows():
        acc = str(row["accession"]).strip() if pd.notna(row["accession"]) else ""
        if not acc:
            continue
        org = str(row["organism"]).strip() if pd.notna(row["organism"]) else ""
        desc = str(row["description"]).strip() if pd.notna(row["description"]) else ""
        desc = desc.replace("\n", " ").replace("\r", " ")
        out[acc] = (org, desc)
    return out


def _parse_versioned_accession(row: pd.Series) -> str:
    """Return NCBI accession.version for metadata lookup (exact match only)."""
    ba = str(row.get("best_accession", "") or "").strip().rstrip("|")
    if ba and ba.lower() != "ref" and _RE_ACCESSION.fullmatch(ba) and "." in ba:
        return ba
    ss = str(row.get("best_sseqid", "") or "").strip()
    m = _RE_REF_PIPE.search(ss)
    if m:
        return m.group(1).strip()
    if ss and not ss.lower().startswith("ref|"):
        head = ss.split("|", 1)[0].strip()
        if head and _RE_ACCESSION.fullmatch(head) and "." in head:
            return head
    return ""


def _fix_placeholder_ncbi_ids(ann: pd.DataFrame) -> pd.DataFrame:
    """BLAST tables sometimes leave best_hit/best_accession as 'ref'; use parsed accession.version."""
    if ann.empty or "best_hit" not in ann.columns:
        return ann
    out = ann.copy()
    if "best_accession" not in out.columns:
        out["best_accession"] = ""
    hits = []
    accs = []
    for _, r in out.iterrows():
        bh = str(r.get("best_hit", "") or "").strip().lower()
        ba = str(r.get("best_accession", "") or "").strip().lower()
        acc = _parse_versioned_accession(r)
        new_hit = str(r.get("best_hit", "") or "")
        new_acc = str(r.get("best_accession", "") or "")
        if acc and bh in ("", "ref", "nan") and new_hit.strip().lower() in ("", "ref"):
            new_hit = acc
        if acc and ba in ("", "ref", "nan") and new_acc.strip().lower() in ("", "ref"):
            new_acc = acc
        hits.append(new_hit)
        accs.append(new_acc)
    out["best_hit"] = hits
    out["best_accession"] = accs
    return out


def _enrich_ncbi_from_metadata(ann: pd.DataFrame, meta_map: dict[str, tuple[str, str]]) -> pd.DataFrame:
    if ann.empty or not meta_map:
        return ann
    out = ann.copy()
    if "best_species" not in out.columns:
        out["best_species"] = ""
    if "best_title" not in out.columns:
        out["best_title"] = ""
    species: list[str] = []
    titles: list[str] = []
    for _, r in out.iterrows():
        acc = _parse_versioned_accession(r)
        if acc in meta_map:
            o, d = meta_map[acc]
            species.append(o)
            titles.append(d)
        else:
            bs, bt = r.get("best_species", ""), r.get("best_title", "")
            species.append("" if pd.isna(bs) else str(bs))
            titles.append("" if pd.isna(bt) else str(bt))
    out["best_species"] = species
    out["best_title"] = titles
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Summarize cluster consensus annotations (NCBI + oligo best hits).")
    ap.add_argument("--project", required=True)
    ap.add_argument("--clusters-tsv", required=True)
    ap.add_argument("--ncbi-best-tsv", required=True)
    ap.add_argument("--oligo-best-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    project = str(args.project)
    cfg = _load_config()
    enabled = bool(resolve_cfg_value(cfg, project=project, path=["post_tarean_params", "graph_report", "enabled"]))
    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    if not enabled:
        out.write_text("")
        return 0

    clusters = pd.read_csv(Path(args.clusters_tsv), sep="\t") if Path(args.clusters_tsv).stat().st_size else pd.DataFrame()
    ncbi = pd.read_csv(Path(args.ncbi_best_tsv), sep="\t") if Path(args.ncbi_best_tsv).stat().st_size else pd.DataFrame()
    oligo = pd.read_csv(Path(args.oligo_best_tsv), sep="\t") if Path(args.oligo_best_tsv).stat().st_size else pd.DataFrame()

    # ncbi/oligo best tables are per-query (qseqid). Here qseqid == cluster_id (from align script).
    keep_ncbi = [c for c in ncbi.columns if c.startswith("best_")] + ["qseqid"]
    keep_oligo = [c for c in oligo.columns if c.startswith("oligo_") or c in ("qseqid", "oligo_fitting")]
    ncbi2 = ncbi[keep_ncbi].copy() if not ncbi.empty else pd.DataFrame(columns=["qseqid"])
    oligo2 = oligo[keep_oligo].copy() if not oligo.empty else pd.DataFrame(columns=["qseqid"])

    # Cluster sizes
    if not clusters.empty:
        sizes = clusters.groupby("cluster_id", as_index=False).size().rename(columns={"size": "n_repeats"})
    else:
        sizes = pd.DataFrame(columns=["cluster_id", "n_repeats"])

    ann = sizes.rename(columns={"cluster_id": "qseqid"}).merge(ncbi2, on="qseqid", how="left").merge(oligo2, on="qseqid", how="left")
    ann = ann.rename(columns={"qseqid": "cluster_id"})
    ann = _fix_placeholder_ncbi_ids(ann)

    meta_csv = _pick_newest_ncbi_metadata_csv(project)
    if meta_csv is not None:
        meta_map = _load_ncbi_metadata_map(meta_csv)
        ann = _enrich_ncbi_from_metadata(ann, meta_map)

    ann.to_csv(out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

