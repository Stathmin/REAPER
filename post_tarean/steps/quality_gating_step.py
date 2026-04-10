"""
Heuristic quality gating / confidence scoring for post-TAREAN results.

This step consumes the merged summary table produced by ``SummaryStep`` and
emits a per-sample JSON file with per-cluster confidence metrics and an
action-oriented recommendation.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple


from .base import PostTareanContext, PostTareanStep


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        if value is None:
            return default
        if isinstance(value, str) and not value.strip():
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def _clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def _pick(row: Mapping[str, Any], keys: List[str]) -> Any:
    for k in keys:
        if k in row:
            return row[k]
    return None


def _annotation_prior(annotation: str) -> Tuple[float, List[str]]:
    a = (annotation or "").lower()
    reasons: List[str] = []

    # Common labels observed in pipeline reports.
    if "putative satellites (high confidence)" in a:
        reasons.append("TAREAN indicates high-confidence satellite")
        return 0.95, reasons
    if "putative satellites (low confidence)" in a:
        reasons.append("TAREAN indicates low-confidence satellite")
        return 0.75, reasons
    if "satellite" in a:
        reasons.append("Satellite-like annotation")
        return 0.80, reasons
    if "ltr" in a:
        reasons.append("LTR element annotation (not a satellite)")
        return 0.55, reasons
    if "mobile" in a or "class_i" in a or "class_ii" in a:
        reasons.append("Mobile element annotation")
        return 0.50, reasons
    if "unclassified" in a:
        reasons.append("Unclassified annotation")
        return 0.35, reasons

    return 0.45, reasons


def _evalue_score(evalue: float, max_evalue: float) -> float:
    """
    Map e-value to [0..1] where <= max_evalue is strong.
    Uses a log-scale decay above the threshold.
    """
    if evalue <= 0.0:
        return 1.0
    if evalue <= max_evalue:
        return 1.0
    # Penalize by log10 ratio; one order of magnitude over threshold -> subtract ~0.2
    ratio = evalue / max_evalue
    penalty = 0.2 * math.log10(max(1.0, ratio))
    return _clamp(1.0 - penalty)


def _pident_score(pident: float, min_pident: float) -> float:
    if pident <= 0.0:
        return 0.0
    if pident >= 100.0:
        return 1.0
    if pident <= min_pident:
        return 0.0
    return _clamp((pident - min_pident) / (100.0 - min_pident))


def _coverage_score(coverage: float) -> float:
    # Simple linear map; coverage is already 0..1 in the pipeline.
    return _clamp(coverage)


def _recommendation(score: float, cutoffs: Mapping[str, float]) -> str:
    # Order matters.
    if score >= float(cutoffs.get("publish", 0.85)):
        return "PUBLISH"
    if score >= float(cutoffs.get("validate", 0.65)):
        return "VALIDATE"
    if score >= float(cutoffs.get("refine", 0.40)):
        return "REFINE"
    return "DISCARD"


class QualityGatingStep(PostTareanStep):
    name = "quality_gating"

    def run(self, analyzer: Any, context: PostTareanContext, state: Dict[str, Any]) -> Dict[str, Any]:
        summary = state.get("summary_data")
        if summary is None:
            # Summary is expected; if missing, do nothing but keep pipeline running.
            return state

        # Pull configuration with safe fallbacks (P1.2 will formalize this).
        gating_cfg = getattr(analyzer.config, "quality_gating", None)
        enabled = getattr(gating_cfg, "enabled", True) if gating_cfg is not None else True
        if not enabled:
            return state

        w_annotation = float(getattr(gating_cfg, "w_annotation", 0.35) if gating_cfg is not None else 0.35)
        w_abundance = float(getattr(gating_cfg, "w_abundance", 0.25) if gating_cfg is not None else 0.25)
        w_blast = float(getattr(gating_cfg, "w_blast", 0.40) if gating_cfg is not None else 0.40)

        min_pident = float(
            getattr(gating_cfg, "min_pident", None) if gating_cfg is not None else None
            or getattr(getattr(analyzer.config, "blast", None), "pident_threshold", 60.0)
        )
        max_evalue = float(
            getattr(gating_cfg, "max_evalue", None) if gating_cfg is not None else None
            or getattr(getattr(analyzer.config, "blast", None), "evalue_max_olig", 1e-3)
        )

        cutoffs = {
            "publish": float(getattr(gating_cfg, "publish", 0.85) if gating_cfg is not None else 0.85),
            "validate": float(getattr(gating_cfg, "validate", 0.65) if gating_cfg is not None else 0.65),
            "refine": float(getattr(gating_cfg, "refine", 0.40) if gating_cfg is not None else 0.40),
        }

        clusters: List[Dict[str, Any]] = []

        # Iterate rows robustly: support pandas DataFrame without importing pandas here.
        records = summary.to_dict(orient="records") if hasattr(summary, "to_dict") else []
        for row in records:
            cluster_id = _pick(row, ["Cluster", "cluster", "qseqid"]) or "unknown"
            annotation = str(_pick(row, ["TAREAN_annotation", "tarean_annotation", "annotation"]) or "")
            percent_genome = _safe_float(_pick(row, ["size, %", "size,%", "percent_genome", "Percent_genome"]), default=0.0)

            reasons: List[str] = []

            ann_score, ann_reasons = _annotation_prior(annotation)
            reasons.extend(ann_reasons)

            # Abundance: saturate at ~5% of genome for score=1.0 (tunable later).
            abundance_score = _clamp(percent_genome / 5.0)
            if percent_genome >= 1.0:
                reasons.append("High abundance cluster")
            elif percent_genome <= 0.1:
                reasons.append("Very low abundance cluster")

            # BLAST-derived score (if available)
            best_hit = str(_pick(row, ["best_hit", "Best_hit"]) or "")
            has_hit = best_hit and "no significant" not in best_hit.lower()
            evalue = _safe_float(_pick(row, ["best_evalue", "Best_evalue", "evalue"]), default=float("inf"))
            pident = _safe_float(_pick(row, ["best_pident", "Best_pident", "pident"]), default=0.0)
            coverage = _safe_float(_pick(row, ["coverage", "Coverage"]), default=0.0)

            if has_hit:
                ev = _evalue_score(evalue, max_evalue=max_evalue)
                pid = _pident_score(pident, min_pident=min_pident)
                cov = _coverage_score(coverage)
                blast_score = _clamp(0.45 * ev + 0.35 * pid + 0.20 * cov)
                if evalue <= max_evalue:
                    reasons.append("Strong BLAST e-value")
                if pident >= min_pident:
                    reasons.append("Meets identity threshold")
                if coverage >= 0.8:
                    reasons.append("High coverage alignment")
            else:
                blast_score = 0.0
                reasons.append("No significant BLAST hits")

            # Composite score
            denom = max(1e-9, (w_annotation + w_abundance + w_blast))
            composite = _clamp((w_annotation * ann_score + w_abundance * abundance_score + w_blast * blast_score) / denom)

            rec = _recommendation(composite, cutoffs=cutoffs)

            clusters.append(
                {
                    "cluster_id": str(cluster_id),
                    "percent_genome": percent_genome,
                    "metrics": {
                        "tarean_annotation_score": ann_score,
                        "abundance_score": abundance_score,
                        "blast_score": blast_score,
                    },
                    "composite_confidence": composite,
                    "recommendation": rec,
                    "reasons": reasons[:8],
                }
            )

        # Sort by abundance desc for convenience.
        clusters.sort(key=lambda x: float(x.get("percent_genome") or 0.0), reverse=True)

        payload: Dict[str, Any] = {
            "project": context.project_id,
            "sample": context.subject,
            "generated_at": datetime.now(tz=timezone.utc).isoformat(),
            "config": {
                "w_annotation": w_annotation,
                "w_abundance": w_abundance,
                "w_blast": w_blast,
                "min_pident": min_pident,
                "max_evalue": max_evalue,
                "cutoffs": cutoffs,
            },
            "clusters": clusters,
        }

        out_path = Path(context.output_dir) / f"{context.subject}_confidence.json"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")

        new_state = dict(state)
        new_state["confidence_json"] = str(out_path)
        new_state["confidence_payload"] = payload
        return new_state

