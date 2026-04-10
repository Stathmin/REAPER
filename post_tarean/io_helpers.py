#!/usr/bin/env python3
"""
Small IO/paths helpers for the post-TAREAN pipeline.

These utilities are intentionally lightweight and side-effect free so they can
be reused by both the legacy `RepeatAnalyzer.run_full_pipeline` method and the
new step-based orchestrator.
"""

from pathlib import Path
from typing import Any, List, Tuple, Optional
import re


def _infer_project_id(config_manager: Any, config: Any) -> Optional[str]:
    project_id = getattr(config, "project_id", None) or getattr(config_manager, "project_id", None)
    if project_id:
        return project_id

    # Last resort: if exactly one project exists in the global config, use it.
    try:
        import os
        import yaml  # type: ignore[import-not-found]

        global_cfg_path = os.environ.get("REPORTR_GLOBAL_CONFIG", "projects/global_config.yaml")
        with open(global_cfg_path, "r") as f:
            global_cfg = yaml.safe_load(f) or {}
        projects = (global_cfg.get("projects", {}) or {})
        if len(projects) == 1:
            return next(iter(projects.keys()))
    except Exception:
        return None

    return None


def resolve_tarean_path(config_manager: Any, config: Any, index: str) -> Tuple[Optional[Path], List[Path]]:
    """
    Locate the TAREAN output directory for a given sample or comparative index.

    Returns a tuple of (resolved_path, tried_paths). If no path exists, the
    resolved_path will be ``None`` and the caller is expected to handle the
    error condition.
    """
    project_id = _infer_project_id(config_manager, config)
    if not project_id:
        # Without a project id, we can only try the legacy relative-path layout.
        # Callers should treat this as a likely misconfiguration.
        candidates = [Path(f"./{index}")]
        for path in candidates:
            if path.exists():
                return path, candidates
        return None, candidates

    candidates: List[Path] = [
        Path(f"./{index}"),  # Original expected path
        Path(f"projects/{project_id}/samples/{index}/tarean"),
        # Comparative analyses (seqclust `-P` mode) store outputs here.
        Path(f"projects/{project_id}/comparative/{index}/tarean"),
        Path(f"../projects/{project_id}/samples/{index}/tarean"),
        Path(f"../projects/{project_id}/comparative/{index}/tarean"),
        Path(f"../../projects/{project_id}/samples/{index}/tarean"),
        Path(f"../../projects/{project_id}/comparative/{index}/tarean"),
    ]

    # Iterative layout (modern Snakemake): choose the newest *completed* iteration.
    # Solo: projects/<p>/samples/<s>/tareanN/...
    # Comparative: projects/<p>/comparative/<a>/tareanN/tarean/...
    try:
        sample_base = Path(f"projects/{project_id}/samples/{index}")
        if sample_base.exists():
            iters_done = []
            iters_partial = []
            for p in sample_base.glob("tarean*/"):
                m = re.match(r"^tarean(\d+)$", p.name)
                if not m:
                    continue
                it = int(m.group(1))
                if (p / "tarean.done").exists():
                    iters_done.append((it, p))
                elif (p / "CLUSTER_TABLE.csv").exists():
                    # Fallback for older/partial runs: keep as a last resort, but
                    # prefer a true completion token when present.
                    iters_partial.append((it, p))
            if iters_done:
                _, p = sorted(iters_done, key=lambda t: t[0], reverse=True)[0]
                candidates.insert(1, p)
            elif iters_partial:
                _, p = sorted(iters_partial, key=lambda t: t[0], reverse=True)[0]
                candidates.insert(1, p)
    except Exception:
        pass

    try:
        comp_base = Path(f"projects/{project_id}/comparative/{index}")
        if comp_base.exists():
            iters_done = []
            iters_partial = []
            for p in comp_base.glob("tarean*/"):
                m = re.match(r"^tarean(\d+)$", p.name)
                if not m:
                    continue
                it = int(m.group(1))
                # Iterative comparative rules write core outputs under tareanN/tarean/.
                core = p / "tarean"
                if (p / "COMPARATIVE_TAREAN_COMPLETE").exists():
                    iters_done.append((it, core))
                elif (core / "CLUSTER_TABLE.csv").exists():
                    iters_partial.append((it, core))
            if iters_done:
                _, core = sorted(iters_done, key=lambda t: t[0], reverse=True)[0]
                candidates.insert(1, core)
            elif iters_partial:
                _, core = sorted(iters_partial, key=lambda t: t[0], reverse=True)[0]
                candidates.insert(1, core)
    except Exception:
        pass

    for path in candidates:
        if path.exists():
            return path, candidates

    # Nothing found; let the caller decide how to report this.
    return None, candidates

