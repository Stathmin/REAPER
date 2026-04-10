"""
Layered configuration resolver for RepOrtR.

Design goals:
- One DRY implementation used by Snakemake (Snakefile/workflows) and Python scripts.
- Deep-merge semantics for dicts: global → project → comparative_analysis → sample.
- Strict comparative mode: ignore sample overrides for operational parameters.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping


def deep_merge_dicts(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge override into base (without mutating inputs)."""
    out: dict[str, Any] = deepcopy(dict(base))
    for k, v in override.items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict):
            out[k] = deep_merge_dicts(out[k], v)
        else:
            out[k] = deepcopy(v)
    return out


def _get_path(root: Mapping[str, Any], path: list[str]) -> Any:
    cur: Any = root
    for k in path:
        if not isinstance(cur, Mapping) or k not in cur:
            raise KeyError(k)
        cur = cur[k]
    return cur


def resolve_config(
    config: Mapping[str, Any],
    *,
    project: str,
    sample: str | None = None,
    comparative_analysis: str | None = None,
    allow_sample_overrides: bool = True,
) -> dict[str, Any]:
    """
    Return a resolved config dict for a given context.

    Layers:
      global
      projects.<project>
      projects.<project>.comparative_analyses.<comparative_analysis> (if provided)
      projects.<project>.samples.<sample> (if provided and allow_sample_overrides)
    """
    global_cfg = dict(config.get("global", {}) or {})
    proj_cfg = dict((config.get("projects", {}) or {}).get(project, {}) or {})
    resolved = deep_merge_dicts(global_cfg, proj_cfg)

    if comparative_analysis is not None:
        ca = (proj_cfg.get("comparative_analyses", {}) or {}).get(comparative_analysis, {}) or {}
        resolved = deep_merge_dicts(resolved, dict(ca))

    if sample is not None and allow_sample_overrides:
        sm = (proj_cfg.get("samples", {}) or {}).get(sample, {}) or {}
        resolved = deep_merge_dicts(resolved, dict(sm))

    return resolved


def resolve_value(
    config: Mapping[str, Any],
    *,
    project: str,
    path: list[str],
    sample: str | None = None,
    comparative_analysis: str | None = None,
    allow_sample_overrides: bool = True,
) -> Any:
    """Resolve a scalar/list/dict from the resolved context."""
    resolved = resolve_config(
        config,
        project=project,
        sample=sample,
        comparative_analysis=comparative_analysis,
        allow_sample_overrides=allow_sample_overrides,
    )
    try:
        return _get_path(resolved, path)
    except KeyError as e:
        ctx = f"project={project}"
        if comparative_analysis is not None:
            ctx += f" comparative_analysis={comparative_analysis}"
        if sample is not None and allow_sample_overrides:
            ctx += f" sample={sample}"
        raise KeyError(f"Missing config key at path={'.'.join(path)} ({ctx})") from e


def sample_metadata(config: Mapping[str, Any], *, project: str, sample: str) -> dict[str, Any]:
    """Access per-sample metadata (never merged)."""
    return (
        dict(((config.get("projects", {}) or {}).get(project, {}) or {}).get("samples", {}) or {})
        .get(sample, {})
        or {}
    )

