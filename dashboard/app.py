from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

import streamlit as st

from dashboard.lib.progress import (
    GLOBAL_CONFIG_PATH,
    LOGS_DIR,
    compute_overall_progress,
    find_progress_files,
    format_eta,
    guess_project_and_sample_from_stem,
    human_time,
    list_projects,
    list_samples_for_project,
    load_global_config,
    load_progress,
    path_exists,
    progress_mtime_seconds,
    sample_status,
    tail_text_file,
    tarean_done_path,
)


st.set_page_config(page_title="RepOrtR Dashboard", layout="wide")


@dataclass(frozen=True)
class UiSample:
    project_id: str
    sample_id: str
    progress_path: Optional[Path]
    progress: Dict[str, Any]


def _cache_dir() -> Path:
    cfg = load_global_config()
    global_cfg = (cfg.get("global") or {}) if isinstance(cfg, dict) else {}
    raw = global_cfg.get("cache_dir") or ".snakemake/cache"
    return Path(str(raw))


def _dir_size_bytes(root: Path, max_files: int = 200_000) -> int:
    total = 0
    count = 0
    if not root.exists():
        return 0
    for base, _dirs, files in os.walk(root):
        for name in files:
            count += 1
            if count > max_files:
                return total
            p = Path(base) / name
            try:
                total += p.stat().st_size
            except FileNotFoundError:
                continue
    return total


def _format_bytes(num: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    value = float(num)
    for u in units:
        if value < 1024.0 or u == units[-1]:
            return f"{value:,.2f} {u}"
        value /= 1024.0
    return f"{num} B"


def _load_samples() -> Dict[str, UiSample]:
    cfg = load_global_config()
    projects = list_projects(cfg)
    progress_files = find_progress_files()

    samples: Dict[str, UiSample] = {}

    # Index progress files by parsed project/sample (best-effort).
    progress_by_key: Dict[tuple[str, str], Path] = {}
    legacy_by_sample: Dict[str, Path] = {}
    for stem, path in progress_files.items():
        p, s = guess_project_and_sample_from_stem(stem)
        if p and s:
            progress_by_key[(p, s)] = path
        elif s:
            legacy_by_sample[s] = path

    for project_id in sorted(projects.keys()):
        project_samples = list_samples_for_project(cfg, project_id)
        for sample_id in sorted(project_samples.keys()):
            p_path = progress_by_key.get((project_id, sample_id)) or legacy_by_sample.get(sample_id)
            pdata = load_progress(p_path) if p_path else {}
            key = f"{project_id}/{sample_id}"
            samples[key] = UiSample(
                project_id=project_id,
                sample_id=sample_id,
                progress_path=p_path,
                progress=pdata,
            )

    return samples


def _sample_header(sample: UiSample) -> str:
    status = sample_status(
        project_id=sample.project_id,
        sample_id=sample.sample_id,
        progress_path=sample.progress_path,
    )
    label = f"{sample.project_id}/{sample.sample_id}"
    if status == "completed":
        return f"✅ {label}"
    if status == "running":
        return f"🔄 {label}"
    if status == "stale":
        return f"⚠️ {label}"
    return f"⏳ {label}"


def _render_sample_card(sample: UiSample, *, stale_after_seconds: int) -> None:
    pdata = sample.progress or {}

    status = sample_status(
        project_id=sample.project_id,
        sample_id=sample.sample_id,
        progress_path=sample.progress_path,
        stale_after_seconds=stale_after_seconds,
    )

    overall = compute_overall_progress(pdata) if pdata else 0.0
    current_stage = str(pdata.get("current_stage") or "unknown")
    eta = pdata.get("eta")
    resource = pdata.get("resource_usage") or {}
    cpu = resource.get("cpu_percent", 0.0)
    mem_mb = resource.get("memory_mb", 0.0)

    meta = []
    if sample.progress_path:
        mtime = progress_mtime_seconds(sample.progress_path)
        meta.append(f"progress_json: `{sample.progress_path}` (updated {human_time(mtime)})")
    done = tarean_done_path(sample.project_id, sample.sample_id)
    if done.exists():
        meta.append(f"tarean_done: `{done}`")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Status", status)
    c2.metric("Overall", f"{overall:.1f}%")
    c3.metric("Stage", current_stage.replace("_", " ").title())
    c4.metric("ETA", format_eta(eta) if eta is not None else "unknown")

    c5, c6 = st.columns(2)
    c5.metric("CPU%", f"{float(cpu):.1f}" if cpu is not None else "0.0")
    c6.metric("MEM (MB)", f"{float(mem_mb):.1f}" if mem_mb is not None else "0.0")

    if meta:
        st.caption(" | ".join(meta))

    stages = pdata.get("progress") if isinstance(pdata, dict) else None
    if isinstance(stages, dict) and stages:
        st.write("Stages")
        for stage_name, info in stages.items():
            if not isinstance(info, dict):
                continue
            prog = float(info.get("progress") or 0.0)
            st.progress(min(1.0, max(0.0, prog / 100.0)), text=f"{stage_name}: {prog:.1f}% ({info.get('status','?')})")
    else:
        st.info("No per-stage progress available yet (missing or unreadable progress JSON).")

    # Heuristic failure / stale warning.
    if status == "stale" and not done.exists():
        st.warning(
            "Progress looks stale (no recent updates) and `tarean.done` is missing. "
            "This may indicate a failed or stuck run."
        )

    # Confidence / quality gating (optional)
    conf_path = Path(
        f"projects/{sample.project_id}/samples/{sample.sample_id}/post_tarean/{sample.sample_id}_confidence.json"
    )
    if conf_path.exists():
        try:
            conf = json.loads(conf_path.read_text())
        except Exception:
            conf = None

        if isinstance(conf, dict):
            st.write("Quality gating")
            clusters = conf.get("clusters") or []
            if isinstance(clusters, list) and clusters:
                # Distribution of recommendations
                counts: Dict[str, int] = {}
                for c in clusters:
                    if not isinstance(c, dict):
                        continue
                    rec = str(c.get("recommendation") or "UNKNOWN")
                    counts[rec] = counts.get(rec, 0) + 1

                cols = st.columns(4)
                for idx, key in enumerate(["PUBLISH", "VALIDATE", "REFINE", "DISCARD"]):
                    cols[idx].metric(key, str(counts.get(key, 0)))

                # Show top clusters by abundance
                top = []
                for c in clusters[:20]:
                    if not isinstance(c, dict):
                        continue
                    top.append(
                        {
                            "cluster_id": c.get("cluster_id"),
                            "percent_genome": c.get("percent_genome"),
                            "confidence": c.get("composite_confidence"),
                            "recommendation": c.get("recommendation"),
                        }
                    )
                st.dataframe(top, use_container_width=True)
                st.caption(f"Loaded `{conf_path}`")
            else:
                st.caption(f"No clusters in `{conf_path}`")


def _render_logs_panel(sample: UiSample, *, max_lines: int) -> None:
    st.write("Logs")
    candidates = [
        LOGS_DIR / f"seqclust_{sample.project_id}_{sample.sample_id}.log",
        LOGS_DIR / "snakemake.log",
        LOGS_DIR / f"progress_tracking_{sample.project_id}_{sample.sample_id}.log",
    ]
    available = [p for p in candidates if p.exists()]
    if not available:
        st.caption("No known log files found for this sample yet.")
        return

    selected = st.selectbox(
        "Select log file",
        options=available,
        format_func=lambda p: str(p),
        key=f"log_select_{sample.project_id}_{sample.sample_id}",
    )
    content = tail_text_file(Path(selected), max_lines=max_lines)
    if not content:
        st.caption("Log file is empty (or missing).")
        return
    st.code(content)


def _render_input_path_warnings(project_cfg: Dict[str, Any], sample_id: str) -> None:
    samples = project_cfg.get("samples") or {}
    sample_cfg = samples.get(sample_id) if isinstance(samples, dict) else None
    if not isinstance(sample_cfg, dict):
        return

    r1 = sample_cfg.get("r1_path")
    r2 = sample_cfg.get("r2_path")
    missing = []
    for label, p in (("r1_path", r1), ("r2_path", r2)):
        if isinstance(p, str) and p and not path_exists(p):
            missing.append((label, p))
    if missing:
        st.warning(
            "Some input paths from `projects/global_config.yaml` are not reachable from the current working directory "
            "(common inside containers unless you bind-mount those locations)."
        )
        for label, p in missing:
            st.code(f"{label}: {p}")


def main() -> None:
    st.title("RepOrtR Dashboard")
    st.caption(f"Config: `{GLOBAL_CONFIG_PATH}` | Logs: `{LOGS_DIR}`")

    with st.sidebar:
        st.header("Controls")
        if st.button("Refresh now"):
            try:
                st.rerun()
            except AttributeError:
                st.experimental_rerun()

        st.caption(f"Last refresh: {datetime.now(tz=timezone.utc).astimezone().strftime('%Y-%m-%d %H:%M:%S')}")
        stale_after_minutes = st.number_input("Stale threshold (minutes)", min_value=1, max_value=240, value=5)
        log_tail_lines = st.number_input("Log tail lines", min_value=50, max_value=2000, value=300, step=50)
        show_cache = st.checkbox("Show cache size (can be slow)", value=False)

        st.divider()
        st.header("About")
        st.write(
            "This dashboard is read-only. It renders progress from `logs/progress_*.json` "
            "generated by `workflows/progress_tracking.smk`."
        )

    cfg = load_global_config()
    projects = list_projects(cfg)
    samples = _load_samples()

    top_left, top_right = st.columns([2, 1])
    with top_left:
        st.subheader("Projects")
        if not projects:
            st.error("No projects found in `projects/global_config.yaml`.")
        else:
            for project_id in sorted(projects.keys()):
                project_cfg = projects[project_id] if isinstance(projects, dict) else {}
                project_samples = list_samples_for_project(cfg, project_id)
                total = len(project_samples)
                running = 0
                completed = 0
                stale = 0
                unknown = 0
                for s_id in project_samples.keys():
                    key = f"{project_id}/{s_id}"
                    ui = samples.get(key)
                    p_path = ui.progress_path if ui else None
                    s = sample_status(
                        project_id=project_id,
                        sample_id=s_id,
                        progress_path=p_path,
                        stale_after_seconds=int(stale_after_minutes * 60),
                    )
                    if s == "running":
                        running += 1
                    elif s == "completed":
                        completed += 1
                    elif s == "stale":
                        stale += 1
                    else:
                        unknown += 1

                st.write(
                    f"**{project_id}** — samples: {total} | "
                    f"✅ {completed} | 🔄 {running} | ⚠️ {stale} | ⏳ {unknown}"
                )
                _render_input_path_warnings(project_cfg if isinstance(project_cfg, dict) else {}, next(iter(project_samples.keys()), ""))

    with top_right:
        st.subheader("Environment")
        st.write(f"cwd: `{Path.cwd()}`")
        st.write(f"cache_dir: `{_cache_dir()}`")
        if show_cache:
            cache_dir = _cache_dir()
            size = _dir_size_bytes(cache_dir)
            st.metric("Cache size", _format_bytes(size))
            st.caption("If this looks too slow, uncheck the option in the sidebar.")

    st.divider()
    st.subheader("Samples")

    if not samples:
        st.info("No samples found (no projects or no samples configured).")
        return

    # Filter controls
    all_keys = list(samples.keys())
    default_project = all_keys[0].split("/", 1)[0] if all_keys else ""
    project_filter = st.selectbox("Project", options=sorted(projects.keys()) if projects else [default_project])
    status_filter = st.multiselect(
        "Status",
        options=["running", "completed", "stale", "unknown"],
        default=["running", "stale", "unknown", "completed"],
    )

    stale_after_seconds = int(stale_after_minutes * 60)
    for key, sample in samples.items():
        if sample.project_id != project_filter:
            continue
        s = sample_status(
            project_id=sample.project_id,
            sample_id=sample.sample_id,
            progress_path=sample.progress_path,
            stale_after_seconds=stale_after_seconds,
        )
        if s not in status_filter:
            continue

        with st.expander(_sample_header(sample), expanded=(s in ("running", "stale"))):
            left, right = st.columns([2, 1])
            with left:
                _render_sample_card(sample, stale_after_seconds=stale_after_seconds)
            with right:
                _render_logs_panel(sample, max_lines=int(log_tail_lines))


if __name__ == "__main__":
    main()

