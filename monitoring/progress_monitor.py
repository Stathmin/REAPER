#!/usr/bin/env python3
"""
Config-driven TAREAN progress monitor for RepOrtR.

This tool reads the JSON progress files produced by the Snakemake
`track_seqclust_progress` rule (see `workflows/progress_tracking.smk`)
and renders a simple terminal dashboard.

It is intended to be the single, optional monitoring entry point, in
preference to older ad-hoc scripts such as `monitoring/holy_monitor.py` and
`monitoring/monitor_pipeline.py`.
"""

import glob
import json
import os
import time
from datetime import datetime
from typing import Any, Dict

import yaml


GLOBAL_CONFIG_PATH = "projects/global_config.yaml"

def load_log_dir() -> str:
    """Resolve `global.log_dir` from projects/global_config.yaml."""
    try:
        with open(GLOBAL_CONFIG_PATH, "r") as f:
            cfg = yaml.safe_load(f) or {}
        return (cfg.get("global") or {}).get("log_dir", "logs")
    except FileNotFoundError:
        return "logs"
    except Exception:
        return "logs"


def load_progress_config() -> Dict[str, Any]:
    """Load `global.progress_tracking` configuration."""
    try:
        with open(GLOBAL_CONFIG_PATH, "r") as f:
            cfg = yaml.safe_load(f) or {}
    except FileNotFoundError:
        return {}

    return (cfg.get("global") or {}).get("progress_tracking") or {}


def find_progress_files() -> Dict[str, str]:
    """Return mapping {sample_id: path_to_progress_json}."""
    log_dir = load_log_dir()
    files = glob.glob(os.path.join(log_dir, "progress_*.json"))
    mapping: Dict[str, str] = {}
    for path in files:
        name = os.path.basename(path)
        # progress_<project>_<sample>.json OR progress_<sample>.json
        stem = os.path.splitext(name)[0].replace("progress_", "", 1)
        mapping[stem] = path
    return mapping


def load_progress(path: str) -> Dict[str, Any]:
    """Load a single progress JSON file, returning {} on error."""
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def format_progress_bar(progress: float, width: int = 20) -> str:
    """Return a simple ASCII progress bar."""
    progress = max(0.0, min(100.0, progress))
    filled = int(width * progress / 100.0)
    bar = "█" * filled + "░" * (width - filled)
    return f"[{bar}] {progress:5.1f}%"


def compute_overall_progress(progress_data: Dict[str, Any]) -> float:
    """Compute overall pipeline progress from per-stage data."""
    stages = progress_data.get("progress") or {}
    if not stages:
        return 0.0

    total = 0.0
    count = 0
    for stage_info in stages.values():
        status = stage_info.get("status")
        if status == "completed":
            total += 100.0
        else:
            total += float(stage_info.get("progress") or 0.0)
        count += 1

    return total / max(1, count)


def format_eta(seconds: Any) -> str:
    """Format ETA seconds into a human-readable string."""
    try:
        value = int(seconds)
    except (TypeError, ValueError):
        return "unknown"

    if value <= 0:
        return "0s"

    minutes, sec = divmod(value, 60)
    hours, minutes = divmod(minutes, 60)
    if hours:
        return f"{hours}h {minutes}m {sec}s"
    if minutes:
        return f"{minutes}m {sec}s"
    return f"{sec}s"


def render_sample_block(progress_data: Dict[str, Any]) -> str:
    """Render a single-sample progress block as text."""
    sample = progress_data.get("sample") or "unknown"
    project = progress_data.get("project") or "-"
    stages = progress_data.get("progress") or {}
    current_stage = progress_data.get("current_stage") or "initialization"
    resource = progress_data.get("resource_usage") or {}
    eta_seconds = progress_data.get("eta")

    overall = compute_overall_progress(progress_data)

    lines = []
    header = f"TAREAN Progress: {project}/{sample}"
    lines.append(header)
    lines.append("=" * len(header))
    lines.append(f"Overall: {format_progress_bar(overall)}")
    lines.append(
        f"Stage : {current_stage.replace('_', ' ').title()} | "
        f"CPU {resource.get('cpu_percent', 0):4.1f}% | "
        f"MEM {resource.get('memory_mb', 0):6.1f} MB"
    )

    if eta_seconds is not None:
        lines.append(f"ETA  : {format_eta(eta_seconds)}")

    lines.append("")
    lines.append("Stages:")
    for stage_name, info in stages.items():
        status = info.get("status")
        progress = float(info.get("progress") or 0.0)
        icon = "✅" if status == "completed" else "🔄" if status == "running" else "⏳"
        label = stage_name.replace("_", " ").title()
        bar = format_progress_bar(progress, width=10)
        lines.append(f"  {icon} {label:<18} {bar}")

    return "\n".join(lines)


def run_dashboard() -> None:
    """Main monitoring loop."""
    cfg = load_progress_config()
    if not cfg.get("enabled", False):
        print(
            "Progress tracking is disabled in `projects/global_config.yaml` "
            "(global.progress_tracking.enabled = false)."
        )
        return

    update_interval = int(cfg.get("display_interval") or cfg.get("update_interval") or 10)

    try:
        while True:
            os.system("clear")
            now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            print(f"RepOrtR TAREAN Progress Dashboard  {now}")
            print("=" * 80)

            mapping = find_progress_files()
            if not mapping:
                print(
                    f"No progress files found in `{load_log_dir()}/`. "
                    "Is `track_seqclust_progress` running?"
                )
            else:
                first = True
                for key, path in sorted(mapping.items()):
                    data = load_progress(path)
                    if not data:
                        continue
                    if not first:
                        print("\n" + "-" * 80 + "\n")
                    print(render_sample_block(data))
                    first = False

            print("\n(Refreshes every", update_interval, "seconds; Ctrl+C to exit.)")
            time.sleep(update_interval)
    except KeyboardInterrupt:
        print("\nStopping progress monitor.")


def main() -> None:
    run_dashboard()


if __name__ == "__main__":
    main()

