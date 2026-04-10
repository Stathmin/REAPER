#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional


try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover
    psutil = None  # type: ignore


@dataclass(frozen=True)
class Paths:
    progress_json: Path
    progress_txt: Path
    seqclust_log: Path
    runner_log: Path


def _log(line: str, runner_log: Path) -> None:
    runner_log.parent.mkdir(parents=True, exist_ok=True)
    with runner_log.open("a", encoding="utf-8") as f:
        f.write(f"{datetime.now().isoformat(timespec='seconds')} {line}\n")


def _initial_progress(project: str, sample: str, stages: Dict[str, Any]) -> Dict[str, Any]:
    return {
        "sample": sample,
        "project": project,
        "start_time": int(time.time()),
        "current_stage": "initialization",
        "stages": stages,
        "progress": {
            name: {"status": "pending", "progress": 0, "start_time": None, "end_time": None}
            for name in stages.keys()
        },
        "resource_usage": {"cpu_percent": 0, "memory_percent": 0, "memory_mb": 0},
        "eta": None,
    }


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except FileNotFoundError:
        return ""


def _seqclust_running(sample: str) -> bool:
    if psutil is None:
        return True  # can't reliably detect; avoid falsely marking complete
    needle = sample
    for proc in psutil.process_iter(["cmdline"]):
        try:
            cmd = " ".join(proc.info.get("cmdline") or [])
            if "seqclust" in cmd and needle in cmd:
                return True
        except Exception:
            continue
    return False


def _update_resources(progress: Dict[str, Any], sample: str) -> None:
    if psutil is None:
        return
    for proc in psutil.process_iter(["cmdline", "cpu_percent", "memory_percent"]):
        try:
            cmd = " ".join(proc.info.get("cmdline") or [])
            if "seqclust" in cmd and sample in cmd:
                p = psutil.Process(proc.pid)
                progress["resource_usage"]["cpu_percent"] = float(proc.info.get("cpu_percent") or 0.0)
                progress["resource_usage"]["memory_percent"] = float(proc.info.get("memory_percent") or 0.0)
                progress["resource_usage"]["memory_mb"] = float(p.memory_info().rss) / 1024 / 1024
                return
        except Exception:
            continue


def _detect_stages(progress: Dict[str, Any], stages: Dict[str, Any], log_content: str) -> None:
    for stage_name, stage_info in stages.items():
        stage_data = progress["progress"].get(stage_name)
        if not stage_data:
            continue

        # start
        for pattern in stage_info.get("patterns", []):
            if re.search(pattern, log_content, re.IGNORECASE):
                if stage_data["status"] == "pending":
                    stage_data["status"] = "running"
                    stage_data["start_time"] = int(time.time())
                    progress["current_stage"] = stage_name
                break

        # completion (heuristic: detect next stage patterns)
        if stage_data["status"] == "running":
            next_stages = list(stages.keys())
            try:
                idx = next_stages.index(stage_name)
            except ValueError:
                idx = -1
            if 0 <= idx < len(next_stages) - 1:
                next_stage = next_stages[idx + 1]
                for pattern in stages.get(next_stage, {}).get("patterns", []):
                    if re.search(pattern, log_content, re.IGNORECASE):
                        stage_data["status"] = "completed"
                        stage_data["end_time"] = int(time.time())
                        stage_data["progress"] = 100
                        break


def _estimate_progress_and_eta(progress: Dict[str, Any], stages: Dict[str, Any]) -> None:
    running_stage: Optional[str] = None
    for name, data in progress["progress"].items():
        if data.get("status") == "running":
            running_stage = name
            break

    if running_stage:
        stage_data = progress["progress"][running_stage]
        start_time = stage_data.get("start_time")
        if start_time:
            elapsed = int(time.time()) - int(start_time)
            est = int(stages.get(running_stage, {}).get("estimated_duration") or 0)
            if est > 0 and elapsed < est:
                stage_progress = min(80.0, (elapsed / est) * 100.0)
            else:
                stage_progress = min(95.0, 80.0 + max(0, elapsed - est) / 300.0)
            stage_data["progress"] = stage_progress

            remaining = max(0, est - elapsed)
            for s in stages.keys():
                if progress["progress"][s]["status"] == "pending":
                    remaining += int(stages[s].get("estimated_duration") or 0)
            progress["eta"] = remaining


def _render_txt(progress: Dict[str, Any]) -> str:
    sample = progress.get("sample") or "unknown"
    project = progress.get("project") or "-"
    current_stage = (progress.get("current_stage") or "initialization").replace("_", " ").title()
    resource = progress.get("resource_usage") or {}
    eta = progress.get("eta")

    # overall
    stages = progress.get("progress") or {}
    if stages:
        total = 0.0
        count = 0
        for info in stages.values():
            if info.get("status") == "completed":
                total += 100.0
            else:
                total += float(info.get("progress") or 0.0)
            count += 1
        overall = min(95.0, total / max(1, count))
    else:
        overall = 0.0

    lines = []
    lines.append(f"🎯 TAREAN Progress: {project}/{sample}")
    lines.append("=" * 50)
    lines.append(f"📊 Overall Progress: {overall:.1f}%")
    lines.append(f"⏱️  Current Stage: {current_stage}")
    lines.append(f"💻 CPU Usage: {float(resource.get('cpu_percent') or 0):.1f}%")
    lines.append(f"💾 Memory Usage: {float(resource.get('memory_mb') or 0):.1f} MB")
    lines.append("")
    lines.append("📋 Stage Status:")
    for stage_name, stage_data in stages.items():
        status = stage_data.get("status")
        icon = "✅" if status == "completed" else "🔄" if status == "running" else "⏳"
        label = stage_name.replace("_", " ").title()
        p = float(stage_data.get("progress") or 0.0)
        bar = "█" * int(p / 10) + "░" * (10 - int(p / 10))
        lines.append(f"  {icon} {label}: [{bar}] {p:.1f}%")
    if eta:
        mins = int(eta) // 60
        secs = int(eta) % 60
        lines.append(f"\n⏰ Estimated Time Remaining: {mins}m {secs}s")
    lines.append(f"\n🕐 Last Updated: {datetime.now().strftime('%H:%M:%S')}")
    return "\n".join(lines) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description="Track seqclust progress and write JSON+TXT progress files.")
    ap.add_argument("--project", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--seqclust-log", required=True)
    ap.add_argument("--progress-json", required=True)
    ap.add_argument("--progress-txt", required=True)
    ap.add_argument("--runner-log", required=True)
    ap.add_argument("--stages-json", required=True, help="JSON dict of stage definitions")
    ap.add_argument("--interval", type=int, default=30)
    args = ap.parse_args()

    paths = Paths(
        progress_json=Path(args.progress_json),
        progress_txt=Path(args.progress_txt),
        seqclust_log=Path(args.seqclust_log),
        runner_log=Path(args.runner_log),
    )

    stages = json.loads(args.stages_json)
    progress = _initial_progress(args.project, args.sample, stages)

    paths.progress_json.parent.mkdir(parents=True, exist_ok=True)
    paths.progress_txt.parent.mkdir(parents=True, exist_ok=True)

    _log(f"Starting progress tracking for {args.project}/{args.sample}", paths.runner_log)

    while True:
        log_content = _read_text(paths.seqclust_log)
        _detect_stages(progress, stages, log_content)
        _estimate_progress_and_eta(progress, stages)
        _update_resources(progress, args.sample)

        paths.progress_json.write_text(json.dumps(progress, indent=2), encoding="utf-8")
        paths.progress_txt.write_text(_render_txt(progress), encoding="utf-8")

        if not _seqclust_running(args.sample):
            progress["eta"] = 0
            if "finalization" in progress["progress"]:
                progress["progress"]["finalization"]["status"] = "completed"
                progress["progress"]["finalization"]["progress"] = 100
            paths.progress_json.write_text(json.dumps(progress, indent=2), encoding="utf-8")
            paths.progress_txt.write_text(_render_txt(progress), encoding="utf-8")
            _log(f"seqclust finished for {args.project}/{args.sample}", paths.runner_log)
            return 0

        time.sleep(max(1, int(args.interval)))


if __name__ == "__main__":
    raise SystemExit(main())

