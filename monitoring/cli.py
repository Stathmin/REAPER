#!/usr/bin/env python3
"""Unified monitoring CLI for RepOrtR."""

from __future__ import annotations

import argparse
import sys


def main() -> None:
    parser = argparse.ArgumentParser(description="RepOrtR monitoring entrypoint")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("progress", help="Run JSON-based progress dashboard")
    sub.add_parser("pipeline", help="Run legacy pipeline monitor")
    sub.add_parser("holy", help="Run legacy holy monitor")

    stability = sub.add_parser("stability", help="Run stability monitor wrapper")
    stability.add_argument("project_id", help="Project ID to monitor")
    stability.add_argument("--target", help="Specific Snakemake target")
    stability.add_argument("--cores", type=int, default=16, help="Number of cores")
    stability.add_argument("--memory", default="80G", help="Memory limit")
    stability.add_argument(
        "--max-memory",
        type=int,
        default=80,
        help="Max memory usage percent before throttling",
    )
    stability.add_argument(
        "--checkpoint-interval",
        type=int,
        default=3600,
        help="Checkpoint interval seconds",
    )

    args = parser.parse_args()

    if args.command == "progress":
        from monitoring.progress_monitor import main as progress_main

        progress_main()
        return

    if args.command == "pipeline":
        from monitoring.monitor_pipeline import main as pipeline_main

        pipeline_main()
        return

    if args.command == "holy":
        from monitoring.holy_monitor import HolyMonitor

        HolyMonitor().run()
        return

    if args.command == "stability":
        from monitoring.stability_monitor import StabilityMonitor

        monitor = StabilityMonitor(
            args.project_id,
            max_memory_percent=args.max_memory,
            checkpoint_interval=args.checkpoint_interval,
        )
        ok = monitor.run_snakemake_with_stability(
            target=args.target,
            cores=args.cores,
            memory=args.memory,
        )
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
