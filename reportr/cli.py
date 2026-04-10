#!/usr/bin/env python3
"""Top-level CLI for RepOrtR."""

from __future__ import annotations

import argparse
import subprocess
import sys


def _run(cmd: list[str]) -> int:
    proc = subprocess.run(cmd, check=False)
    return proc.returncode


def main() -> None:
    parser = argparse.ArgumentParser(description="RepOrtR command line interface")
    sub = parser.add_subparsers(dest="command", required=True)

    lock = sub.add_parser("lock", help="Generate conda-lock files")
    lock.add_argument(
        "--env",
        choices=["reportr", "repeatexplorer", "all"],
        default="all",
        help="Environment lock target",
    )
    lock.add_argument(
        "--platform",
        default="linux-64",
        help="Conda platform (default: linux-64)",
    )

    monitor = sub.add_parser("monitor", help="Run unified monitoring CLI")
    monitor.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="Arguments forwarded to `reportr-monitor`",
    )

    args = parser.parse_args()

    if args.command == "monitor":
        code = _run(["python3", "-m", "monitoring.cli", *args.args])
        sys.exit(code)

    if args.command == "lock":
        targets = ["reportr", "repeatexplorer"] if args.env == "all" else [args.env]
        for env_name in targets:
            env_file = f"envs/{env_name}.yaml"
            lock_path = f"envs/locks/{env_name}-{args.platform}.lock.yml"
            cmd = [
                "conda-lock",
                "lock",
                "-f",
                env_file,
                "-p",
                args.platform,
                "--lockfile",
                lock_path,
            ]
            code = _run(cmd)
            if code != 0:
                sys.exit(code)
        return


if __name__ == "__main__":
    main()
