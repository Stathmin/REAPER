#!/usr/bin/env python3
"""
NCBI freshness checker used by the Snakemake NCBI gathering rules.

Given an NCBI metadata CSV file and a max age threshold (in days), decide
whether the cached data is still fresh enough for post-TAREAN BLAST runs.
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Check NCBI cached data freshness")
    parser.add_argument("--metadata", required=True, help="Path to ncbi_metadata_*.csv")
    parser.add_argument("--max-age-days", required=True, type=int, help="Max age in days")
    parser.add_argument("--output", required=True, help="Write freshness report to this path")
    args = parser.parse_args()

    metadata_path = Path(args.metadata)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    now = datetime.now(timezone.utc)
    if not metadata_path.exists():
        payload = {
            "metadata": str(metadata_path),
            "exists": False,
            "is_fresh": False,
            "max_age_days": args.max_age_days,
            "now_utc": now.isoformat(),
            "reason": "metadata file missing",
        }
        output_path.write_text(json.dumps(payload, indent=2) + "\n")
        return

    # Use filesystem mtime as a simple, robust freshness signal.
    mtime = datetime.fromtimestamp(metadata_path.stat().st_mtime, tz=timezone.utc)
    age_days = (now - mtime).total_seconds() / 86400.0
    is_fresh = age_days <= float(args.max_age_days)

    payload = {
        "metadata": str(metadata_path),
        "exists": True,
        "metadata_mtime_utc": mtime.isoformat(),
        "age_days": age_days,
        "max_age_days": args.max_age_days,
        "is_fresh": is_fresh,
        "now_utc": now.isoformat(),
        "reason": "freshness based on metadata file mtime",
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()

