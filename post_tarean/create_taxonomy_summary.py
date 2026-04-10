#!/usr/bin/env python3
"""
Minimal taxonomy summary creator for RepOrtR NCBI gathering.

This matches the CLI used by the Snakemake rule:
  python3 post_tarean/create_taxonomy_summary.py --input-files ... --output ...

It currently:
- Reads one or more metadata CSV files (as produced by ncbi_data_gatherer).
- Concatenates them and writes a compact summary CSV with per-accession metadata.
"""

import argparse
from pathlib import Path
from typing import List

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create NCBI taxonomy summary from metadata CSV files.")
    parser.add_argument(
        "--input-files",
        nargs="+",
        required=True,
        help="Input NCBI metadata CSV files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output taxonomy summary CSV.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_paths: List[Path] = [Path(p) for p in args.input_files]
    out_path = Path(args.output)

    out_path.parent.mkdir(parents=True, exist_ok=True)

    frames = []
    for p in input_paths:
        if p.exists():
            frames.append(pd.read_csv(p))

    if not frames:
        # Nothing to summarize; write an empty file with a header.
        pd.DataFrame(
            columns=[
                "accession",
                "organism",
                "taxonomy",
                "length",
                "gc_content",
                "source_file",
            ]
        ).to_csv(out_path, index=False)
        return

    df = pd.concat(frames, ignore_index=True)
    df["source_file"] = df.get("source_file", None)

    # Keep key columns; tolerate missing ones.
    cols = [c for c in ["accession", "organism", "taxonomy", "length", "gc_content", "source_file"] if c in df.columns]
    df[cols].to_csv(out_path, index=False)


if __name__ == "__main__":
    main()

