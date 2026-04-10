#!/usr/bin/env python3
"""
Python3 adapter for satMiner sat_subfam2fam.py:
Replace subfamily names in a RepeatMasker .align file according to a pattern table.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--align", required=True, help="RepeatMasker .align file")
    ap.add_argument("--pattern", required=True, help="Tab-separated pattern file: subfam<TAB>fam")
    ap.add_argument("--out", required=True, help="Output .fam align file")
    args = ap.parse_args()

    align = Path(args.align)
    pat = Path(args.pattern)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    patterns = {}
    for line in pat.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        patterns[parts[0]] = parts[1]

    with align.open() as r, out.open("w") as w:
        for line in r:
            for k, v in patterns.items():
                line = line.replace(k, v)
            w.write(line)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

