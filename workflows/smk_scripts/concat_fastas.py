#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--in", dest="inputs", nargs="+", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as w:
        for p in args.inputs:
            path = Path(p)
            if not path.exists():
                continue
            txt = path.read_text()
            w.write(txt)
            if txt and not txt.endswith("\n"):
                w.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

