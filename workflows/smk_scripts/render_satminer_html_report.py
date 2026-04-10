#!/usr/bin/env python3

from __future__ import annotations

import argparse
import html
from pathlib import Path

import pandas as pd


def _read_tsv(path: Path, n: int = 30) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()
    return df.head(n)


def _df_to_html(df: pd.DataFrame) -> str:
    if df.empty:
        return "<p><em>(no rows)</em></p>"
    return df.to_html(index=False, escape=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--blast-top-tsv", required=True)
    ap.add_argument("--hits-tsv", required=True)
    ap.add_argument("--div-tsv", required=True)
    ap.add_argument("--plots-dir", required=True)
    args = ap.parse_args()

    sample = args.sample
    out_html = Path(args.out_html)
    out_html.parent.mkdir(parents=True, exist_ok=True)

    blast_top = Path(args.blast_top_tsv)
    hits = Path(args.hits_tsv)
    div = Path(args.div_tsv)
    plots_dir = Path(args.plots_dir)

    top_hits_png = plots_dir / "top_hits.png"
    pdiv_png = plots_dir / "pdiv_by_iter.png"

    df_blast = _read_tsv(blast_top, 30)
    df_hits = _read_tsv(hits, 30)
    df_div = _read_tsv(div, 30)

    def img_tag(p: Path, alt: str) -> str:
        if not p.exists():
            return "<p><em>(plot not generated)</em></p>"
        # Link to file; keep HTML lightweight (no inline base64).
        rel = html.escape(str(p))
        return f'<a href="{rel}"><img src="{rel}" alt="{html.escape(alt)}" style="max-width: 100%; height: auto;"/></a>'

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>RepOrtR satMiner report: {html.escape(sample)}</title>
  <style>
    body {{ font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial; margin: 24px; }}
    h1,h2 {{ margin: 0.8em 0 0.3em; }}
    .grid {{ display: grid; grid-template-columns: 1fr; gap: 18px; }}
    .card {{ border: 1px solid #ddd; border-radius: 10px; padding: 14px 16px; }}
    table {{ border-collapse: collapse; width: 100%; font-size: 12px; }}
    th, td {{ border: 1px solid #eee; padding: 6px 8px; text-align: left; }}
    th {{ background: #fafafa; }}
    code {{ background: #f6f6f6; padding: 2px 6px; border-radius: 6px; }}
  </style>
</head>
<body>
  <h1>satMiner comparisons: {html.escape(sample)}</h1>
  <p>Inputs are produced under <code>post_tarean/satminer/</code>.</p>

  <div class="grid">
    <div class="card">
      <h2>BLAST: tagged consensus vs project NCBI repeats DB (top hit per query)</h2>
      {_df_to_html(df_blast)}
      <p>Full table: <code>{html.escape(str(blast_top))}</code></p>
    </div>

    <div class="card">
      <h2>Abundance proxy (RepeatMasker hit counts per consensus tag)</h2>
      {img_tag(top_hits_png, "Top hits")}
      {_df_to_html(df_hits)}
      <p>Table: <code>{html.escape(str(hits))}</code></p>
    </div>

    <div class="card">
      <h2>Divergence proxy (RepeatMasker perc div.)</h2>
      {img_tag(pdiv_png, "pdiv by iter")}
      {_df_to_html(df_div)}
      <p>Table: <code>{html.escape(str(div))}</code></p>
    </div>
  </div>
</body>
</html>
"""

    out_html.write_text(html_doc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

