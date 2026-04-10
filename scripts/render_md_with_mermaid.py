#!/usr/bin/env python3
from __future__ import annotations

import argparse
import base64
import mimetypes
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class RenderOpts:
    in_md: Path
    out_html: Path
    base_dir: Path
    embed_images: bool


_MERMAID_FENCE_RE = re.compile(r"(^```mermaid\s*$)([\s\S]*?)(^```\s*$)", re.MULTILINE)


def _read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8")


def _sniff_mime(p: Path) -> str:
    mt, _ = mimetypes.guess_type(str(p))
    if mt:
        return mt
    # Fallbacks for common figure types
    ext = p.suffix.lower()
    if ext == ".svg":
        return "image/svg+xml"
    if ext == ".pdf":
        return "application/pdf"
    if ext in {".png", ".jpg", ".jpeg", ".gif", ".webp"}:
        return f"image/{ext.lstrip('.')}"
    return "application/octet-stream"


def _data_uri_for_file(p: Path) -> str:
    raw = p.read_bytes()
    mime = _sniff_mime(p)
    b64 = base64.b64encode(raw).decode("ascii")
    return f"data:{mime};base64,{b64}"


def _rewrite_mermaid_fences(md: str) -> str:
    def _repl(m: re.Match[str]) -> str:
        code = m.group(2).strip("\n")
        # Many of our Mermaid sources use literal "\n" inside quoted labels.
        # Mermaid prefers <br/> for line breaks in flowchart node labels.
        code = code.replace("\\n", "<br/>")
        # Mermaid reads from <pre class="mermaid">...</pre>
        return f'<pre class="mermaid">\n{code}\n</pre>'

    return _MERMAID_FENCE_RE.sub(_repl, md)


def _iter_image_srcs(html: str) -> Iterable[str]:
    # Minimal HTML-level rewrite (produced by markdown engines): <img ... src="...">
    for m in re.finditer(r'<img[^>]+src="([^"]+)"', html):
        yield m.group(1)


def _embed_images_in_html(html: str, base_dir: Path) -> str:
    # Only embed relative paths; keep http(s), data: unchanged.
    def _embed_src(src: str) -> str:
        if src.startswith(("http://", "https://", "data:")):
            return src
        # Allow paths like "figures/foo.svg" relative to base_dir
        p = (base_dir / src).resolve()
        if not p.exists() or not p.is_file():
            return src
        return _data_uri_for_file(p)

    def _replace(match: re.Match[str]) -> str:
        full = match.group(0)
        quote = match.group("q")
        src = match.group("src")
        new_src = _embed_src(src)
        if new_src == src:
            return full
        return full.replace(f"src={quote}{src}{quote}", f"src={quote}{new_src}{quote}")

    # Handles src="..." and src='...'
    img_src_re = re.compile(r'<img[^>]+src=(?P<q>["\'])(?P<src>[^"\']+)(?P=q)', re.IGNORECASE)
    return img_src_re.sub(_replace, html)


def _markdown_to_html(md: str) -> str:
    # Import locally so the script can print a clean error if deps missing.
    try:
        import markdown  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "Missing renderer dependencies. Install with:\n"
            "  python -m pip install markdown\n"
        ) from e

    # Use a small set of built-in extensions; keep behavior predictable.
    return markdown.markdown(
        md,
        extensions=[
            "tables",
            "fenced_code",
        ],
        output_format="html5",
    )


def _wrap_html_document(body_html: str, title: str) -> str:
    # Mermaid is rendered client-side.
    mermaid_js = "https://cdn.jsdelivr.net/npm/mermaid@10.9.1/dist/mermaid.min.js"
    css = """
    :root { color-scheme: light; }
    body { margin: 0; padding: 40px 24px; font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, \"Apple Color Emoji\", \"Segoe UI Emoji\"; line-height: 1.45; color: #111; background: #fff; }
    main { max-width: 900px; margin: 0 auto; }
    h1,h2,h3 { line-height: 1.2; }
    h1 { font-size: 28px; margin: 0 0 16px; }
    h2 { font-size: 22px; margin-top: 28px; }
    h3 { font-size: 18px; margin-top: 22px; }
    code { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, \"Liberation Mono\", \"Courier New\", monospace; font-size: 0.95em; }
    pre { background: #f6f8fa; padding: 12px 14px; overflow: auto; border-radius: 8px; }
    table { border-collapse: collapse; width: 100%; margin: 12px 0; }
    th, td { border: 1px solid #ddd; padding: 6px 8px; vertical-align: top; }
    th { background: #f3f4f6; text-align: left; }
    img { max-width: 100%; height: auto; }
    blockquote { margin: 12px 0; padding-left: 14px; border-left: 3px solid #ddd; color: #444; }
    a { color: #0b57d0; text-decoration: none; }
    a:hover { text-decoration: underline; }
    """
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{title}</title>
  <style>{css}</style>
</head>
<body>
  <main>
{body_html}
  </main>
  <script src="{mermaid_js}"></script>
  <script>
    // Render Mermaid diagrams produced from fenced blocks.
    mermaid.initialize({{ startOnLoad: true, theme: "default" }});
  </script>
</body>
</html>
"""


def render(opts: RenderOpts) -> None:
    md_raw = _read_text(opts.in_md)
    # Inline the lifecycle schematic if referenced, so Mermaid content is present
    # in the rendered HTML even when the manuscript links to a separate doc.
    def _inline_appendix(rel_path: str, heading: str) -> None:
        nonlocal md_raw
        if rel_path not in md_raw:
            return
        p = (opts.base_dir / rel_path).resolve()
        if not p.exists():
            return
        md_raw = md_raw + f"\n\n## {heading}\n\n" + _read_text(p)

    _inline_appendix("project_management_lifecycle.md", "Appendix: Project management lifecycle schematic")
    _inline_appendix("project_management_entities_outputs.md", "Appendix: Project management entities and outputs")
    md_pre = _rewrite_mermaid_fences(md_raw)
    body_html = _markdown_to_html(md_pre)
    if opts.embed_images:
        body_html = _embed_images_in_html(body_html, opts.base_dir)

    title = opts.in_md.name
    html = _wrap_html_document(body_html, title=title)
    opts.out_html.parent.mkdir(parents=True, exist_ok=True)
    opts.out_html.write_text(html, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="Render Markdown to self-contained HTML with Mermaid + embedded images.")
    ap.add_argument("--in", dest="in_md", required=True, help="Input Markdown file")
    ap.add_argument("--out", dest="out_html", required=True, help="Output HTML file")
    ap.add_argument(
        "--base-dir",
        default=".",
        help="Base directory to resolve relative image paths (default: .)",
    )
    ap.add_argument(
        "--embed-images",
        action="store_true",
        default=True,
        help="Embed relative <img src> as data URIs (default: enabled)",
    )
    args = ap.parse_args()

    opts = RenderOpts(
        in_md=Path(args.in_md),
        out_html=Path(args.out_html),
        base_dir=Path(args.base_dir),
        embed_images=bool(args.embed_images),
    )
    render(opts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

