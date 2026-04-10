#!/usr/bin/env bash
set -euo pipefail

pass() { echo "PASS: $*"; }
fail() { echo "FAIL: $*" >&2; exit 1; }

cd /opt/reportr || fail "missing /opt/reportr (wrong WORKDIR or mount?)"

echo "=== RepOrtR Docker full test ==="
echo "pwd: $(pwd)"

echo "--- reportr env: pytest core suite (post_tarean) ---"
conda run -n reportr pytest -q tests/test_post_tarean_pipeline.py
pass "pytest post_tarean suite passed"

echo "--- reportr env: modular Snakemake dry-run ---"
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml --dry-run --cores 1
pass "modular Snakemake dry-run completed"

echo "=== FULL TEST OK ==="

