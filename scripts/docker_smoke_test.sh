#!/usr/bin/env bash
set -euo pipefail

pass() { echo "PASS: $*"; }
fail() { echo "FAIL: $*" >&2; exit 1; }

cd /opt/reportr || fail "missing /opt/reportr (wrong WORKDIR or mount?)"

echo "=== RepOrtR Docker smoke test ==="
echo "pwd: $(pwd)"

echo "--- reportr env: python/yaml ---"
conda run -n reportr python -c "import yaml; print('yaml_ok')"
pass "reportr python can import yaml"

echo "--- reportr env: snakemake ---"
conda run -n reportr snakemake --version
pass "snakemake runs in reportr"

echo "--- reportr env: streamlit/dashboard import ---"
conda run -n reportr python -c "import streamlit; import dashboard.app; print('dashboard_import_ok')"
pass "dashboard imports in reportr"

echo "--- repeatexplorer env: seqclust help ---"
conda run -n repeatexplorer /opt/reportr/repex_tarean/seqclust --help >/dev/null
pass "seqclust --help runs in repeatexplorer"

echo "=== OK ==="

