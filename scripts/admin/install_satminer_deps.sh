#!/usr/bin/env bash
set -euo pipefail

# Install satMiner-related dependencies into the existing `reportr` conda env.
# This script is intentionally best-effort: some satMiner deps (DeconSeq, UCSC faSplit,
# full RepeatMasker stack) may not be available via conda on all platforms.

env_name="${1:-reportr}"

if ! command -v mamba >/dev/null 2>&1; then
  echo "ERROR: mamba not found in PATH."
  echo "Install it first (example): conda install -n base -c conda-forge mamba"
  exit 1
fi

echo "Installing satMiner dependencies into conda env: ${env_name}"

# Prefer bioconda/conda-forge for bio tools.
channels=(-c conda-forge -c bioconda)

# Tools commonly available in conda and used by satMiner-style workflows.
pkgs=(
  # satDNA analysis after RepeatExplorer2
  repeatmasker
)

optional_pkgs=()

set -x
mamba install -n "$env_name" -y "${channels[@]}" "${pkgs[@]}" || true
set +x

echo "Ensuring RepeatMasker util script is available: calcDivergenceFromAlign.pl"
conda_prefix="$(conda run -n "$env_name" python3 -c 'import sys; print(sys.prefix)')"
if [[ -z "${conda_prefix}" ]]; then
  echo "ERROR: Could not determine CONDA_PREFIX for env ${env_name}" >&2
  exit 1
fi

rm_util_dir="${conda_prefix}/share/RepeatMasker/util"
mkdir -p "${rm_util_dir}"
target_util="${rm_util_dir}/calcDivergenceFromAlign.pl"
target_bin="${conda_prefix}/bin/calcDivergenceFromAlign.pl"

need_fetch=0
if [[ ! -f "${target_util}" ]]; then
  need_fetch=1
elif ! head -n 1 "${target_util}" | grep -q '^#!/usr/bin/perl'; then
  # If something replaced it (e.g. an old wrapper), refetch.
  need_fetch=1
fi

if [[ "${need_fetch}" -eq 0 ]]; then
  echo "OK: ${target_util}"
else
  tmpdir="$(mktemp -d)"
  cleanup() { rm -rf "${tmpdir}" || true; }
  trap cleanup EXIT

  url="https://raw.githubusercontent.com/Dfam-consortium/RepeatMasker/master/util/calcDivergenceFromAlign.pl"
  echo "Fetching ${url}"
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "${url}" -o "${tmpdir}/calcDivergenceFromAlign.pl"
  elif command -v wget >/dev/null 2>&1; then
    wget -qO "${tmpdir}/calcDivergenceFromAlign.pl" "${url}"
  else
    echo "ERROR: need curl or wget to fetch RepeatMasker util script" >&2
    exit 1
  fi

  install -m 0755 "${tmpdir}/calcDivergenceFromAlign.pl" "${target_util}"
  echo "Installed: ${target_util}"
fi

# If an old symlink exists, remove it so we don't overwrite the util script.
if [[ -L "${target_bin}" ]]; then
  rm -f "${target_bin}"
fi

cat > "${target_bin}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
rm_root="${conda_prefix}/share/RepeatMasker"
export PERL5LIB="\${rm_root}:\${PERL5LIB:-}"
exec perl "${target_util}" "\$@"
EOF
chmod 0755 "${target_bin}"

echo "Verifying tools inside env ${env_name}"
conda run -n "$env_name" bash -c 'command -v RepeatMasker && command -v calcDivergenceFromAlign.pl && calcDivergenceFromAlign.pl -version >/dev/null'

echo "Done."
echo "NOTE: DeconSeq and UCSC faSplit may require manual installation."

