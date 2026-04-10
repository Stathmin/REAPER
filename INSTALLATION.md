# Installation (host / conda)

This repo uses **Snakemake + rule-level conda environments** and expects host execution (no container required).

## What gets installed

Running `install_reportr.py` ensures:

- **Conda envs** (from repo YAMLs in `envs/`):
  - `reportr` (`envs/reportr.yaml`): Snakemake + pipeline tooling
  - `repeatexplorer` (`envs/repeatexplorer.yaml`): RepeatExplorer2/TAREAN runtime
  - `reportr_graph` (`envs/reportr_graph.yaml`): R/igraph/visNetwork for graph reports
- **RepeatExplorer wrapper checkout**:
  - `./repex_tarean` cloned into repo root (if missing)
  - built with `make` inside the `repeatexplorer` env

Snakemake must be run with `--use-conda` so rules pick the correct env.

## Prerequisites

- **Conda** (Miniconda or Anaconda)
- **Git**
- A basic build toolchain for `make` (to compile `repex_tarean` utilities)

## Install / update (recommended)

From repo root:

```bash
python3 install_reportr.py
```

The installer will:

1. Find `conda`
2. Create or update envs in-place
3. Clone `repex_tarean` if needed
4. Build and verify `repex_tarean/seqclust`
5. Verify key tools in each env

## Verify manually (quick)

```bash
# Snakemake is in the reportr env
conda run -n reportr snakemake --version

# seqclust is built in the repo, run under repeatexplorer
conda run -n repeatexplorer repex_tarean/seqclust --help

# graph reporting stack (R)
conda run -n reportr_graph R -q -e 'library(igraph); library(visNetwork); quit(save="no")'
```

## Common issues

### Conda not found

- Ensure `conda` is on `PATH` (or run `conda init`, restart the shell).

### Snakemake runs but rules fail with missing tools

- You probably forgot `--use-conda`.

### `make` fails in `repex_tarean/`

- Install a build toolchain on the host (compiler + `make`).
- Then rerun `python3 install_reportr.py`.

### Stale Snakemake lock

If you interrupted a run:

```bash
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml --unlock
```

