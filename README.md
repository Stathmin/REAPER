# REAPER: Repeatome Extended Analysis Pipeline — Execution and Reporting. Project-centric repeatome workflows
REAPER is a **project-centric layer** for running **RepeatExplorer2/TAREAN** repeatome analyses in a reproducible, incremental way:
- **Project manager** (`project_manager.py`) registers projects/samples/comparatives and enforces a stable on-disk layout under `projects/`.
- **Modular Snakemake workflow** (`Snakefile_modular` + `workflows/*.smk`) runs QC → deterministic prep → `seqclust` (RepeatExplorer2/TAREAN) → post‑TAREAN tables → optional reports/graphs.
## Start here
- **Workflow structure**: `workflows/README.md`
- **Snakemake profiles**: `profiles/README.md`
## Quickstart (host / conda)
### Install
```bash
python3 install_reportr.py
```
This ensures three conda envs and a built `repex_tarean/` checkout:
- `reportr`: Snakemake + pipeline Python/tooling
- `repeatexplorer`: RepeatExplorer2/TAREAN runtime for `repex_tarean/seqclust`
- `reportr_graph`: R stack for interactive graph reporting (igraph/visNetwork)
### Create a project + add a sample
```bash
conda run -n reportr python3 project_manager.py create-project \
  --project-id my_project \
  --taxonomy Triticeae \
  --taxonomy-id 147389 \
  --description "My repeatome project" \
  --ncbi-repeats data/ncbi_repeats_triticeae.fa \
  --total-reads-per-assembly 50000
conda run -n reportr python3 project_manager.py add-sample \
  --project-id my_project \
  --sample-id S1 \
  --taxonomy Triticeae \
  --r1-path /path/to/S1_R1.fastq.gz \
  --r2-path /path/to/S1_R2.fastq.gz \
  --genome-size 1.0
```
### Dry-run one target (recommended)
```bash
CORES=8
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  -n --use-conda --rerun-incomplete --cores "$CORES" -p \
  projects/my_project/samples/S1/tarean/tarean.done
```
### Run
```bash
CORES=8
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --cores "$CORES" -p \
  projects/my_project/samples/S1/tarean/tarean.done
```
## Documentation index (root)
- `INSTALLATION.md`: install + verify environments and `repex_tarean/`
- `PROJECT_MANAGER.md`: `project_manager.py` + `projects/global_config.yaml` explained
- `RUN_GUIDE.md`: standard Snakemake flag sets, profiles, and end‑to‑end run lifecycle

## Study and example dataset

REAPER is described in:

> Ulyanov, Daniil S., Anna I. Yurkina, Viktoria S. Voronezhskaya, Alana A. Ulyanova, Pavel Yu. Kroupin, and Mikhail G. Divashuk. "REAPER: a project-centric workflow layer for comparative repeatome analysis." *Frontiers in Plant Science* (Plant Bioinformatics section). Manuscript ID 1855960. *[Under review; DOI/public article link to be added once published.]*

The worked *Triticeae* example referenced in the manuscript (`triticeae_F21FTSEUHT1241`; samples KA1–KA5) is shared at:
**https://izba-memes.ru/share/y6TpVoNv**

See that dataset's own README for sample metadata, genome-size sources, and the comparative analyses it was used for; `projects/global_config.yaml` in this repository holds the matching `triticeae_F21FTSEUHT1241` project entry.
