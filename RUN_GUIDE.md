# Running the workflow (standard flag sets)

The canonical commands live in `docs/TRITICEAE_RUN_COMMANDS.md`. This file provides **copy/paste wrappers**: consistent “flag sets” and an end‑to‑end lifecycle from project init → reports/graphs.

## Always / never

- **Always** run from the repo root.
- **Always** use `Snakefile_modular` and `--configfile projects/global_config.yaml`.
- **Always** include `--use-conda`.
- **Avoid** `--forceall` unless you really want to recompute large parts of the DAG.

## Standard base command

```bash
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml
```

## Profiles (recommended)

Profiles capture stable defaults for long runs (keep-going, rerun-incomplete, resource slots, etc.).

- Local: `--profile profiles/local`
- HPC: `--profile profiles/hpc`
- Quieter output: `--profile profiles/quiet`

Example:

```bash
conda run -n reportr snakemake --profile profiles/local \
  -s Snakefile_modular --configfile projects/global_config.yaml --use-conda \
  -p <targets...>
```

## Flag sets you’ll use often

### 1) Dry-run (validate DAG before spending hours)

```bash
CORES=20
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  -n --use-conda --rerun-incomplete --cores "$CORES" -p \
  <targets...>
```

### 2) Foreground run (interactive)

```bash
CORES=20
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --cores "$CORES" \
  --keep-going --show-failed-logs \
  -p <targets...>
```

### 3) Long run over SSH (nohup + line-buffered output)

```bash
CORES=28
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"
nohup stdbuf -oL -eL conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --rerun-triggers mtime --cores "$CORES" \
  --scheduler greedy \
  --keep-going --show-failed-logs \
  -p <targets...> \
  > "nohup_reportr_${RUN_ID}.out" 2>&1 &
```

### 4) Resource “slot” caps (avoid I/O/RAM contention)

The workflow uses Snakemake resources as **concurrency caps** (see `profiles/*/config*.yaml` and `docs/TRITICEAE_RUN_COMMANDS.md`):

- `seqclust_slots`: keep TAREAN/seqclust serial
- `bbmap_slots`: keep bbmap/bbduk iterative filtering serial
- `repeatmasker_slots`: keep RepeatMasker serial (if enabled)
- `bbduk_slots`, `fastqc_slots`: allow a small amount of overlap for light steps
- `prep_slots`: serialize comparative prep when needed

Example (conservative for slow disks/NFS):

```bash
CORES=28
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --cores "$CORES" \
  --keep-going --show-failed-logs \
  --resources seqclust_slots=1 bbmap_slots=1 repeatmasker_slots=1 prep_slots=1 fastqc_slots=1 bbduk_slots=1 \
  -p <targets...>
```

## End-to-end lifecycle (project → samples → reports/graphs)

### A) Initiate a project

```bash
conda run -n reportr python3 project_manager.py create-project \
  --project-id <project_id> \
  --taxonomy <taxonomy_name> \
  --taxonomy-id <ncbi_taxid> \
  --description "<free text>" \
  --ncbi-repeats <path/to/repeats.fasta> \
  --total-reads-per-assembly 50000
```

### B) Register samples (creates directories + symlinks + config entries)

```bash
conda run -n reportr python3 project_manager.py add-sample \
  --project-id <project_id> --sample-id <sample_id> --taxonomy <taxonomy_name> \
  --r1-path /abs/path/to/R1.fastq.gz \
  --r2-path /abs/path/to/R2.fastq.gz \
  --genome-size 1.0
```

### C) Register comparative analyses (explicitly)

```bash
conda run -n reportr python3 project_manager.py add-comparative \
  --project-id <project_id> \
  --analysis-id <analysis_id> \
  --samples <sample_id_1> <sample_id_2> \
  --analysis-description "comparison description"
```

### D) Validate quickly (config + filesystem coherence)

```bash
conda run -n reportr python3 project_manager.py validate --project-id <project_id>

# Snakemake-side config validation target (fast fail)
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --cores 4 -p logs/config_validation.txt
```

### E) Run sample TAREAN targets

Targets are **file paths** in the `projects/` tree.

```bash
CORES=20
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --cores "$CORES" -p \
  projects/<project_id>/samples/<sample_id>/tarean/tarean.done
```

If iterative assembly is enabled (`iterative_assembly.depth: 3`), target:

```bash
projects/<project_id>/samples/<sample_id>/tarean3/tarean.done
```

### F) Run post‑TAREAN “rich report” (per sample)

```bash
CORES=20
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --cores "$CORES" -p \
  projects/<project_id>/samples/<sample_id>/post_tarean/satminer/report/satminer_rich_report.done
```

### G) Run comparative assemblies

```bash
CORES=20
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --cores "$CORES" -p \
  projects/<project_id>/comparative/<analysis_id>/COMPARATIVE_TAREAN_COMPLETE
```

For iterative comparative (depth 3 example):

```bash
projects/<project_id>/comparative/<analysis_id>/tarean3/COMPARATIVE_TAREAN_COMPLETE
```

### H) Project-level graph/reporting (debug without re-running assemblies)

If you only want the graph report and do **not** want to re-run seqclust/TAREAN, keep the rule set narrow (this is in the canonical runbook):

```bash
conda run -n reportr snakemake -s Snakefile_modular --configfile projects/global_config.yaml \
  --use-conda --rerun-incomplete --keep-going --show-failed-logs \
  --cores 16 \
  --allowed-rules project_repeat_interaction_graph_report project_localdb_from_rich_tables project_ncbi_repeats_blastdb project_oligo_blastdb \
  -p \
  projects/<project_id>/reports/repeat_interaction_graph/graph_report.done
```

## Rerun patterns (safe defaults)

- **After interruption**: add `--rerun-incomplete`.
- **After code/config change**:
  - for stability-sensitive long runs, many examples use `--rerun-triggers mtime`
  - for provenance-driven recompute, prefer `--rerun-triggers input params code software-env`

See `docs/TRITICEAE_RUN_COMMANDS.md` for the repo’s current recommended choices.

