# Project manager + global config

The project lifecycle is driven by two things:

- **CLI**: `project_manager.py` (creates/updates project structure and writes config)
- **Config**: `projects/global_config.yaml` (source of truth for global + per-project + per-sample + comparative settings)

Users should **not** create `projects/...` directories manually.

## Repository model (what exists where)

At a high level:

- `projects/global_config.yaml`: global defaults + a map of projects
- `projects/<project_id>/`: on-disk project area created/managed by the CLI
  - `samples/<sample_id>/`: per-sample layout (symlinked raw reads, filtered reads, TAREAN outputs, post‑TAREAN outputs)
  - `comparative/<analysis_id>/`: comparative layouts

## `projects/global_config.yaml` structure

Top-level keys (typical):

- `global`: execution defaults and shared parameters
- `projects`: mapping of `project_id -> project config`

Within a project:

- `taxonomy`, `taxonomy_id` / `taxonomy_ids`
- `description`
- `ncbi_repeats`: path to a FASTA used as the project repeat database
- `samples`: mapping of sample IDs to sample metadata (`r1_path`, `r2_path`, `genome_size`, `prefix`, …)
- `comparative_analyses`: mapping of `analysis_id -> {samples: [...], description: ...}`

## `project_manager.py` (what it does)

`project_manager.py` enforces:

- stable directory layout under `projects/<project_id>/`
- symlinks to raw FASTQs under `projects/<project_id>/samples/<sample_id>/raw_reads/`
- unique sample prefixes (used by RepeatExplorer/seqclust labelling)
- config updates in `projects/global_config.yaml`

Key commands (see the canonical examples in `docs/TRITICEAE_RUN_COMMANDS.md`):

### Create a project

```bash
conda run -n reportr python3 project_manager.py create-project \
  --project-id triticeae_F21FTSEUHT1241 \
  --taxonomy Triticeae \
  --taxonomy-id 147389 \
  --description "Triticeae samples from F21FTSEUHT1241_PLAsynvR" \
  --ncbi-repeats data/ncbi_repeats_triticeae.fa \
  --total-reads-per-assembly 50000 \
  --tarean-options ILLUMINA_SENSITIVE_BLASTPLUS \
  --tarean-mincl 0.0001 \
  --tarean-assembly-min 2 \
  --tarean-domain-search DIAMOND \
  --tarean-cleanup \
  --tarean-automatic-filtering
```

### Add a sample

```bash
conda run -n reportr python3 project_manager.py add-sample \
  --project-id triticeae_F21FTSEUHT1241 \
  --sample-id KA1 \
  --taxonomy Triticeae \
  --r1-path /storage/.../KA1_R1.fq.gz \
  --r2-path /storage/.../KA1_R2.fq.gz \
  --genome-size 1.5
```

This creates (among other paths):

- `projects/<project>/samples/<sample>/raw_reads/R1.fq[.gz]` (symlink)
- `projects/<project>/samples/<sample>/raw_reads/R2.fq[.gz]` (symlink)
- `projects/<project>/samples/<sample>/tarean/`
- `projects/<project>/samples/<sample>/post_tarean/`

### Add a comparative

```bash
conda run -n reportr python3 project_manager.py add-comparative \
  --project-id triticeae_F21FTSEUHT1241 \
  --analysis-id KA1_KA2 \
  --samples KA1 KA2 \
  --analysis-description "KA1 vs KA2 comparative"
```

### Validate config + filesystem consistency

```bash
conda run -n reportr python3 project_manager.py validate --project-id triticeae_F21FTSEUHT1241
```

## What you typically edit by hand

Even when using the CLI, you will often adjust `projects/global_config.yaml` directly.

High-impact knobs:

- `global.total_reads_per_assembly` (default read budget)
- `global.tarean_params.*` and project overrides at `projects.<id>.tarean_params.*`
- `global.post_tarean_params.graph_report.*` (graph report filters + limits)
- `projects.<id>.iterative_assembly.*` (enable/depth)
- `projects.<id>.samples.<sid>.(r1_path,r2_path,genome_size,prefix)`

For details on override resolution and comparative strictness, see `docs/TRITICEAE_RUN_COMMANDS.md` (it documents the effective hierarchy used by the workflow).

