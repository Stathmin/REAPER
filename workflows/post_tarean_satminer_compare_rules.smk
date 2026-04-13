# Post-tarean satMiner comparisons and plots
#
# Produces TSVs + plots + HTML report under:
#   projects/{project}/samples/{sample}/post_tarean/satminer/

import os
from pathlib import Path


def _pick_newest_repeats_fasta(project: str) -> str:
    db_dir = Path(f"projects/{project}/blast_db")
    cands = []
    if db_dir.exists():
        cands.extend(sorted(db_dir.glob("ncbi_repeats_*.fasta"), key=lambda p: p.stat().st_mtime, reverse=True))
    if not cands:
        raise FileNotFoundError(f"No ncbi_repeats_*.fasta found under {db_dir}")
    return str(cands[0])

def _sample_attr(project: str, sample: str, key: str) -> str:
    proj = config.get("projects", {}).get(project, {})
    smpl = proj.get("samples", {}).get(sample, {})
    v = smpl.get(key, "")
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return ",".join(str(x) for x in v)
    return str(v)

def _project_oligo_fasta(project: str) -> str:
    proj = config.get("projects", {}).get(project, {}) or {}
    v = proj.get("oligo_fasta", None)
    if not v:
        raise KeyError(f"Missing required config: projects.{project}.oligo_fasta")
    return str(v)


rule project_ncbi_repeats_blastdb:
    """Ensure the project's NCBI repeats BLAST DB exists (makeblastdb from FASTA)."""
    input:
        fasta=lambda w: _pick_newest_repeats_fasta(w.project)
    output:
        ndb="projects/{project}/blast_db/ncbi_repeats.ndb",
        nhr="projects/{project}/blast_db/ncbi_repeats.nhr",
        nin="projects/{project}/blast_db/ncbi_repeats.nin",
        not_="projects/{project}/blast_db/ncbi_repeats.not",
        nsq="projects/{project}/blast_db/ncbi_repeats.nsq",
        ntf="projects/{project}/blast_db/ncbi_repeats.ntf",
        nto="projects/{project}/blast_db/ncbi_repeats.nto",
        njs="projects/{project}/blast_db/ncbi_repeats.njs",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/project_ncbi_repeats_blastdb_{{project}}.log"
    benchmark:
        "benchmarks/project_ncbi_repeats_blastdb_{project}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the project (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/tmp/jobs/project_ncbi_repeats_blastdb"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): PROJECT_NCBI_BLASTDB_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "projects/{wildcards.project}/blast_db"
        python3 workflows/smk_scripts/build_project_blastdb.py \
          --fasta "{input.fasta}" \
          --out-prefix "projects/{wildcards.project}/blast_db/ncbi_repeats"
        test -f "{output.nsq}"
        """

rule project_oligo_blastdb:
    """Ensure the project's oligo BLAST DB exists (makeblastdb from FASTA)."""
    input:
        fasta=lambda w: _project_oligo_fasta(w.project)
    output:
        ndb="projects/{project}/blast_db/oligo.ndb",
        nhr="projects/{project}/blast_db/oligo.nhr",
        nin="projects/{project}/blast_db/oligo.nin",
        not_="projects/{project}/blast_db/oligo.not",
        nsq="projects/{project}/blast_db/oligo.nsq",
        ntf="projects/{project}/blast_db/oligo.ntf",
        nto="projects/{project}/blast_db/oligo.nto",
        njs="projects/{project}/blast_db/oligo.njs",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/project_oligo_blastdb_{{project}}.log"
    benchmark:
        "benchmarks/project_oligo_blastdb_{project}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the project (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/tmp/jobs/project_oligo_blastdb"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): PROJECT_OLIGO_BLASTDB_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "projects/{wildcards.project}/blast_db"
        python3 workflows/smk_scripts/build_project_blastdb.py \
          --fasta "{input.fasta}" \
          --out-prefix "projects/{wildcards.project}/blast_db/oligo" \
          --no-parse-seqids
        test -f "{output.nsq}"
        """


def _project_rich_tables(project: str):
    """Return full rich_table.tsv paths for ALL configured samples + comparatives."""
    proj_cfg = config["projects"][project]
    sample_names = sorted((proj_cfg.get("samples") or {}).keys())
    comparative_names = sorted((proj_cfg.get("comparative_analyses") or {}).keys())

    paths = []
    paths.extend(
        f"projects/{project}/samples/{s}/post_tarean/satminer/report/rich_table.tsv"
        for s in sample_names
    )
    paths.extend(
        f"projects/{project}/comparative/{c}/post_tarean/satminer/report/rich_table.tsv"
        for c in comparative_names
    )
    return paths


rule project_localdb_from_rich_tables:
    """Build canonical local BLAST DB (projects/<project>/blast_db/multifasta.fasta) from sample+comparative rich tables."""
    input:
        rich_tables=lambda w: _project_rich_tables(w.project)
    output:
        fasta="projects/{project}/blast_db/multifasta.fasta",
        nhr="projects/{project}/blast_db/multifasta.fasta.nhr",
        token=touch("projects/{project}/blast_db/localdb.ready"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/project_localdb_from_rich_tables_{{project}}.log"
    benchmark:
        "benchmarks/project_localdb_from_rich_tables_{project}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the project (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/tmp/jobs/project_localdb_from_rich_tables"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): PROJECT_LOCALDB_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "projects/{wildcards.project}/blast_db"
        python3 workflows/smk_scripts/build_localdb_from_rich_tables.py \
          --project "{wildcards.project}" \
          --rich-tsv {input.rich_tables} \
          --out-fasta "{output.fasta}" \
          --out-manifest "projects/{wildcards.project}/blast_db/multifasta.manifest.tsv"

        python3 workflows/smk_scripts/build_project_blastdb.py \
          --fasta "{output.fasta}" \
          --out-prefix "projects/{wildcards.project}/blast_db/multifasta.fasta" \
          --no-parse-seqids

        test -f "{output.nhr}"
        touch "{output.token}"
        """


def _subject_from_rich_top(rich_top: str) -> str:
    p = Path(str(rich_top))
    parts = list(p.parts)
    if "samples" in parts:
        i = parts.index("samples")
        return parts[i + 1]
    if "comparative" in parts:
        i = parts.index("comparative")
        return parts[i + 1]
    return p.stem


def _rich_top_for_subject(project: str, subject: str) -> str:
    # Prefer sample layout; fall back to comparative.
    s = Path(f"projects/{project}/samples/{subject}/post_tarean/satminer/report/rich_table.core.top.tsv")
    if s.exists():
        return str(s)
    c = Path(f"projects/{project}/comparative/{subject}/post_tarean/satminer/report/rich_table.core.top.tsv")
    return str(c)


def _consensus_fasta_for_subject(project: str, subject: str) -> str:
    s = Path(f"projects/{project}/samples/{subject}/post_tarean/satminer/consensus_all_iters.tagged.fasta")
    if s.exists():
        return str(s)
    c = Path(f"projects/{project}/comparative/{subject}/post_tarean/satminer/consensus_all_iters.tagged.fasta")
    return str(c)


rule project_central_blast_subject:
    """Centralized BLAST per ready subject: query latest-iter consensus vs local+ncbi+oligo DBs."""
    input:
        rich_top=lambda w: _rich_top_for_subject(w.project, w.subject),
        consensus=lambda w: _consensus_fasta_for_subject(w.project, w.subject),
        localdb="projects/{project}/blast_db/localdb.ready",
        ncbi_nsq="projects/{project}/blast_db/ncbi_repeats.nsq",
        oligo_nsq="projects/{project}/blast_db/oligo.nsq",
    output:
        query="projects/{project}/blast_db/central_blast/{subject}/queries.latest.fasta",
        local_full="projects/{project}/blast_db/central_blast/{subject}/local.full.tsv",
        local_best="projects/{project}/blast_db/central_blast/{subject}/local.best.tsv",
        ncbi_full="projects/{project}/blast_db/central_blast/{subject}/ncbi.full_x3.tsv",
        ncbi_best="projects/{project}/blast_db/central_blast/{subject}/ncbi.best_x3.tsv",
        token=touch("projects/{project}/blast_db/central_blast/{subject}/central_blast.done"),
    conda:
        "../envs/reportr.yaml"
    threads: 8
    log:
        f"{LOG_DIR}/project_central_blast_subject_{{project}}_{{subject}}.log"
    benchmark:
        "benchmarks/project_central_blast_subject_{project}_{subject}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        tmp_base="projects/{wildcards.project}/tmp/jobs/project_central_blast_subject"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): CENTRAL_BLAST_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "$(dirname "{output.query}")"

        python3 workflows/smk_scripts/extract_latest_iter_queries_from_rich_table.py \
          --rich-top-tsv "{input.rich_top}" \
          --consensus-fasta "{input.consensus}" \
          --out-fasta "{output.query}"

        # LocalDB: explicit prefix (avoid “newest *.nsq” ambiguity in blast_db dir).
        python3 workflows/smk_scripts/blast_tagged_consensus_vs_explicit_db.py \
          --query-fasta "{output.query}" \
          --db-prefix "projects/{wildcards.project}/blast_db/multifasta.fasta" \
          --out-tsv "{output.local_full}" \
          --out-top-tsv "{output.local_best}" \
          --threads "{threads}" \
          --task megablast

        # NCBI repeats + oligo evidence (x1+x3 + task ladder) using satMiner engine.
        python3 workflows/smk_scripts/satminer_blast_x3.py \
          --project "{wildcards.project}" \
          --query-fasta "{output.query}" \
          --blast-db-prefix "projects/{wildcards.project}/blast_db/ncbi_repeats" \
          --oligo-db-prefix "projects/{wildcards.project}/blast_db/oligo" \
          --out-full-tsv "{output.ncbi_full}" \
          --out-best-tsv "{output.ncbi_best}" \
          --threads "{threads}"

        touch "{output.token}"
        """


def _central_blast_subject_tokens(project: str):
    # Discover ready subjects from existing rich tables, then depend on their per-subject tokens.
    base = Path(f"projects/{project}")
    if not base.exists():
        return []
    rich_paths = []
    rich_paths.extend(sorted(base.glob("samples/*/post_tarean/satminer/report/rich_table.core.top.tsv")))
    rich_paths.extend(sorted(base.glob("comparative/*/post_tarean/satminer/report/rich_table.core.top.tsv")))
    subjects = sorted({_subject_from_rich_top(str(p)) for p in rich_paths})
    return [f"projects/{project}/blast_db/central_blast/{s}/central_blast.done" for s in subjects]


rule project_central_blast_done:
    """Project-wide token: centralized BLAST completed for all ready subjects."""
    input:
        localdb="projects/{project}/blast_db/localdb.ready",
        per_subject=lambda w: _central_blast_subject_tokens(w.project),
    output:
        token=touch("projects/{project}/blast_db/central_blast/central_blast.done"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/project_central_blast_done_{{project}}.log"
    benchmark:
        "benchmarks/project_central_blast_done_{project}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): centralized BLAST complete token" >&2
        touch "{output.token}"
        """


rule project_repeat_interaction_graph_report:
    """Project-scoped repeat interaction graph report (LocalDB clusters + NCBI/oligo annotation)."""
    input:
        localdb="projects/{project}/blast_db/localdb.ready",
        ncbi_nsq="projects/{project}/blast_db/ncbi_repeats.nsq",
        oligo_nsq="projects/{project}/blast_db/oligo.nsq",
        rich_tables=lambda w: _project_rich_tables(w.project),
    output:
        queries="projects/{project}/reports/repeat_interaction_graph/queries.latest.fasta",
        index="projects/{project}/reports/repeat_interaction_graph/queries.index.tsv",
        nodes="projects/{project}/reports/repeat_interaction_graph/nodes.tsv",
        hsps="projects/{project}/reports/repeat_interaction_graph/localdb.hsps.tsv",
        edges_full="projects/{project}/reports/repeat_interaction_graph/edges.full.tsv",
        edges="projects/{project}/reports/repeat_interaction_graph/edges.tsv",
        clusters="projects/{project}/reports/repeat_interaction_graph/clusters.tsv",
        cluster_sizes="projects/{project}/reports/repeat_interaction_graph/cluster_sizes.tsv",
        cluster_consensus="projects/{project}/reports/repeat_interaction_graph/clusters/cluster_consensus.fasta",
        ncbi_full="projects/{project}/reports/repeat_interaction_graph/cluster_consensus_vs_ncbi.full_x3.tsv",
        ncbi_best="projects/{project}/reports/repeat_interaction_graph/cluster_consensus_vs_ncbi.best_x3.tsv",
        cluster_ann="projects/{project}/reports/repeat_interaction_graph/cluster_annotations.tsv",
        nodes_enriched="projects/{project}/reports/repeat_interaction_graph/nodes.enriched.tsv",
        html="projects/{project}/reports/repeat_interaction_graph/graph.html",
        token=touch("projects/{project}/reports/repeat_interaction_graph/graph_report.done"),
    params:
        enabled=lambda w: bool(get_param(w.project, "post_tarean_params", "graph_report", "enabled")),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    conda:
        "../envs/reportr.yaml"
    threads: 8
    log:
        f"{LOG_DIR}/project_repeat_interaction_graph_report_{{project}}.log"
    benchmark:
        "benchmarks/repeat_interaction_graph_report_{project}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1

        if [[ "{params.enabled}" != "True" ]]; then
          mkdir -p "$(dirname "{output.token}")"
          : > "{output.queries}"
          : > "{output.index}"
          : > "{output.nodes}"
          : > "{output.hsps}"
          : > "{output.edges_full}"
          : > "{output.edges}"
          : > "{output.clusters}"
          : > "{output.cluster_sizes}"
          mkdir -p "$(dirname "{output.cluster_consensus}")"
          : > "{output.cluster_consensus}"
          : > "{output.ncbi_full}"
          : > "{output.ncbi_best}"
          : > "{output.cluster_ann}"
          : > "{output.nodes_enriched}"
          : > "{output.html}"
          touch "{output.token}"
          exit 0
        fi

        # Per-job tmp under the project (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/tmp/jobs/project_repeat_interaction_graph_report"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): GRAPH_REPORT_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "projects/{wildcards.project}/reports/repeat_interaction_graph"
        mkdir -p "projects/{wildcards.project}/reports/repeat_interaction_graph/clusters"

        python3 workflows/smk_scripts/graph_report_build_queries_latest.py \
          --project "{wildcards.project}" \
          --out-fasta "{output.queries}" \
          --out-index-tsv "{output.index}" \
          --rich-table {input.rich_tables}

        # If there are no ready rich tables, produce an empty report and exit cleanly.
        if [[ ! -s "{output.queries}" ]]; then
          : > "{output.nodes}"
          : > "{output.hsps}"
          : > "{output.edges_full}"
          : > "{output.edges}"
          : > "{output.clusters}"
          : > "{output.cluster_sizes}"
          : > "{output.cluster_consensus}"
          : > "{output.ncbi_full}"
          : > "{output.ncbi_best}"
          : > "{output.cluster_ann}"
          : > "{output.nodes_enriched}"
          : > "{output.html}"
          touch "{output.token}"
          exit 0
        fi

        python3 workflows/smk_scripts/graph_report_build_nodes_table.py \
          --index-tsv "{output.index}" \
          --out-nodes-tsv "{output.nodes}"

        python3 workflows/smk_scripts/graph_report_build_local_edges_cached.py \
          --project "{wildcards.project}" \
          --query-fasta "{output.queries}" \
          --query-index-tsv "{output.index}" \
          --localdb-prefix "projects/{wildcards.project}/blast_db/multifasta.fasta" \
          --out-hsps-tsv "{output.hsps}" \
          --out-edges-full-tsv "{output.edges_full}" \
          --out-edges-filtered-tsv "{output.edges}" \
          --threads "{threads}"

        python3 workflows/smk_scripts/graph_report_cluster_components.py \
          --edges-tsv "{output.edges}" \
          --nodes-fasta "{output.queries}" \
          --out-clusters-tsv "{output.clusters}" \
          --out-cluster-sizes-tsv "{output.cluster_sizes}"

        python3 workflows/smk_scripts/graph_report_align_clusters.py \
          --project "{wildcards.project}" \
          --queries-fasta "{output.queries}" \
          --clusters-tsv "{output.clusters}" \
          --out-dir "projects/{wildcards.project}/reports/repeat_interaction_graph/clusters" \
          --out-consensus-fasta "{output.cluster_consensus}"

        # Cluster consensus annotation (x1+x3 + task ladder) against NCBI repeats, with oligo evidence.
        python3 workflows/smk_scripts/satminer_blast_x3.py \
          --project "{wildcards.project}" \
          --query-fasta "{output.cluster_consensus}" \
          --blast-db-prefix "projects/{wildcards.project}/blast_db/ncbi_repeats" \
          --oligo-db-prefix "projects/{wildcards.project}/blast_db/oligo" \
          --out-full-tsv "{output.ncbi_full}" \
          --out-best-tsv "{output.ncbi_best}" \
          --threads "{threads}"

        python3 workflows/smk_scripts/graph_report_annotate_clusters.py \
          --project "{wildcards.project}" \
          --clusters-tsv "{output.clusters}" \
          --ncbi-best-tsv "{output.ncbi_best}" \
          --oligo-best-tsv "{output.ncbi_best}" \
          --out-tsv "{output.cluster_ann}"

        conda run -n reportr_graph Rscript workflows/smk_scripts/graph_report_render_visnetwork.R \
          "{output.nodes}" \
          "{output.edges_full}" \
          "{output.clusters}" \
          "{output.cluster_ann}" \
          "{output.html}" \
          "{output.nodes_enriched}"

        touch "{output.token}"
        """


rule post_tarean_satminer_blast_x3:
    """BLAST tagged consensus (x1+x3 + task ladder) against the project's NCBI repeats BLAST DB."""
    input:
        query="projects/{project}/samples/{sample}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        db_nsq="projects/{project}/blast_db/ncbi_repeats.nsq",
        oligo_nsq="projects/{project}/blast_db/oligo.nsq",
    output:
        tsv="projects/{project}/samples/{sample}/post_tarean/satminer/compare/consensus_vs_project_ncbi.blast_x3.tsv",
        best_tsv="projects/{project}/samples/{sample}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/compare/BLAST_DONE",
    params:
        raw_dir="projects/{project}/samples/{sample}/post_tarean/satminer/compare/raw",
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    conda:
        "../envs/reportr.yaml"
    threads: 8
    log:
        f"{LOG_DIR}/post_tarean_satminer_blast_x3_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_blast_x3_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_satminer_blast_x3"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): POST_TAREAN_BLAST_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "$(dirname "{output.tsv}")"
        mkdir -p "{params.raw_dir}"
        python3 workflows/smk_scripts/satminer_blast_x3.py \
          --project "{wildcards.project}" \
          --query-fasta "{input.query}" \
          --blast-db-prefix "projects/{wildcards.project}/blast_db/ncbi_repeats" \
          --oligo-db-prefix "projects/{wildcards.project}/blast_db/oligo" \
          --out-full-tsv "{output.tsv}" \
          --out-best-tsv "{output.best_tsv}" \
          --raw-hsps-dir "{params.raw_dir}" \
          --threads "{threads}"
        touch "{output.token}"
        """


rule post_tarean_satminer_parse_rm_out:
    """Parse reads-vs-consensus RepeatMasker output into abundance/divergence TSVs and plots."""
    input:
        rm_out="projects/{project}/samples/{sample}/post_tarean/satminer/repeatmasker_reads_vs_consensus/prepared_forRE.fasta.out",
    output:
        hits="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/samples/{sample}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/RM_OUT_PARSED",
    params:
        plots_dir="projects/{project}/samples/{sample}/post_tarean/satminer/plots",
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_satminer_parse_rm_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_parse_rm_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_satminer_parse_rm_out"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): POST_TAREAN_PARSE_RM_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        python3 workflows/smk_scripts/parse_repeatmasker_out_hits.py \
          --rm-out "{input.rm_out}" \
          --hits-tsv "{output.hits}" \
          --divergence-tsv "{output.div}" \
          --plots-dir "{params.plots_dir}"
        touch "{output.token}"
        """


rule post_tarean_satminer_render_report:
    """Render a lightweight HTML report (tables + plots)."""
    input:
        blast_top="projects/{project}/samples/{sample}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        hits="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/samples/{sample}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        rm_token="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/RM_OUT_PARSED",
        blast_token="projects/{project}/samples/{sample}/post_tarean/satminer/compare/BLAST_DONE",
    output:
        html="projects/{project}/samples/{sample}/post_tarean/satminer/report/{sample}_satminer_report.html",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/report/satminer_compare.done",
    params:
        plots_dir="projects/{project}/samples/{sample}/post_tarean/satminer/plots",
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_satminer_render_report_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_render_report_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_satminer_render_report"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): POST_TAREAN_RENDER_REPORT_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        python3 workflows/smk_scripts/render_satminer_html_report.py \
          --sample "{wildcards.sample}" \
          --out-html "{output.html}" \
          --blast-top-tsv "{input.blast_top}" \
          --hits-tsv "{input.hits}" \
          --div-tsv "{input.div}" \
          --plots-dir "{params.plots_dir}"
        touch "{output.token}"
        """


rule post_tarean_satminer_build_rich_table:
    """Build a unified per-consensus TSV (tags + image path + abundance/divergence + BLAST)."""
    input:
        fasta="projects/{project}/samples/{sample}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        reads_fasta="projects/{project}/samples/{sample}/tarean1/prepared_forRE.fasta",
        hits="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/samples/{sample}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        blast_best="projects/{project}/samples/{sample}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        legacy_summary=lambda w: (
            f"projects/{w.project}/samples/{w.sample}/post_tarean/legacy/{w.sample}_legacy_blast_summary.tsv"
            if (
                get_param(w.project, "post_tarean_params", "legacy_enabled")
                and get_param(w.project, "post_tarean_params", "legacy_merge_into_satminer")
            )
            else []
        ),
        rm_token="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/RM_OUT_PARSED",
        blast_token="projects/{project}/samples/{sample}/post_tarean/satminer/compare/BLAST_DONE",
    output:
        table="projects/{project}/samples/{sample}/post_tarean/satminer/report/rich_table.tsv",
        top="projects/{project}/samples/{sample}/post_tarean/satminer/report/rich_table.top.tsv",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/report/RICH_TABLE_DONE",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_satminer_build_rich_table_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_build_rich_table_{project}_{sample}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_satminer_build_rich_table"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): POST_TAREAN_BUILD_RICH_TABLE_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        python3 workflows/smk_scripts/satminer_build_rich_report_table.py \
          --project "{wildcards.project}" \
          --sample "{wildcards.sample}" \
          --consensus-fasta "{input.fasta}" \
          --reads-fasta "{input.reads_fasta}" \
          --hits-tsv "{input.hits}" \
          --divergence-tsv "{input.div}" \
          --blast-best-tsv "{input.blast_best}" \
          --legacy-blast-summary-tsv "{input.legacy_summary}" \
          --out-tsv "{output.table}" \
          --out-top-tsv "{output.top}" \
          --top-n 30
        touch "{output.token}"
        """


rule post_tarean_satminer_build_rich_table_core:
    """Build rich table without legacy merge (for project LocalDB construction)."""
    input:
        fasta="projects/{project}/samples/{sample}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        reads_fasta="projects/{project}/samples/{sample}/tarean1/prepared_forRE.fasta",
        hits="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/samples/{sample}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        blast_best="projects/{project}/samples/{sample}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        rm_token="projects/{project}/samples/{sample}/post_tarean/satminer/abundance/RM_OUT_PARSED",
        blast_token="projects/{project}/samples/{sample}/post_tarean/satminer/compare/BLAST_DONE",
    output:
        table="projects/{project}/samples/{sample}/post_tarean/satminer/report/rich_table.core.tsv",
        top="projects/{project}/samples/{sample}/post_tarean/satminer/report/rich_table.core.top.tsv",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/report/RICH_TABLE_CORE_DONE",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_satminer_build_rich_table_core_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_build_rich_table_core_{project}_{sample}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_satminer_build_rich_table_core"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): POST_TAREAN_BUILD_RICH_TABLE_CORE_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        python3 workflows/smk_scripts/satminer_build_rich_report_table.py \
          --project "{wildcards.project}" \
          --sample "{wildcards.sample}" \
          --consensus-fasta "{input.fasta}" \
          --reads-fasta "{input.reads_fasta}" \
          --hits-tsv "{input.hits}" \
          --divergence-tsv "{input.div}" \
          --blast-best-tsv "{input.blast_best}" \
          --out-tsv "{output.table}" \
          --out-top-tsv "{output.top}" \
          --top-n 30
        touch "{output.token}"
        """


rule post_tarean_satminer_rich_reports:
    """Render XLSX+DOCX reports with embedded graph images."""
    input:
        table="projects/{project}/samples/{sample}/post_tarean/satminer/report/rich_table.tsv",
        top="projects/{project}/samples/{sample}/post_tarean/satminer/report/rich_table.top.tsv",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/report/RICH_TABLE_DONE",
    output:
        xlsx="projects/{project}/samples/{sample}/post_tarean/satminer/report/{sample}_satminer_report.xlsx",
        docx="projects/{project}/samples/{sample}/post_tarean/satminer/report/{sample}_satminer_report.docx",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/report/satminer_rich_report.done",
    params:
        organism=lambda w: _sample_attr(w.project, w.sample, "organism"),
        genomes=lambda w: _sample_attr(w.project, w.sample, "genomes"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_satminer_rich_reports_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_rich_reports_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_satminer_rich_reports"
        mkdir -p "$tmp_base"
        job_tmp="$(mktemp -d "$tmp_base/run_XXXXXXXX")"
        _cleanup_job_tmp() {{
          local exit_code="$?"
          if [[ "$exit_code" -eq 0 ]]; then
            rm -rf "$job_tmp" || true
            if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
              rmdir "$tmp_base" 2>/dev/null || true
              rmdir "$(dirname "$tmp_base")" 2>/dev/null || true
              rmdir "$(dirname "$(dirname "$tmp_base")")" 2>/dev/null || true
            fi
          else
            echo "$(date): POST_TAREAN_RICH_REPORTS_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        python3 workflows/smk_scripts/satminer_write_xlsx.py \
          --table-tsv "{input.table}" \
          --out-xlsx "{output.xlsx}" \
          --sample "{wildcards.sample}" \
          --org "{params.organism}" \
          --genomes "{params.genomes}"

        python3 workflows/smk_scripts/satminer_write_docx.py \
          --table-tsv "{input.table}" \
          --out-docx "{output.docx}" \
          --sample "{wildcards.sample}" \
          --org "{params.organism}" \
          --genomes "{params.genomes}" \
          --embed-top-images 10
        touch "{output.token}"
        """
