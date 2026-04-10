# Post-tarean satMiner for comparative analyses (iteration-aware)
#
# Produces satMiner-style TSVs + plots + XLSX/DOCX report under:
#   projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/
#
# Key difference vs sample pipeline:
# - Abundance is split into per-sample columns using read ID prefixes in comparative_reads.fasta
#   (e.g. KA1read123_f -> sample prefix KA1).

import os
from pathlib import Path


def _iter_depth(project: str) -> int:
    return int(get_param(project, "iterative_assembly", "depth"))


def _comparative_samples_satminer(project: str, comparative_analysis: str) -> list[str]:
    proj = config.get("projects", {}).get(project, {}) or {}
    ca = (proj.get("comparative_analyses", {}) or {}).get(comparative_analysis, {}) or {}
    samples = ca.get("samples", []) or []
    return list(samples) if isinstance(samples, (list, tuple)) else []


def _sample_prefix(project: str, sample: str) -> str:
    return str(get_sample_prefix(project, sample))


def _comparative_prefixes_csv(project: str, comparative_analysis: str) -> str:
    m = comparative_prefix_map(project, comparative_analysis)
    return ",".join(m[s] for s in _comparative_samples_satminer(project, comparative_analysis))


def _comparative_consensus_inputs(project: str, comparative_analysis: str):
    # Collect only completed iterations. This avoids hard-failing when depth is N
    # but some iterations were not produced yet.
    base = Path(f"projects/{project}/comparative/{comparative_analysis}")
    # Critical: do not return an empty list. An empty input makes the rule runnable
    # immediately (and it will then fail if no TAREAN consensus files exist yet).
    #
    # Gate consensus collection on at least one successful comparative assembly (iter1).
    required = base / "tarean1" / "COMPARATIVE_TAREAN_COMPLETE"
    if required.exists():
        return [str(required)]
    # If it doesn't exist yet, keep the dependency explicit so Snakemake schedules
    # upstream work rather than running this rule early.
    return [str(required)]


rule comparative_post_tarean_collect_tarean_consensus:
    """Collect and tag all TAREAN consensus FASTAs across tarean1..N for a comparative analysis."""
    input:
        lambda w: _comparative_consensus_inputs(w.project, w.comparative_analysis)
    output:
        fasta="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/CONSENSUS_READY",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_collect_consensus_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_collect_consensus_{project}_{comparative_analysis}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        start_ts="$(date +%s)"
        python3 workflows/smk_scripts/smk_log.py header \
          --rule comparative_post_tarean_collect_tarean_consensus --category wrapper \
          --kv project="{wildcards.project}" --kv comparative_analysis="{wildcards.comparative_analysis}"
        # Per-job tmp under the comparative analysis (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tmp/jobs/comparative_post_tarean_collect_tarean_consensus"
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
            echo "$(date): COMP_POST_TAREAN_CONSENSUS_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        mkdir -p "$(dirname "{output.fasta}")"
        consensus_files=()
        for d in projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tarean*/tarean/ ; do
          if [[ ! -d "$d" ]]; then
            continue
          fi
          for r in 1 2 3 4; do
            f="${{d}}/TAREAN_consensus_rank_${{r}}.fasta"
            if [[ -s "$f" ]]; then
              consensus_files+=("$f")
            fi
          done
        done
        if [[ "${{#consensus_files[@]}}" -eq 0 ]]; then
          echo "ERROR: no non-empty TAREAN_consensus_rank_*.fasta found under projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tarean*/tarean/" >&2
          exit 2
        fi
        python3 workflows/smk_scripts/concat_tarean_consensus_tagged.py \
          --out "{output.fasta}" \
          --sample "{wildcards.comparative_analysis}" \
          --organism "" \
          --genomes "" \
          --inputs "${{consensus_files[@]}}"
        python3 workflows/smk_scripts/smk_log.py summarize-fasta --path "{output.fasta}"
        touch "{output.token}"
        elapsed="$(( $(date +%s) - start_ts ))"
        python3 workflows/smk_scripts/smk_log.py footer --ok 1 --elapsed-s "$elapsed" \
          --kv out="{output.fasta}"
        """


rule comparative_post_tarean_repeatmasker_reads_vs_consensus:
    """RepeatMasker on comparative reads using combined comparative consensus as custom library."""
    input:
        genome="projects/{project}/comparative/{comparative_analysis}/tarean1/comparative_reads.fasta",
        lib="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/CONSENSUS_READY",
    output:
        rm_out="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/repeatmasker_reads_vs_consensus/comparative_reads.fasta.out",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/REPEATMASKER_READY",
    params:
        threads=lambda w: int(
            get_param(w.project, "satminer", "repeatmasker_threads", comparative_context=True)
        ),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    threads:
        lambda w: int(
            get_param(w.project, "satminer", "repeatmasker_threads", comparative_context=True)
        )
    resources:
        repeatmasker_slots=1
    conda:
        "../envs/reportr.yaml"
    log:
        f"{LOG_DIR}/comparative_post_tarean_repeatmasker_reads_vs_consensus_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_repeatmasker_reads_vs_consensus_{project}_{comparative_analysis}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the comparative analysis (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tmp/jobs/comparative_post_tarean_repeatmasker_reads_vs_consensus"
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
            echo "$(date): COMP_POST_TAREAN_RM_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        outdir="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/post_tarean/satminer/repeatmasker_reads_vs_consensus"
        mkdir -p "$outdir"
        python3 workflows/smk_scripts/satminer_repeatmasker_customlib.py \
          --genome-fasta "{input.genome}" \
          --lib-fasta "{input.lib}" \
          --outdir "$outdir" \
          --threads "{params.threads}"
        test -f "{output.rm_out}"
        touch "{output.token}"
        """


rule comparative_post_tarean_satminer_parse_rm_out:
    """Parse RepeatMasker output into abundance/divergence TSVs, split by read prefix (per-sample columns)."""
    input:
        rm_out="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/repeatmasker_reads_vs_consensus/comparative_reads.fasta.out",
        rm_token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/REPEATMASKER_READY",
    output:
        hits="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/abundance/RM_OUT_PARSED",
    params:
        plots_dir="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/plots",
        prefixes_csv=lambda w: _comparative_prefixes_csv(w.project, w.comparative_analysis),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_parse_rm_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_parse_rm_{project}_{comparative_analysis}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        python3 workflows/smk_scripts/parse_repeatmasker_out_hits_prefix.py \
          --rm-out "{input.rm_out}" \
          --hits-tsv "{output.hits}" \
          --divergence-tsv "{output.div}" \
          --plots-dir "{params.plots_dir}" \
          --prefixes "{params.prefixes_csv}"
        touch "{output.token}"
        """


rule comparative_post_tarean_satminer_blast_x3:
    """BLAST comparative tagged consensus against the project's NCBI repeats BLAST DB."""
    input:
        query="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        db_nsq="projects/{project}/blast_db/ncbi_repeats.nsq",
        oligo_nsq="projects/{project}/blast_db/oligo.nsq",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/CONSENSUS_READY",
    output:
        tsv="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/consensus_vs_project_ncbi.blast_x3.tsv",
        best_tsv="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/BLAST_DONE",
    params:
        raw_dir="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/raw",
    conda:
        "../envs/reportr.yaml"
    threads: 8
    log:
        f"{LOG_DIR}/comparative_post_tarean_satminer_blast_x3_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_blast_x3_{project}_{comparative_analysis}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
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


rule comparative_post_tarean_satminer_build_rich_table:
    """Build a unified per-consensus TSV with per-sample abundance columns (read prefixes)."""
    input:
        fasta="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        reads_fasta="projects/{project}/comparative/{comparative_analysis}/tarean1/comparative_reads.fasta",
        hits="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        blast_best="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        legacy_summary=lambda w: (
            f"projects/{w.project}/comparative/{w.comparative_analysis}/post_tarean/legacy/{w.comparative_analysis}_legacy_blast_summary.tsv"
            if (
                get_param(w.project, "post_tarean_params", "legacy_enabled")
                and get_param(w.project, "post_tarean_params", "legacy_merge_into_satminer")
            )
            else []
        ),
        rm_token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/abundance/RM_OUT_PARSED",
        blast_token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/BLAST_DONE",
    output:
        table="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/rich_table.tsv",
        top="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/rich_table.top.tsv",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/RICH_TABLE_DONE",
    params:
        prefixes_csv=lambda w: _comparative_prefixes_csv(w.project, w.comparative_analysis),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_satminer_build_rich_table_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_build_rich_table_{project}_{comparative_analysis}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        python3 workflows/smk_scripts/satminer_build_rich_report_table_comparative.py \
          --project "{wildcards.project}" \
          --comparative-analysis "{wildcards.comparative_analysis}" \
          --consensus-fasta "{input.fasta}" \
          --reads-fasta "{input.reads_fasta}" \
          --hits-tsv "{input.hits}" \
          --divergence-tsv "{input.div}" \
          --blast-best-tsv "{input.blast_best}" \
          --legacy-blast-summary-tsv "{input.legacy_summary}" \
          --prefixes "{params.prefixes_csv}" \
          --out-tsv "{output.table}" \
          --out-top-tsv "{output.top}" \
          --top-n 30
        touch "{output.token}"
        """


rule comparative_post_tarean_satminer_build_rich_table_core:
    """Build rich table without legacy merge (for project LocalDB construction)."""
    input:
        fasta="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        reads_fasta="projects/{project}/comparative/{comparative_analysis}/tarean1/comparative_reads.fasta",
        hits="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/abundance/reads_vs_consensus_hits.tsv",
        div="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/divergence/reads_vs_consensus_divergence.tsv",
        blast_best="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/consensus_vs_project_ncbi.best_x3.tsv",
        rm_token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/abundance/RM_OUT_PARSED",
        blast_token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/compare/BLAST_DONE",
    output:
        table="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/rich_table.core.tsv",
        top="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/rich_table.core.top.tsv",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/RICH_TABLE_CORE_DONE",
    params:
        prefixes_csv=lambda w: _comparative_prefixes_csv(w.project, w.comparative_analysis),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_satminer_build_rich_table_core_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_build_rich_table_core_{project}_{comparative_analysis}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        python3 workflows/smk_scripts/satminer_build_rich_report_table_comparative.py \
          --project "{wildcards.project}" \
          --comparative-analysis "{wildcards.comparative_analysis}" \
          --consensus-fasta "{input.fasta}" \
          --reads-fasta "{input.reads_fasta}" \
          --hits-tsv "{input.hits}" \
          --divergence-tsv "{input.div}" \
          --blast-best-tsv "{input.blast_best}" \
          --prefixes "{params.prefixes_csv}" \
          --out-tsv "{output.table}" \
          --out-top-tsv "{output.top}" \
          --top-n 30
        touch "{output.token}"
        """


rule comparative_post_tarean_satminer_rich_reports:
    """Render XLSX+DOCX reports for comparative satMiner table (includes per-sample columns)."""
    input:
        table="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/rich_table.tsv",
        top="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/rich_table.top.tsv",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/RICH_TABLE_DONE",
    output:
        xlsx="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/{comparative_analysis}_satminer_report.xlsx",
        docx="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/{comparative_analysis}_satminer_report.docx",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/report/satminer_rich_report.done",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_satminer_rich_reports_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_rich_reports_{project}_{comparative_analysis}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        python3 workflows/smk_scripts/satminer_write_xlsx.py \
          --table-tsv "{input.table}" \
          --out-xlsx "{output.xlsx}" \
          --sample "{wildcards.comparative_analysis}" \
          --org "" \
          --genomes ""

        python3 workflows/smk_scripts/satminer_write_docx.py \
          --table-tsv "{input.table}" \
          --out-docx "{output.docx}" \
          --sample "{wildcards.comparative_analysis}" \
          --org "" \
          --genomes "" \
          --embed-top-images 10
        touch "{output.token}"
        """


rule comparative_post_tarean_legacy_reports:
    """Optional legacy post_tarean reports for comparative analysis (config-driven)."""
    input:
        "projects/{project}/comparative/{comparative_analysis}/post_tarean/satminer/CONSENSUS_READY",
        "projects/{project}/blast_db/localdb.ready"
    output:
        blast_csv="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_blast_results.csv",
        analysis_csv="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_repeat_analysis.csv",
        analysis_xlsx="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_repeat_analysis.xlsx",
        analysis_docx="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_repeat_analysis.docx",
        confidence_json="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_confidence.json",
        token=touch("projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/legacy.done"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_legacy_reports_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_legacy_reports_{project}_{comparative_analysis}.tsv"
    params:
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
        legacy_enabled=lambda w: get_param(w.project, "post_tarean_params", "legacy_enabled"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        start_ts="$(date +%s)"
        python3 workflows/smk_scripts/smk_log.py header \
          --rule comparative_post_tarean_legacy_reports --category wrapper \
          --kv project="{wildcards.project}" --kv comparative_analysis="{wildcards.comparative_analysis}"

        if [[ "{params.legacy_enabled}" != "True" ]]; then
          echo "Legacy post_tarean disabled (post_tarean_params.legacy_enabled=false); skipping." >&2
          mkdir -p "$(dirname "{output.blast_csv}")"
          : > "{output.blast_csv}"
          : > "{output.analysis_csv}"
          : > "{output.confidence_json}"
          : > "{output.analysis_xlsx}"
          : > "{output.analysis_docx}"
          touch "{output.token}"
          exit 0
        fi

        tmp_base="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tmp/jobs/comparative_post_tarean_legacy_reports"
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
            echo "$(date): POST_TAREAN_LEGACY_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"
        export PYTHONHASHSEED="{params.pythonhashseed}"
        export REPORTR_GLOBAL_CONFIG="projects/global_config.yaml"

        outdir="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/post_tarean/legacy"
        mkdir -p "$outdir"

        python3 -m post_tarean.pipeline \
          "{wildcards.comparative_analysis}" \
          --project-id "{wildcards.project}" \
          --output-dir "$outdir"

        test -s "$outdir/{wildcards.comparative_analysis}_blast_results.csv"
        test -s "$outdir/{wildcards.comparative_analysis}_repeat_analysis.csv"
        test -s "$outdir/{wildcards.comparative_analysis}_confidence.json"

        touch "{output.token}"
        elapsed="$(( $(date +%s) - start_ts ))"
        python3 workflows/smk_scripts/smk_log.py footer --ok 1 --elapsed-s "$elapsed" \
          --kv out="$outdir"
        """


rule comparative_post_tarean_legacy_blast_summary:
    """Convert legacy blast_results.csv into a satMiner-mergeable per-cluster TSV (comparative)."""
    input:
        blast_csv="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_blast_results.csv",
        token="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/legacy.done",
    output:
        tsv="projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/{comparative_analysis}_legacy_blast_summary.tsv",
        token=touch("projects/{project}/comparative/{comparative_analysis}/post_tarean/legacy/LEGACY_BLAST_SUMMARY_DONE"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/comparative_post_tarean_legacy_blast_summary_{{project}}_{{comparative_analysis}}.log"
    benchmark:
        "benchmarks/post_tarean_legacy_blast_summary_{project}_{comparative_analysis}.tsv"
    params:
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
        legacy_enabled=lambda w: get_param(w.project, "post_tarean_params", "legacy_enabled"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        if [[ "{params.legacy_enabled}" != "True" ]]; then
          mkdir -p "$(dirname "{output.tsv}")"
          printf "legacy_cluster_key\tlegacy_best_sseqid\tlegacy_best_evalue\tlegacy_best_pident\tlegacy_best_length\tlegacy_best_coverage\tlegacy_coverage_type\tlegacy_task_type\tlegacy_hits_topN\n" > "{output.tsv}"
          touch "{output.token}"
          exit 0
        fi

        tmp_base="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tmp/jobs/comparative_post_tarean_legacy_blast_summary"
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
            echo "$(date): POST_TAREAN_LEGACY_SUMMARY_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"

        python3 workflows/smk_scripts/legacy_post_tarean_blast_to_best_tsv.py \
          --in-blast-results-csv "{input.blast_csv}" \
          --out-tsv "{output.tsv}"

        touch "{output.token}"
        """
