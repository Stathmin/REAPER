# Post-tarean satMiner rules
#
# Write all outputs under:
#   projects/<project>/samples/<sample>/post_tarean/satminer/
#
# This module is iteration-aware: it collects TAREAN consensus FASTAs from
# tarean1..N and tags headers with the iteration depth.

import os
from pathlib import Path


def _iter_depth(project: str) -> int:
    return int(get_param(project, "iterative_assembly", "depth"))

def _sample_consensus_inputs(project: str, sample: str):
    # Collect only completed iterations. This avoids hard-failing when depth is N
    # but some intermediate iterations were not produced yet.
    base = Path(f"projects/{project}/samples/{sample}")
    # Critical: do not return an empty list. An empty input makes the rule runnable
    # immediately (and it will then fail if no TAREAN consensus files exist yet).
    #
    # Gate consensus collection on at least one successful assembly (iter1 in the
    # iterative layout; legacy `tarean/` in non-iterative projects).
    required = base / "tarean1" / "tarean.done"
    legacy = base / "tarean" / "tarean.done"
    if required.exists():
        return [str(required)]
    if legacy.exists():
        return [str(legacy)]
    # If neither exists yet, keep the dependency explicit so Snakemake schedules
    # upstream work rather than running this rule early.
    return [str(required)]


def _sample_attr(project: str, sample: str, key: str) -> str:
    proj = config.get("projects", {}).get(project, {})
    smpl = proj.get("samples", {}).get(sample, {})
    v = smpl.get(key, "")
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return ",".join(str(x) for x in v)
    return str(v)


rule post_tarean_collect_tarean_consensus:
    """Collect and tag all TAREAN consensus FASTAs across tarean1..N for a sample."""
    input:
        lambda w: _sample_consensus_inputs(w.project, w.sample)
    output:
        fasta="projects/{project}/samples/{sample}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/CONSENSUS_READY",
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_collect_consensus_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_collect_consensus_{project}_{sample}.tsv"
    params:
        organism=lambda w: _sample_attr(w.project, w.sample, "organism"),
        genomes=lambda w: _sample_attr(w.project, w.sample, "genomes"),
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        start_ts="$(date +%s)"
        python3 workflows/smk_scripts/smk_log.py header \
          --rule post_tarean_collect_tarean_consensus --category wrapper \
          --kv project="{wildcards.project}" --kv sample="{wildcards.sample}"
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_collect_tarean_consensus"
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
            echo "$(date): POST_TAREAN_CONSENSUS_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"
        export PYTHONHASHSEED="{params.pythonhashseed}"

        mkdir -p "$(dirname "{output.fasta}")"
        consensus_files=()
        for d in projects/{wildcards.project}/samples/{wildcards.sample}/tarean*/ ; do
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
          echo "ERROR: no non-empty TAREAN_consensus_rank_*.fasta found under projects/{wildcards.project}/samples/{wildcards.sample}/tarean*/" >&2
          exit 2
        fi
        python3 workflows/smk_scripts/concat_tarean_consensus_tagged.py \
          --out "{output.fasta}" \
          --sample "{wildcards.sample}" \
          --organism "{params.organism}" \
          --genomes "{params.genomes}" \
          --inputs "${{consensus_files[@]}}"
        # Record output counts for provenance.
        python3 workflows/smk_scripts/smk_log.py summarize-fasta --path "{output.fasta}"
        touch "{output.token}"
        elapsed="$(( $(date +%s) - start_ts ))"
        python3 workflows/smk_scripts/smk_log.py footer --ok 1 --elapsed-s "$elapsed" \
          --kv out="{output.fasta}"
        """


rule post_tarean_filtered_reads_fasta:
    """Convert filtered_reads FASTQs to a single FASTA for RepeatMasker."""
    input:
        r1="projects/{project}/samples/{sample}/filtered_reads/R1.fq",
        r2="projects/{project}/samples/{sample}/filtered_reads/R2.fq",
    output:
        fasta=temp("projects/{project}/samples/{sample}/post_tarean/satminer/genome/filtered_reads.fasta"),
        token=temp("projects/{project}/samples/{sample}/post_tarean/satminer/genome/FILTERED_READS_FASTA_READY"),
    conda:
        "../envs/reportr.yaml"
    threads: 2
    log:
        f"{LOG_DIR}/post_tarean_filtered_reads_fasta_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_filtered_reads_fasta_{project}_{sample}.tsv"
    params:
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_filtered_reads_fasta"
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
            echo "$(date): POST_TAREAN_FILTERED_READS_FASTA_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"
        export PYTHONHASHSEED="{params.pythonhashseed}"

        outdir="$(dirname "{output.fasta}")"
        mkdir -p "$outdir"
        # IMPORTANT:
        # - Use .fa extension so reformat writes FASTA (not FASTQ).
        # - Trim read descriptions to avoid very long IDs that RepeatMasker rejects.
        tmp1="$outdir/R1.fa"
        tmp2="$outdir/R2.fa"
        reformat.sh in="{input.r1}" out="$tmp1" fastawrap=0 overwrite=t trimreaddescription=t
        reformat.sh in="{input.r2}" out="$tmp2" fastawrap=0 overwrite=t trimreaddescription=t
        cat "$tmp1" "$tmp2" > "{output.fasta}"
        rm -f "$tmp1" "$tmp2"
        touch "{output.token}"
        """


rule post_tarean_repeatmasker_filtered_reads_vs_consensus:
    """RepeatMasker on full filtered_reads (FASTA) using combined consensus as custom library."""
    input:
        genome="projects/{project}/samples/{sample}/post_tarean/satminer/genome/filtered_reads.fasta",
        genome_token="projects/{project}/samples/{sample}/post_tarean/satminer/genome/FILTERED_READS_FASTA_READY",
        lib="projects/{project}/samples/{sample}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/CONSENSUS_READY",
    output:
        rm_out="projects/{project}/samples/{sample}/post_tarean/satminer/repeatmasker_filtered_reads_vs_consensus/filtered_reads.fasta.out",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/REPEATMASKER_FILTERED_READY",
    params:
        threads=lambda w: int(get_param(w.project, "satminer", "repeatmasker_threads", sample_id=w.sample)),  # type: ignore[name-defined]
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    threads:
        lambda w: int(get_param(w.project, "satminer", "repeatmasker_threads", sample_id=w.sample))  # type: ignore[name-defined]
    resources:
        repeatmasker_slots=1
    conda:
        "../envs/reportr.yaml"
    log:
        f"{LOG_DIR}/post_tarean_repeatmasker_filtered_reads_vs_consensus_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_repeatmasker_filtered_reads_vs_consensus_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_repeatmasker_filtered_reads_vs_consensus"
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
            echo "$(date): POST_TAREAN_RM_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"
        export PYTHONHASHSEED="{params.pythonhashseed}"

        outdir="projects/{wildcards.project}/samples/{wildcards.sample}/post_tarean/satminer/repeatmasker_filtered_reads_vs_consensus"
        mkdir -p "$outdir"
        python3 workflows/smk_scripts/satminer_repeatmasker_customlib.py \
          --genome-fasta "{input.genome}" \
          --lib-fasta "{input.lib}" \
          --outdir "$outdir" \
          --threads "{params.threads}"
        test -f "{output.rm_out}"
        touch "{output.token}"
        """


rule post_tarean_repeatmasker_reads_vs_consensus:
    """RepeatMasker on tarean1 prepared reads using combined consensus as custom library."""
    input:
        genome="projects/{project}/samples/{sample}/tarean1/prepared_forRE.fasta",
        lib="projects/{project}/samples/{sample}/post_tarean/satminer/consensus_all_iters.tagged.fasta",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/CONSENSUS_READY",
    output:
        rm_out="projects/{project}/samples/{sample}/post_tarean/satminer/repeatmasker_reads_vs_consensus/prepared_forRE.fasta.out",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/REPEATMASKER_READY",
    params:
        threads=lambda w: int(get_param(w.project, "satminer", "repeatmasker_threads", sample_id=w.sample)),  # type: ignore[name-defined]
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    threads:
        lambda w: int(get_param(w.project, "satminer", "repeatmasker_threads", sample_id=w.sample))  # type: ignore[name-defined]
    resources:
        repeatmasker_slots=1
    conda:
        "../envs/reportr.yaml"
    log:
        f"{LOG_DIR}/post_tarean_repeatmasker_reads_vs_consensus_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_repeatmasker_reads_vs_consensus_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_repeatmasker_reads_vs_consensus"
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
            echo "$(date): POST_TAREAN_RM_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"
        export PYTHONHASHSEED="{params.pythonhashseed}"

        outdir="projects/{wildcards.project}/samples/{wildcards.sample}/post_tarean/satminer/repeatmasker_reads_vs_consensus"
        mkdir -p "$outdir"
        python3 workflows/smk_scripts/satminer_repeatmasker_customlib.py \
          --genome-fasta "{input.genome}" \
          --lib-fasta "{input.lib}" \
          --outdir "$outdir" \
          --threads "{params.threads}"
        test -f "{output.rm_out}"
        touch "{output.token}"
        """


rule post_tarean_rm_getseq_reads_vs_consensus:
    """Extract RepeatMasker hit sequences from reads-vs-consensus run."""
    input:
        genome="projects/{project}/samples/{sample}/tarean1/prepared_forRE.fasta",
        rm_out="projects/{project}/samples/{sample}/post_tarean/satminer/repeatmasker_reads_vs_consensus/prepared_forRE.fasta.out",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/REPEATMASKER_READY",
    output:
        hits="projects/{project}/samples/{sample}/post_tarean/satminer/rm_getseq/prepared_forRE.fasta.out.fas",
        token="projects/{project}/samples/{sample}/post_tarean/satminer/RM_GETSEQ_DONE",
    params:
        min_len=lambda w: int(get_param(w.project, "satminer", "rm_getseq_min_len", sample_id=w.sample)),  # type: ignore[name-defined]
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_rm_getseq_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_rm_getseq_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        # Per-job tmp under the sample (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_rm_getseq_reads_vs_consensus"
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
            echo "$(date): POST_TAREAN_RM_GETSEQ_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT
        export TMPDIR="$job_tmp"
        export TMP="$job_tmp"
        export TEMP="$job_tmp"
        export PYTHONHASHSEED="{params.pythonhashseed}"

        mkdir -p "$(dirname "{output.hits}")"
        python3 workflows/smk_scripts/satminer_rm_getseq.py \
          --fasta "{input.genome}" \
          --rm-out "{input.rm_out}" \
          --min-len "{params.min_len}" \
          --out "{output.hits}"
        touch "{output.token}"
        """


rule post_tarean_satminer_done:
    """Convenience per-sample token for post-tarean satMiner outputs."""
    input:
        "projects/{project}/samples/{sample}/post_tarean/satminer/RM_GETSEQ_DONE"
    output:
        touch("projects/{project}/samples/{sample}/post_tarean/satminer/satminer.done")
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_satminer_done_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_satminer_done_{project}_{sample}.tsv"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): satMiner token ready" >&2
        touch "{output}"
        """


rule post_tarean_legacy_reports:
    """Optional legacy post_tarean reports (config-driven; tmp cleaned on success only)."""
    input:
        # Ensure at least one iteration exists (and that the sample is real).
        "projects/{project}/samples/{sample}/post_tarean/satminer/CONSENSUS_READY",
        # Ensure canonical LocalDB exists so BLAST 'local' is discoverable.
        "projects/{project}/blast_db/localdb.ready"
    output:
        blast_csv="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_blast_results.csv",
        analysis_csv="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_repeat_analysis.csv",
        analysis_xlsx="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_repeat_analysis.xlsx",
        analysis_docx="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_repeat_analysis.docx",
        confidence_json="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_confidence.json",
        token=touch("projects/{project}/samples/{sample}/post_tarean/legacy/legacy.done"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_legacy_reports_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_legacy_reports_{project}_{sample}.tsv"
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
          --rule post_tarean_legacy_reports --category wrapper \
          --kv project="{wildcards.project}" --kv sample="{wildcards.sample}"

        if [[ "{params.legacy_enabled}" != "True" ]]; then
          echo "Legacy post_tarean disabled (post_tarean_params.legacy_enabled=false); skipping." >&2
          # Produce empty placeholders to satisfy downstream optional merges.
          mkdir -p "$(dirname "{output.blast_csv}")"
          : > "{output.blast_csv}"
          : > "{output.analysis_csv}"
          : > "{output.confidence_json}"
          : > "{output.analysis_xlsx}"
          : > "{output.analysis_docx}"
          touch "{output.token}"
          exit 0
        fi

        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_legacy_reports"
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

        outdir="projects/{wildcards.project}/samples/{wildcards.sample}/post_tarean/legacy"
        mkdir -p "$outdir"

        python3 -m post_tarean.pipeline \
          "{wildcards.sample}" \
          --project-id "{wildcards.project}" \
          --output-dir "$outdir"

        # Legacy pipeline uses index-based filenames.
        test -s "$outdir/{wildcards.sample}_blast_results.csv"
        test -s "$outdir/{wildcards.sample}_repeat_analysis.csv"
        test -s "$outdir/{wildcards.sample}_confidence.json"

        touch "{output.token}"
        elapsed="$(( $(date +%s) - start_ts ))"
        python3 workflows/smk_scripts/smk_log.py footer --ok 1 --elapsed-s "$elapsed" \
          --kv out="$outdir"
        """


rule post_tarean_legacy_blast_summary:
    """Convert legacy blast_results.csv into a satMiner-mergeable per-cluster TSV."""
    input:
        blast_csv="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_blast_results.csv",
        token="projects/{project}/samples/{sample}/post_tarean/legacy/legacy.done",
    output:
        tsv="projects/{project}/samples/{sample}/post_tarean/legacy/{sample}_legacy_blast_summary.tsv",
        token=touch("projects/{project}/samples/{sample}/post_tarean/legacy/LEGACY_BLAST_SUMMARY_DONE"),
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_legacy_blast_summary_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_legacy_blast_summary_{project}_{sample}.tsv"
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

        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tmp/jobs/post_tarean_legacy_blast_summary"
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
