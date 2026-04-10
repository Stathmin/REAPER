# Iterative comparative TAREAN assemblies (optional; config-gated)

import json
import os
import subprocess
from pathlib import Path


def _iterative_enabled(project_id: str) -> bool:
    return bool(get_param(project_id, "iterative_assembly", "enabled"))


def _iter_depth(project_id: str) -> int:
    return int(get_param(project_id, "iterative_assembly", "depth"))


def _deconseq_db(project_id: str) -> str:
    return str(get_param(project_id, "deconseq", "db"))


def _deconseq_threads(project_id: str) -> int:
    return int(get_param(project_id, "deconseq", "threads", comparative_context=True))


def _deconseq_removeifeitherbad_t(project_id: str) -> str:
    v = bool(get_param(project_id, "deconseq", "removeifeitherbad", comparative_context=True))
    return "t" if v else "f"


def _deconseq_minkmerhits(project_id: str) -> int:
    return int(get_param(project_id, "deconseq", "minkmerhits", comparative_context=True))


def _deconseq_iterative_reference_coverage_fraction(project_id: str) -> float:
    return float(get_param(project_id, "deconseq", "iterative_reference_coverage_fraction"))


def _cleanup_after_prepare(project_id: str) -> bool:
    return bool(get_param(project_id, "cleanup_after_prepare"))


def _guard_iterative(project_id: str, iter_s: str, *, min_iter: int = 1) -> str:
    if not _iterative_enabled(project_id):
        raise ValueError("Iterative assembly disabled for project (iterative_assembly.enabled=false)")
    it = int(iter_s)
    if it < 1 or it > _iter_depth(project_id):
        raise ValueError("Requested iteration outside configured depth")
    if it < min_iter:
        raise ValueError(f"Rule only valid for iter>={min_iter}")
    return "ok"


def _iterative_comparative_samples(project_id: str, comparative_analysis: str) -> list[str]:
    proj = config.get("projects", {}).get(project_id, {}) or {}
    ca = (proj.get("comparative_analyses", {}) or {}).get(comparative_analysis, {}) or {}
    samples = ca.get("samples", []) or []
    return list(samples) if isinstance(samples, (list, tuple)) else []


def _comparative_prefix_length(project_id: str, comparative_analysis: str) -> int:
    base = int(get_tarean_param(project_id, "prefix_length"))
    # If analysis exists, allow analysis-level override for tarean_params.prefix_length.
    if comparative_analysis in (config.get("projects", {}).get(project_id, {}) or {}).get("comparative_analyses", {}):
        base = int(get_tarean_param_for_comparative(project_id, comparative_analysis, "prefix_length"))
    smpls = _iterative_comparative_samples(project_id, comparative_analysis)
    if not smpls:
        return base
    longest = max(
        len(str(config["projects"][project_id]["samples"][s].get("prefix") or str(s).upper()))
        for s in smpls
    )
    return max(base, int(longest))


rule prepare_comparative_reads_iter:
    """Prepare comparative reads for a given iteration under tarean{iter}/."""
    input:
        # Ensures concat_prev_comparative_iter_contigs runs before prep when iter>=2 (same path as params.ref_fasta).
        concat_ref=lambda w: (
            []
            if int(w.iter) < 2
            else f"projects/{w.project}/comparative/{w.comparative_analysis}/tarean{w.iter}/iterative/reference_contigs.fasta"
        ),
    output:
        combined_reads="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/comparative_reads.fasta",
        summary="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/comparative_summary.txt",
        token="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/PREPARATION_COMPLETE",
    params:
        project="{project}",
        comparative_analysis="{comparative_analysis}",
        iter="{iter}",
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=1),
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        filter_threads=lambda w: _deconseq_threads(w.project),
        kmer_k=lambda w: int(get_param(w.project, "deconseq", "kmer_k", comparative_context=True)),
        removeifeitherbad_t=lambda w: _deconseq_removeifeitherbad_t(w.project),
        minkmerhits=lambda w: _deconseq_minkmerhits(w.project),
        ref_fasta=lambda w: (
            ""
            if int(w.iter) == 1
            else f"projects/{w.project}/comparative/{w.comparative_analysis}/tarean{w.iter}/iterative/reference_contigs.fasta"
        ),
        cleanup_after_prepare=lambda w: _cleanup_after_prepare(w.project),
    log:
        f"{LOG_DIR}/comparative_reads_{{project}}_{{comparative_analysis}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/prepare_comparative_reads_{project}_{comparative_analysis}_iter{iter}.tsv"
    threads:
        lambda w: _deconseq_threads(w.project)
    resources:
        bbmap_slots=1,
        # This rule also runs FastQC indirectly via sampling utilities; keep it slot-bounded if enabled.
        fastqc_slots=1,
        prep_slots=1,
    conda:
        "../envs/reportr.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=prepare_comparative_reads_iter CATEGORY=heavy project={wildcards.project} analysis={wildcards.comparative_analysis} iter={wildcards.iter} seed={params.pythonhashseed}"
        outdir="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/tarean{wildcards.iter}"
        mkdir -p "$outdir"

        extra_args=()
        if [[ "{params.ref_fasta}" != "" ]]; then
          extra_args+=(--filter-unmapped-against "{params.ref_fasta}" --filter-threads "{params.filter_threads}")
          extra_args+=(--minkmerhits "{params.minkmerhits}" --removeifeitherbad "{params.removeifeitherbad_t}")
        fi

        python3 workflows/smk_scripts/prepare_comparative_reads.py \
          "{wildcards.project}" \
          "{wildcards.comparative_analysis}" \
          --seed "{params.pythonhashseed}" \
          --outdir "$outdir" \
          --write-summary "{output.summary}" \
          --write-token "{output.token}" \
          --kmer-k "{params.kmer_k}" \
          --keep-sampled-fastq \
          "${{extra_args[@]}}"

        if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
          rm -rf "$outdir/tmp" || true
          rm -rf "$outdir/tarean/tmp" || true
          rm -rf "$outdir/tarean/Rserv" || true
          rmdir "$outdir/tmp" 2>/dev/null || true
        fi
        """


rule run_tarean_comparative_iter:
    """Run TAREAN on comparative_reads.fasta for iteration under tarean{iter}/."""
    input:
        comparative_reads="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/comparative_reads.fasta",
        token="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/PREPARATION_COMPLETE",
    output:
        comparative_complete="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/COMPARATIVE_TAREAN_COMPLETE",
        tarean_log=f"{LOG_DIR}/tarean_comparative_{{project}}_{{comparative_analysis}}_iter{{iter}}.log",
        cluster_table="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/tarean/CLUSTER_TABLE.csv",
    params:
        project="{project}",
        sample="{comparative_analysis}",
        comparative_analysis="{comparative_analysis}",
        reads_per_assembly=lambda w: get_reads_per_assembly(
            w.project,
            comparative_analysis=w.comparative_analysis,
            comparative_context=True,
        ),
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        temp_dir=lambda w: get_param(w.project, "read_preparation", "temp_dir", comparative_context=True),
        assembly_min=lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "assembly_min"),
        mincl_percent=lambda w: get_mincl_percent(w.project),
        min_lcov=lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "min_lcov"),
        merge_threshold=lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "merge_threshold"),
        r_value=lambda w: get_tarean_r_value(w.project),
        options=lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "options"),
        paired=lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "paired")),
        automatic_filtering=lambda w: bool(
            get_tarean_param_for_comparative(w.project, w.comparative_analysis, "automatic_filtering")
        ),
        tarean_mode=lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "tarean_mode")),
        keep_names=lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "keep_names")),
        cleanup=lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "cleanup")),
        domain_search=lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "domain_search"),
        proportions_json=lambda w: json.dumps(
            {
                s: config["projects"][w.project]["samples"][s].get("genome_size", 1.0)
                for s in _iterative_comparative_samples(w.project, w.comparative_analysis)
            }
        ),
        prefix_length=lambda w: _comparative_prefix_length(w.project, w.comparative_analysis),
        tarean_dir="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/tarean",
        rserv_port=lambda w: stable_rserv_port(43113, f"{w.comparative_analysis}_iter{w.iter}"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    log:
        f"{LOG_DIR}/seqclust_{{project}}_{{comparative_analysis}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/run_tarean_{project}_{comparative_analysis}_iter{iter}.tsv"
    threads:
        lambda w: get_seqclust_threads(
            w.project,
            get_tarean_param_for_comparative(w.project, w.comparative_analysis, "threads"),
        )
    resources:
        seqclust_slots=1,
        mem_mb=lambda w: max(1024, int(get_tarean_r_value(w.project) / 1024)),
    conda:
        "../envs/reportr.yaml"
    script:
        "smk_scripts/run_tarean_step.py"


rule select_contigs_comparative_iter_re2:
    """Select contigs from comparative iteration clusters for next iteration filtering."""
    input:
        done="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/COMPARATIVE_TAREAN_COMPLETE",
    output:
        selected_fasta="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/iterative/selected_contigs.fasta",
        token="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/iterative/SELECTED_CONTIGS_READY",
    params:
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=1),
        clusters_dir="projects/{project}/comparative/{comparative_analysis}/tarean{iter}/tarean/seqclust/clustering/clusters",
        coverage_fraction=lambda w: _deconseq_iterative_reference_coverage_fraction(w.project),
    log:
        f"{LOG_DIR}/select_contigs_{{project}}_{{comparative_analysis}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/select_contigs_{project}_{comparative_analysis}_iter{iter}.tsv"
    conda:
        "../envs/reportr.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=select_contigs_comparative_iter_re2 CATEGORY=wrapper project={wildcards.project} analysis={wildcards.comparative_analysis} iter={wildcards.iter}"
        mkdir -p "$(dirname "{output.selected_fasta}")"
        tmp_concat="$(dirname "{output.selected_fasta}")/contigs_info_concat.fasta"
        tmp_list="$(dirname "{output.selected_fasta}")/selected_contigs.list"
        test -d "{params.clusters_dir}"
        python3 workflows/smk_scripts/satminer_select_contigs_re2.py --clusters-dir "{params.clusters_dir}" --out-fasta "$tmp_concat" --out-list "$tmp_list" --coverage-fraction "{params.coverage_fraction}"
        python3 workflows/smk_scripts/satminer_extract_selected_contigs.py --fasta "$tmp_concat" --list "$tmp_list" --out "{output.selected_fasta}"
        rm -f "$tmp_concat" "$tmp_list" || true
        touch "{output.token}"
        """


rule concat_prev_comparative_iter_contigs:
    """Concatenate comparative selected contigs from iterations 1..(iter-1)."""
    input:
        prev_selected=lambda w: expand(
            "projects/{project}/comparative/{comparative_analysis}/tarean{iter}/iterative/selected_contigs.fasta",
            project=w.project,
            comparative_analysis=w.comparative_analysis,
            iter=list(range(1, int(w.iter))),
        )
    output:
        ref_fasta=temp("projects/{project}/comparative/{comparative_analysis}/tarean{iter}/iterative/reference_contigs.fasta"),
    log:
        f"{LOG_DIR}/concat_contigs_{{project}}_{{comparative_analysis}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/concat_contigs_{project}_{comparative_analysis}_iter{iter}.tsv"
    conda:
        "../envs/reportr.yaml"
    params:
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=2)
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=concat_prev_comparative_iter_contigs CATEGORY=wrapper project={wildcards.project} analysis={wildcards.comparative_analysis} iter={wildcards.iter}"
        python3 workflows/smk_scripts/concat_fastas.py --out "{output.ref_fasta}" --in {input.prev_selected}
        """

