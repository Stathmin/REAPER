# Iterative RepeatExplorer2/TAREAN assemblies (optional; config-gated)
#
# This module implements tarean1..tareanN runs where only read preparation changes
# between iterations; seqclust runner (run_tarean_step.py) remains unchanged.


def _iterative_enabled(project_id: str) -> bool:
    return bool(get_param(project_id, "iterative_assembly", "enabled"))


def _iter_depth(project_id: str) -> int:
    return int(get_param(project_id, "iterative_assembly", "depth"))


def _solo_deconseq_db(project_id: str) -> str:
    return str(get_param(project_id, "deconseq", "db"))


def _solo_deconseq_threads(project_id: str, sample_id: str | None = None) -> int:
    return int(get_param(project_id, "deconseq", "threads", sample_id=sample_id))


def _solo_deconseq_removeifeitherbad_t(project_id: str, sample_id: str | None = None) -> str:
    v = bool(get_param(project_id, "deconseq", "removeifeitherbad", sample_id=sample_id))
    return "t" if v else "f"


def _solo_deconseq_minkmerhits(project_id: str, sample_id: str | None = None) -> int:
    return int(get_param(project_id, "deconseq", "minkmerhits", sample_id=sample_id))


def _deconseq_iterative_reference_coverage_fraction(project_id: str) -> float:
    return float(get_param(project_id, "deconseq", "iterative_reference_coverage_fraction"))


def _cleanup_after_prepare(project_id: str) -> bool:
    return bool(get_param(project_id, "cleanup_after_prepare"))


def _guard_iterative(project_id: str, iter_s: str, *, min_iter: int = 1) -> str:
    """Evaluate iterative gating in `params:` for shell-based rules."""
    if not _iterative_enabled(project_id):
        raise ValueError("Iterative assembly disabled for project (iterative_assembly.enabled=false)")
    it = int(iter_s)
    if it < 1 or it > _iter_depth(project_id):
        raise ValueError("Requested iteration outside configured depth")
    if it < min_iter:
        raise ValueError(f"Rule only valid for iter>={min_iter}")
    return "ok"


rule prepare_reads_iter:
    """Iteration-aware read prep: writes into tarean{iter}/."""
    input:
        r1=lambda w: (
            f"projects/{w.project}/samples/{w.sample}/filtered_reads/R1.fq"
            if int(w.iter) == 1
            else f"projects/{w.project}/samples/{w.sample}/tarean{w.iter}/iterative/unmapped_R1.fq"
        ),
        r2=lambda w: (
            f"projects/{w.project}/samples/{w.sample}/filtered_reads/R2.fq"
            if int(w.iter) == 1
            else f"projects/{w.project}/samples/{w.sample}/tarean{w.iter}/iterative/unmapped_R2.fq"
        ),
        unmapped_token=lambda w: (
            []
            if int(w.iter) == 1
            else f"projects/{w.project}/samples/{w.sample}/tarean{w.iter}/iterative/DECONSEQ_UNMAPPED_READY"
        ),
    output:
        prepared="projects/{project}/samples/{sample}/tarean{iter}/prepared_forRE.fasta",
        token="projects/{project}/samples/{sample}/tarean{iter}/PREPARATION_COMPLETE",
        sampled_fastqc_html1="projects/{project}/samples/{sample}/tarean{iter}/fastqc/SAMPLED_R1_fastqc.html",
        sampled_fastqc_html2="projects/{project}/samples/{sample}/tarean{iter}/fastqc/SAMPLED_R2_fastqc.html",
    params:
        project="{project}",
        sample="{sample}",
        iter="{iter}",
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=1),
        sample_prefix=lambda w: get_sample_prefix(w.project, w.sample),
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        # IMPORTANT: keep tmp inside the sample directory (no repo-root tmp/).
        temp_dir="projects/{project}/samples/{sample}/tmp",
        reads_per_assembly=lambda w: get_reads_per_assembly(w.project, w.sample),
        max_retries=3,
        fastqc_threads=lambda w: get_param(w.project, "read_cleaning", "fastqc_threads", sample_id=w.sample),
        cleanup_after_prepare=lambda w: _cleanup_after_prepare(w.project),
    log:
        f"{LOG_DIR}/prepare_reads_{{project}}_{{sample}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/prepare_reads_{project}_{sample}_iter{iter}.tsv"
    resources:
        fastqc_slots=1
    conda:
        "../envs/reportr.yaml"
    shell:
        r"""
        set -euo pipefail
        # Iteration 1 uses filtered reads directly; iterations >1 will be overridden by upstream filter rules.
        bash scripts/prepare_reads.sh \
          "{params.project}" \
          "{params.sample}" \
          "{params.sample_prefix}" \
          "{params.pythonhashseed}" \
          "{params.temp_dir}" \
          "{params.reads_per_assembly}" \
          "{params.max_retries}" \
          "{input.r1}" \
          "{input.r2}" \
          "{output.prepared}" \
          "{output.token}" \
          "{output.sampled_fastqc_html1}" \
          "{output.sampled_fastqc_html2}" \
          "{params.fastqc_threads}" \
          "{log}"

        if [[ "{params.cleanup_after_prepare}" == "True" ]]; then
          rm -rf "projects/{wildcards.project}/samples/{wildcards.sample}/tarean{wildcards.iter}/tmp" || true
          rm -rf "projects/{wildcards.project}/samples/{wildcards.sample}/tarean{wildcards.iter}/tarean/tmp" || true
          rm -rf "projects/{wildcards.project}/samples/{wildcards.sample}/tarean{wildcards.iter}/tarean/Rserv" || true
          rmdir "projects/{wildcards.project}/samples/{wildcards.sample}/tarean{wildcards.iter}/tmp" 2>/dev/null || true
        fi
        """


rule run_tarean_iter:
    """Iteration-aware seqclust run: uses tarean{iter}/ as run directory."""
    input:
        prepared="projects/{project}/samples/{sample}/tarean{iter}/prepared_forRE.fasta",
        token="projects/{project}/samples/{sample}/tarean{iter}/PREPARATION_COMPLETE",
    output:
        tarean_done="projects/{project}/samples/{sample}/tarean{iter}/tarean.done",
        tarean_log=f"{LOG_DIR}/tarean_{{project}}_{{sample}}_iter{{iter}}.log",
        cluster_table="projects/{project}/samples/{sample}/tarean{iter}/CLUSTER_TABLE.csv",
    params:
        project="{project}",
        sample="{sample}",
        reads_per_assembly=lambda w: get_reads_per_assembly(w.project, w.sample),
        sample_prefix=lambda w: get_sample_prefix(w.project, w.sample),
        pythonhashseed=lambda w: get_param(w.project, "pythonhashseed"),
        temp_dir=lambda w: get_param(w.project, "read_preparation", "temp_dir", sample_id=w.sample),
        assembly_min=lambda w: get_tarean_param_for_sample(w.project, w.sample, "assembly_min"),
        mincl_percent=lambda w: get_mincl_percent(w.project),
        min_lcov=lambda w: get_tarean_param_for_sample(w.project, w.sample, "min_lcov"),
        merge_threshold=lambda w: get_tarean_param_for_sample(w.project, w.sample, "merge_threshold"),
        r_value=lambda w: get_tarean_r_value(w.project),
        options=lambda w: get_tarean_param_for_sample(w.project, w.sample, "options"),
        paired=lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "paired")),
        automatic_filtering=lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "automatic_filtering")),
        tarean_mode=lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "tarean_mode")),
        keep_names=lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "keep_names")),
        cleanup=lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "cleanup")),
        domain_search=lambda w: get_tarean_param_for_sample(w.project, w.sample, "domain_search"),
        tarean_dir="projects/{project}/samples/{sample}/tarean{iter}",
        rserv_port=lambda w: stable_rserv_port(42113, f"{w.sample}_iter{w.iter}"),
        cleanup_after_prepare=lambda w: get_param(w.project, "cleanup_after_prepare"),
    log:
        f"{LOG_DIR}/seqclust_{{project}}_{{sample}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/run_tarean_{project}_{sample}_iter{iter}.tsv"
    threads:
        lambda w: get_seqclust_threads(
            w.project,
            get_tarean_param_for_sample(w.project, w.sample, "threads"),
        )
    resources:
        seqclust_slots=1,
        mem_mb=lambda w: max(1024, int(get_tarean_r_value(w.project) / 1024)),
    conda:
        "../envs/reportr.yaml"
    script:
        "smk_scripts/run_tarean_step.py"


rule select_contigs_iter_re2:
    """Select contigs from a specific iteration's clusters (RE2)."""
    input:
        done="projects/{project}/samples/{sample}/tarean{iter}/tarean.done",
    output:
        concat_fasta=temp("projects/{project}/samples/{sample}/tarean{iter}/iterative/contigs_info_concat.fasta"),
        contig_list=temp("projects/{project}/samples/{sample}/tarean{iter}/iterative/selected_contigs.list"),
        token="projects/{project}/samples/{sample}/tarean{iter}/iterative/CONTIG_LIST_READY",
    log:
        f"{LOG_DIR}/select_contigs_{{project}}_{{sample}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/select_contigs_{project}_{sample}_iter{iter}.tsv"
    conda:
        "../envs/reportr.yaml"
    params:
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=1),
        clusters_dir="projects/{project}/samples/{sample}/tarean{iter}/seqclust/clustering/clusters",
        coverage_fraction=lambda w: _deconseq_iterative_reference_coverage_fraction(w.project),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=select_contigs_iter_re2 CATEGORY=wrapper project={wildcards.project} sample={wildcards.sample} iter={wildcards.iter}"
        test -d "{params.clusters_dir}"
        python3 workflows/smk_scripts/satminer_select_contigs_re2.py \
          --clusters-dir "{params.clusters_dir}" \
          --out-fasta "{output.concat_fasta}" \
          --out-list "{output.contig_list}" \
          --coverage-fraction "{params.coverage_fraction}"
        touch "{output.token}"
        """


rule extract_selected_contigs_iter:
    """Extract selected contigs FASTA for a specific iteration."""
    input:
        concat_fasta="projects/{project}/samples/{sample}/tarean{iter}/iterative/contigs_info_concat.fasta",
        contig_list="projects/{project}/samples/{sample}/tarean{iter}/iterative/selected_contigs.list",
        token="projects/{project}/samples/{sample}/tarean{iter}/iterative/CONTIG_LIST_READY",
    output:
        selected_fasta="projects/{project}/samples/{sample}/tarean{iter}/iterative/selected_contigs.fasta",
        token="projects/{project}/samples/{sample}/tarean{iter}/iterative/SELECTED_CONTIGS_READY",
    log:
        f"{LOG_DIR}/extract_contigs_{{project}}_{{sample}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/extract_contigs_{project}_{sample}_iter{iter}.tsv"
    conda:
        "../envs/reportr.yaml"
    params:
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=1)
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=extract_selected_contigs_iter CATEGORY=wrapper project={wildcards.project} sample={wildcards.sample} iter={wildcards.iter}"
        python3 workflows/smk_scripts/satminer_extract_selected_contigs.py \
          --fasta "{input.concat_fasta}" \
          --list "{input.contig_list}" \
          --out "{output.selected_fasta}"
        touch "{output.token}"
        """


rule concat_prev_iter_contigs:
    """Concatenate selected contigs from iterations 1..(iter-1) for filtering reference."""
    input:
        prev_selected=lambda w: expand(
            "projects/{project}/samples/{sample}/tarean{iter}/iterative/selected_contigs.fasta",
            project=w.project,
            sample=w.sample,
            iter=list(range(1, int(w.iter))),
        )
    output:
        ref_fasta=temp("projects/{project}/samples/{sample}/tarean{iter}/iterative/reference_contigs.fasta"),
    log:
        f"{LOG_DIR}/concat_contigs_{{project}}_{{sample}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/concat_contigs_{project}_{sample}_iter{iter}.tsv"
    conda:
        "../envs/reportr.yaml"
    params:
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=2)
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=concat_prev_iter_contigs CATEGORY=wrapper project={wildcards.project} sample={wildcards.sample} iter={wildcards.iter}"
        python3 workflows/smk_scripts/concat_fastas.py --out "{output.ref_fasta}" --in {input.prev_selected}
        """


rule deconseq_filter_reads_iter:
    """Filter reads against previous-iteration contigs; keep unmapped reads (iter>=2)."""
    input:
        r1="projects/{project}/samples/{sample}/filtered_reads/R1.fq",
        r2="projects/{project}/samples/{sample}/filtered_reads/R2.fq",
        ref_fasta="projects/{project}/samples/{sample}/tarean{iter}/iterative/reference_contigs.fasta",
    output:
        unmapped_r1=temp("projects/{project}/samples/{sample}/tarean{iter}/iterative/unmapped_R1.fq"),
        unmapped_r2=temp("projects/{project}/samples/{sample}/tarean{iter}/iterative/unmapped_R2.fq"),
        token="projects/{project}/samples/{sample}/tarean{iter}/iterative/DECONSEQ_UNMAPPED_READY",
    log:
        f"{LOG_DIR}/deconseq_filter_{{project}}_{{sample}}_iter{{iter}}.log"
    benchmark:
        "benchmarks/deconseq_filter_{project}_{sample}_iter{iter}.tsv"
    threads:
        lambda w: _solo_deconseq_threads(w.project)
    resources:
        bbmap_slots=1
    conda:
        "../envs/reportr.yaml"
    params:
        __guard=lambda w: _guard_iterative(w.project, w.iter, min_iter=2),
        threads=lambda w: _solo_deconseq_threads(w.project, w.sample),
        kmer_k=lambda w: int(get_param(w.project, "deconseq", "kmer_k", sample_id=w.sample)),
        removeifeitherbad_t=lambda w: _solo_deconseq_removeifeitherbad_t(w.project, w.sample),
        minkmerhits=lambda w: _solo_deconseq_minkmerhits(w.project, w.sample),
        cleanup_after_prepare=lambda w: _cleanup_after_prepare(w.project),
    shell:
        r"""
        set -euo pipefail
        exec > "{log}" 2>&1
        echo "$(date): RULE=deconseq_filter_reads_iter CATEGORY=heavy project={wildcards.project} sample={wildcards.sample} iter={wildcards.iter} k={params.kmer_k} threads={threads}"
        mkdir -p "$(dirname "{output.unmapped_r1}")"

        # Per-job tmp under the iteration directory (delete on success; keep on failure).
        tmp_base="projects/{wildcards.project}/samples/{wildcards.sample}/tarean{wildcards.iter}/iterative/tmp/jobs/deconseq_filter_reads_iter"
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
            echo "$(date): DECONSEQ_TMPDIR_RETAINED=$job_tmp" >&2 || true
          fi
        }}
        trap _cleanup_job_tmp EXIT

        # Implementation note: Use k-mer filtering (bbduk) for speed on large FASTQs.
        # Pair policy: global.deconseq.removeifeitherbad (default true: drop pair if either mate matches).

        # Ensure pairs are synchronized.
        # Some upstream FASTQs have mate markers after whitespace; normalize by trimming
        # descriptions, then repair by (possibly identical) names.
        repair.sh \
          in="{input.r1}" in2="{input.r2}" \
          out1="$job_tmp/paired_R1.fq" out2="$job_tmp/paired_R2.fq" \
          outs="$job_tmp/singletons.fq" \
          repair=t \
          trimreaddescription=t \
          ain=t \
          overwrite=t \
          1>/dev/null

        # Run BBDuk on an interleaved stream to avoid paired-stream assertions.
        reformat.sh \
          in1="$job_tmp/paired_R1.fq" in2="$job_tmp/paired_R2.fq" \
          out="$job_tmp/paired_int.fq" \
          ow=t \
          1>/dev/null

        bbduk.sh \
          in="$job_tmp/paired_int.fq" \
          interleaved=t \
          ref="{input.ref_fasta}" \
          k="{params.kmer_k}" \
          threads="{threads}" \
          rcomp=t \
          minkmerhits="{params.minkmerhits}" \
          removeifeitherbad="{params.removeifeitherbad_t}" \
          out="$job_tmp/unmapped_int.fq" \
          overwrite=t \
          1>/dev/null

        reformat.sh \
          in="$job_tmp/unmapped_int.fq" \
          interleaved=t \
          out1="{output.unmapped_r1}" out2="{output.unmapped_r2}" \
          ow=t \
          1>/dev/null

        touch "{output.token}"
        """


