# Analysis Rules - Post-TAREAN and Comparative Analysis
# HOLY WORKFLOW COMPLIANCE: All rules use the reportr environment exclusively
# This file contains rules for post-processing and analysis

import json


def _comparative_samples(project_id: str, comparative_analysis: str) -> list[str]:
    proj = config.get("projects", {}).get(project_id, {}) or {}
    ca = (proj.get("comparative_analyses", {}) or {}).get(comparative_analysis, {}) or {}
    samples = ca.get("samples", []) or []
    return list(samples) if isinstance(samples, (list, tuple)) else []


def _comparative_prefix_length(project_id: str, comparative_analysis: str) -> int:
    base = int(get_tarean_param(project_id, "prefix_length"))
    # If analysis exists, allow analysis-level override for tarean_params.prefix_length.
    if comparative_analysis in (config.get("projects", {}).get(project_id, {}) or {}).get("comparative_analyses", {}):
        base = int(get_tarean_param_for_comparative(project_id, comparative_analysis, "prefix_length"))
    smpls = _comparative_samples(project_id, comparative_analysis)
    if not smpls:
        return base
    longest = max(len(str(config["projects"][project_id]["samples"][s].get("prefix") or str(s).upper())) for s in smpls)
    return max(base, int(longest))

# =============================================================================
# POST-TAREAN ANALYSIS RULES
# =============================================================================

rule post_tarean_blast:
    """Run post-TAREAN analysis for a single sample using the step-based pipeline."""
    input:
        cluster_table = "projects/{project}/samples/{sample}/tarean/CLUSTER_TABLE.csv",
        tarean_done = "projects/{project}/samples/{sample}/tarean/tarean.done"
    output:
        blast_results = "projects/{project}/samples/{sample}/post_tarean/{sample}_blast_results.csv",
        excel_report = "projects/{project}/samples/{sample}/post_tarean/{sample}_repeat_analysis.xlsx",
        word_report = "projects/{project}/samples/{sample}/post_tarean/{sample}_repeat_analysis.docx",
        csv_report = "projects/{project}/samples/{sample}/post_tarean/{sample}_repeat_analysis.csv",
        confidence_json = "projects/{project}/samples/{sample}/post_tarean/{sample}_confidence.json"
    params:
        project = "{project}",
        sample = "{sample}",
        enabled_steps = lambda w: config["projects"][w.project].get("post_tarean_params", {}).get(
            "analysis", {}
        ).get("enabled_steps", ["blast", "summary"])
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_blast_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_blast_{project}_{sample}.tsv"
    shell:
        """
        set -euo pipefail
        exec > "{log}" 2>&1
        mkdir -p projects/{params.project}/samples/{params.sample}/post_tarean/
        
        # Decide whether to run BLAST based on enabled steps from configuration.
        ENABLED_STEPS='{params.enabled_steps}'
        if echo "$ENABLED_STEPS" | grep -q "blast"; then
            REPORT_ONLY=""
        else
            REPORT_ONLY="--report-only"
        fi

        # Run the post-TAREAN pipeline as a module; internal orchestration
        # respects enabled_steps from configuration.
        python3 -m post_tarean.pipeline {params.sample} --project-id {params.project} --output-dir projects/{params.project}/samples/{params.sample}/post_tarean/ $REPORT_ONLY
        """

rule post_tarean_pipeline:
    """Run real post-TAREAN pipeline"""
    input:
        blast_results = "projects/{project}/samples/{sample}/post_tarean/{sample}_blast_results.csv",
        excel_report = "projects/{project}/samples/{sample}/post_tarean/{sample}_repeat_analysis.xlsx",
        word_report = "projects/{project}/samples/{sample}/post_tarean/{sample}_repeat_analysis.docx",
        csv_report = "projects/{project}/samples/{sample}/post_tarean/{sample}_repeat_analysis.csv"
    output:
        analysis_complete = "projects/{project}/samples/{sample}/post_tarean/analysis_complete.txt"
    params:
        project = "{project}",
        sample = "{sample}"
    conda:
        "../envs/reportr.yaml"
    threads: 1
    log:
        f"{LOG_DIR}/post_tarean_pipeline_{{project}}_{{sample}}.log"
    benchmark:
        "benchmarks/post_tarean_pipeline_{project}_{sample}.tsv"
    shell:
        """
        set -euo pipefail
        exec > "{log}" 2>&1
        # The pipeline script already generates all reports
        # Just ensure the completion marker exists
        touch {output.analysis_complete}
        """

# =============================================================================
# COMPARATIVE ANALYSIS RULES
# =============================================================================

rule prepare_comparative_reads:
    """Prepare reads for comparative analysis with genome size-based proportions"""
    input:
        # Comparative preparation must use cleaned reads only.
        # Do not depend on per-sample PREPARATION_COMPLETE (that is solo TAREAN prep).
        filtered_r1 = lambda w: expand(
            "projects/{project}/samples/{sample}/filtered_reads/R1.fq",
            project=w.project,
            sample=config["projects"][w.project]["comparative_analyses"][w.comparative_analysis]["samples"],
        ),
        filtered_r2 = lambda w: expand(
            "projects/{project}/samples/{sample}/filtered_reads/R2.fq",
            project=w.project,
            sample=config["projects"][w.project]["comparative_analyses"][w.comparative_analysis]["samples"],
        )
    output:
        combined_reads = "projects/{project}/comparative/{comparative_analysis}/comparative_reads.fasta",
        comparative_summary = "projects/{project}/comparative/{comparative_analysis}/comparative_summary.txt",
        # Comparative-specific sampled FASTQs (per sample) are written under this directory.
        sampled_fastq_dir = directory("projects/{project}/comparative/{comparative_analysis}/sampled_fastq")
    params:
        project = "{project}",
        comparative_analysis = "{comparative_analysis}",
        reads_per_assembly = lambda w: get_reads_per_assembly(w.project),
        # Configuration-driven proportions calculation
        proportions_json = lambda w: json.dumps({
            sample: config["projects"][w.project]["samples"][sample].get("genome_size", 1.0)
            for sample in config["projects"][w.project]["comparative_analyses"][w.comparative_analysis]["samples"]
        })
        ,
        # Deterministic RNG seed for comparative pair sampling.
        pythonhashseed = lambda w: get_param(w.project, "pythonhashseed")
    log:
        f'{LOG_DIR}/comparative_reads_{{project}}_{{comparative_analysis}}.log'
    benchmark:
        "benchmarks/prepare_comparative_reads_{project}_{comparative_analysis}.tsv"
    shell:
        """
        # Setup logging with enhanced error handling
        echo "$(date): Starting comparative reads preparation for {params.project}, analysis {params.comparative_analysis}" >> {log}
        
        # Create directories
        mkdir -p projects/{params.project}/comparative/{params.comparative_analysis}/
        
        # Log configuration-driven proportions
        echo "$(date): Using configuration-driven proportions: {params.proportions_json}" >> {log}
        
        # Prepare comparative reads with genome size-based proportions
        python3 workflows/smk_scripts/prepare_comparative_reads.py {params.project} {params.comparative_analysis} --seed {params.pythonhashseed} \
            2>> {log}
        
        if [ $? -eq 0 ]; then
            echo "$(date): Comparative reads preparation completed successfully" >> {log}
            echo "Configuration-driven proportions used successfully" >> {log}
            
            # Create summary file with detailed information
            echo "Comparative Analysis Summary" > {output.comparative_summary}
            echo "===========================" >> {output.comparative_summary}
            echo "Project: {params.project}" >> {output.comparative_summary}
            echo "Analysis: {params.comparative_analysis}" >> {output.comparative_summary}
            echo "Date: $(date)" >> {output.comparative_summary}
            echo "" >> {output.comparative_summary}
            echo "Samples included:" >> {output.comparative_summary}
            echo "Total reads: $(grep -c '^>' {output.combined_reads})" >> {output.comparative_summary}
            echo "File: {output.combined_reads}" >> {output.comparative_summary}
            echo "" >> {output.comparative_summary}
            echo "Configuration-driven proportions:" >> {output.comparative_summary}
            echo "{params.proportions_json}" >> {output.comparative_summary}

            # Cleanup: comparative prep may create large temporary artifacts.
            # Always clear tmp after a successful run.
            rm -rf "projects/{params.project}/comparative/{params.comparative_analysis}/tmp" || true
        else
            echo "$(date): Error in comparative reads preparation" >> {log}
            echo "ERROR: Comparative analysis failed. Check configuration and input files." >> {log}
            exit 1
        fi
        """


rule comparative_fastqc_sampled:
    """FastQC on comparative-specific sampled reads (per sample)."""
    input:
        sampled_fastq_dir = "projects/{project}/comparative/{comparative_analysis}/sampled_fastq"
    output:
        html1 = "projects/{project}/comparative/{comparative_analysis}/fastqc/{sample}/SAMPLED_R1_fastqc.html",
        html2 = "projects/{project}/comparative/{comparative_analysis}/fastqc/{sample}/SAMPLED_R2_fastqc.html"
    params:
        threads = lambda w: get_param(w.project, "read_cleaning", "fastqc_threads", comparative_context=True),
        r1 = lambda w: f"projects/{w.project}/comparative/{w.comparative_analysis}/sampled_fastq/{w.sample}/R1.fq",
        r2 = lambda w: f"projects/{w.project}/comparative/{w.comparative_analysis}/sampled_fastq/{w.sample}/R2.fq",
    log:
        f'{LOG_DIR}/fastqc_comparative_sampled_{{project}}_{{comparative_analysis}}_{{sample}}.log'
    shell:
        r"""
        set -euo pipefail
        outdir="projects/{wildcards.project}/comparative/{wildcards.comparative_analysis}/fastqc/{wildcards.sample}"
        mkdir -p "$outdir"
        fastqc -t {params.threads} "{params.r1}" "{params.r2}" -o "$outdir" --quiet > "{log}" 2>&1

        # FastQC names outputs after the input basenames (R1_fastqc.html, R2_fastqc.html, ...).
        mv -f "$outdir/R1_fastqc.html" "{output.html1}"
        mv -f "$outdir/R2_fastqc.html" "{output.html2}"
        if [ -f "$outdir/R1_fastqc.zip" ]; then mv -f "$outdir/R1_fastqc.zip" "$outdir/SAMPLED_R1_fastqc.zip"; fi
        if [ -f "$outdir/R2_fastqc.zip" ]; then mv -f "$outdir/R2_fastqc.zip" "$outdir/SAMPLED_R2_fastqc.zip"; fi
        """


rule comparative_ready:
    """Comparative preparation completion token (reads + per-sample sampled FastQC)."""
    input:
        combined_reads = "projects/{project}/comparative/{comparative_analysis}/comparative_reads.fasta",
        comparative_summary = "projects/{project}/comparative/{comparative_analysis}/comparative_summary.txt",
        sampled_fastqc = lambda w: expand(
            "projects/{project}/comparative/{comparative_analysis}/fastqc/{sample}/SAMPLED_R1_fastqc.html",
            project=w.project,
            comparative_analysis=w.comparative_analysis,
            sample=_comparative_samples(w.project, w.comparative_analysis),
        )
    output:
        comparative_token = "projects/{project}/comparative/{comparative_analysis}/COMPARATIVE_READY"
    shell:
        """
        touch {output.comparative_token}
        """

rule comparative_tarean_analysis:
    """Run TAREAN analysis on comparative reads"""
    input:
        comparative_reads = "projects/{project}/comparative/{comparative_analysis}/comparative_reads.fasta",
        comparative_token = "projects/{project}/comparative/{comparative_analysis}/COMPARATIVE_READY"
    output:
        comparative_complete = "projects/{project}/comparative/{comparative_analysis}/COMPARATIVE_TAREAN_COMPLETE",
        tarean_log = f'{LOG_DIR}/tarean_comparative_{{project}}_{{comparative_analysis}}.log',
        cluster_table = "projects/{project}/comparative/{comparative_analysis}/tarean/CLUSTER_TABLE.csv"
    params:
        project = "{project}",
        sample = "{comparative_analysis}",
        comparative_analysis = "{comparative_analysis}",
        reads_per_assembly = lambda w: get_reads_per_assembly(
            w.project,
            comparative_analysis=w.comparative_analysis,
            comparative_context=True,
        ),
        proportions_json = lambda w: json.dumps(
            {
                sample: config["projects"][w.project]["samples"][sample].get("genome_size", 1.0)
                for sample in _comparative_samples(w.project, w.comparative_analysis)
            }
        ),
        assembly_min = lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "assembly_min"),
        mincl_percent = lambda w: get_mincl_percent(w.project),
        min_lcov = lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "min_lcov"),
        merge_threshold = lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "merge_threshold"),
        r_value = lambda w: get_tarean_r_value(w.project),
        # Propagate parity fields used by run_tarean_step.py for determinism/logging.
        pythonhashseed = lambda w: get_param(w.project, "pythonhashseed"),
        temp_dir = lambda w: get_param(w.project, "read_preparation", "temp_dir", comparative_context=True),
        threads = lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "threads"),
        options = lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "options"),
        paired = lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "paired")),
        automatic_filtering = lambda w: bool(
            get_tarean_param_for_comparative(w.project, w.comparative_analysis, "automatic_filtering")
        ),
        tarean_mode = lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "tarean_mode")),
        keep_names = lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "keep_names")),
        cleanup = lambda w: bool(get_tarean_param_for_comparative(w.project, w.comparative_analysis, "cleanup")),
        domain_search = lambda w: get_tarean_param_for_comparative(w.project, w.comparative_analysis, "domain_search"),
        # seqclust comparative mode uses `-P` (prefix length) to decide how to
        # compare read identifiers. `-P` must be at least as long as the
        # longest per-sample header prefix created by `scripts/prepare_reads.sh`.
        #
        # If `-P` is too small (e.g. `KA1` vs `KA12`), both samples share the
        # same truncated prefix and seqclust cannot reliably separate them.
        prefix_length = lambda w: _comparative_prefix_length(w.project, w.comparative_analysis),
        tarean_dir = "projects/{project}/comparative/{comparative_analysis}/tarean",
        rserv_port = lambda w: stable_rserv_port(43113, w.comparative_analysis),
        cleanup_after_prepare = lambda w: get_param(w.project, "cleanup_after_prepare"),
    log:
        f'{LOG_DIR}/seqclust_{{project}}_{{comparative_analysis}}.log'
    benchmark:
        "benchmarks/run_tarean_{project}_{comparative_analysis}.tsv"
    threads:
        lambda w: get_seqclust_threads(
            w.project,
            get_tarean_param(w.project, "threads"),
        )
    resources:
        seqclust_slots = 1,
        mem_mb = lambda w: max(1024, int(get_tarean_r_value(w.project) / 1024))
    conda:
        "../envs/reportr.yaml"
    script:
        "smk_scripts/run_tarean_step.py"

rule comparative_post_tarean_analysis:
    """Run post-TAREAN analysis on comparative results"""
    input:
        cluster_table = "projects/{project}/comparative/{comparative_analysis}/tarean/CLUSTER_TABLE.csv",
        comparative_complete = "projects/{project}/comparative/{comparative_analysis}/COMPARATIVE_TAREAN_COMPLETE",
        ncbi_freshness_ok = lambda w: (
            f"projects/{w.project}/metadata/ncbi_freshness_ok.txt"
            if config.get("global", {}).get("ncbi_gathering", {}).get("enabled", False)
            else []
        )
    output:
        blast_results = "projects/{project}/comparative/{comparative_analysis}/post_tarean/{comparative_analysis}_blast_results.csv",
        excel_report = "projects/{project}/comparative/{comparative_analysis}/post_tarean/{comparative_analysis}_repeat_analysis.xlsx",
        word_report = "projects/{project}/comparative/{comparative_analysis}/post_tarean/{comparative_analysis}_repeat_analysis.docx",
        csv_report = "projects/{project}/comparative/{comparative_analysis}/post_tarean/{comparative_analysis}_repeat_analysis.csv",
        analysis_complete = "projects/{project}/comparative/{comparative_analysis}/post_tarean/analysis_complete.txt"
    params:
        project = "{project}",
        comparative_analysis = "{comparative_analysis}"
    shell:
        """
        mkdir -p projects/{params.project}/comparative/{params.comparative_analysis}/post_tarean/
        
        # Run post-TAREAN pipeline for comparative analysis
        python3 post_tarean/pipeline.py {params.comparative_analysis} --project-id {params.project} --output-dir projects/{params.project}/comparative/{params.comparative_analysis}/post_tarean/
        
        # Create completion marker
        touch {output.analysis_complete}
        """

# =============================================================================
# QUALITY CONTROL RULES
# =============================================================================

# FastQC is now integrated into the clean_reads rule in core_rules.smk 