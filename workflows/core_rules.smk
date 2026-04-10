# Core Rules - Stable Components
# HOLY WORKFLOW COMPLIANCE: All rules use the reportr environment exclusively
# This file contains well-tested, stable rules that should rarely change

# =============================================================================
# CORE READ PREPARATION RULES
# =============================================================================

rule clean_reads:
    """Clean reads using BBDuk - STABLE RULE"""
    input:
        R1 = lambda w: get_sample_reads(w.project, w.sample)[0],
        R2 = lambda w: get_sample_reads(w.project, w.sample)[1]
    output:
        R1 = "projects/{project}/samples/{sample}/filtered_reads/R1.fq",
        R2 = "projects/{project}/samples/{sample}/filtered_reads/R2.fq",
        cleaned = "projects/{project}/samples/{sample}/filtered_reads/cleaned.done",
        raw_fastqc_html1 = "projects/{project}/samples/{sample}/fastqc/RAW_R1_fastqc.html",
        raw_fastqc_html2 = "projects/{project}/samples/{sample}/fastqc/RAW_R2_fastqc.html",
        cleaned_fastqc_html1 = "projects/{project}/samples/{sample}/fastqc/CLEANED_R1_fastqc.html",
        cleaned_fastqc_html2 = "projects/{project}/samples/{sample}/fastqc/CLEANED_R2_fastqc.html"
    params:
        adapters = lambda w: get_param(w.project, "read_cleaning", "adapters", sample_id=w.sample),
        bbduk_params = lambda w: get_param(w.project, "read_cleaning", "bbduk_params", sample_id=w.sample),
        fastqc_threads = lambda w: get_param(w.project, "read_cleaning", "fastqc_threads", sample_id=w.sample),
        # IMPORTANT: keep tmp inside the sample directory (no repo-root tmp/).
        temp_dir = "projects/{project}/samples/{sample}/tmp",
        pythonhashseed = lambda w: get_param(w.project, "pythonhashseed"),
    log:
        f'{LOG_DIR}/clean_reads_{{project}}_{{sample}}.log'
    benchmark:
        "benchmarks/clean_reads_{project}_{sample}.tsv"
    threads:
        lambda w: get_param(w.project, "read_cleaning", "fastqc_threads", sample_id=w.sample)
    resources:
        bbduk_slots = 1,
        fastqc_slots = 1,
        mem_mb = 80000
    shell:
        """
        bash scripts/clean_reads.sh \
          {wildcards.project} \
          {wildcards.sample} \
          {params.pythonhashseed} \
          {input.R1} \
          {input.R2} \
          {output.R1} \
          {output.R2} \
          {output.cleaned} \
          {output.raw_fastqc_html1} \
          {output.raw_fastqc_html2} \
          {output.cleaned_fastqc_html1} \
          {output.cleaned_fastqc_html2} \
          {params.adapters} \
          "{params.bbduk_params}" \
          {params.fastqc_threads} \
          {params.temp_dir} \
          {log}
        """

rule prepare_reads:
    """Prepare reads for RepeatExplorer analysis - STABLE RULE"""
    input:
        R1 = "projects/{project}/samples/{sample}/filtered_reads/R1.fq",
        R2 = "projects/{project}/samples/{sample}/filtered_reads/R2.fq"
    output:
        prepared = "projects/{project}/samples/{sample}/tarean/prepared_forRE.fasta",
        tarean_token = "projects/{project}/samples/{sample}/tarean/PREPARATION_COMPLETE",
        sampled_fastqc_html1 = "projects/{project}/samples/{sample}/fastqc/SAMPLED_R1_fastqc.html",
        sampled_fastqc_html2 = "projects/{project}/samples/{sample}/fastqc/SAMPLED_R2_fastqc.html"
    params:
        project = "{project}",
        sample = "{sample}",
        sample_prefix = lambda w: get_sample_prefix(w.project, w.sample),
        pythonhashseed = lambda w: get_param(w.project, "pythonhashseed"),
        # IMPORTANT: keep tmp inside the sample directory (no repo-root tmp/).
        temp_dir = "projects/{project}/samples/{sample}/tmp",
        reads_per_assembly = lambda w: get_reads_per_assembly(w.project, w.sample),
        max_retries = 3,
        cleanup_after_prepare = lambda w: get_param(w.project, "cleanup_after_prepare"),
        fastqc_threads = lambda w: get_param(w.project, "read_cleaning", "fastqc_threads", sample_id=w.sample)
    log:
        f'{LOG_DIR}/prepare_reads_{{project}}_{{sample}}.log'
    benchmark:
        "benchmarks/prepare_reads_{project}_{sample}.tsv"
    resources:
        fastqc_slots = 1,
        bbduk_slots = 0,
        mem_mb = 8000
    shell:
        """
        bash scripts/prepare_reads.sh \
          {params.project} \
          {params.sample} \
          {params.sample_prefix} \
          {params.pythonhashseed} \
          {params.temp_dir} \
          {params.reads_per_assembly} \
          {params.max_retries} \
          {input.R1} \
          {input.R2} \
          {output.prepared} \
          {output.tarean_token} \
          {output.sampled_fastqc_html1} \
          {output.sampled_fastqc_html2} \
          {params.fastqc_threads} \
          {log}
        """

# =============================================================================
# CORE TAREAN RULES
# =============================================================================

rule run_tarean:
    """Run TAREAN repeat analysis - STABLE RULE"""
    input:
        prepared = "projects/{project}/samples/{sample}/tarean/prepared_forRE.fasta",
        tarean_token = "projects/{project}/samples/{sample}/tarean/PREPARATION_COMPLETE"
    output:
        tarean_done = "projects/{project}/samples/{sample}/tarean/tarean.done",
        tarean_log = f'{LOG_DIR}/tarean_{{project}}_{{sample}}.log',
        cluster_table = "projects/{project}/samples/{sample}/tarean/CLUSTER_TABLE.csv"
    params:
        project = "{project}",
        sample = "{sample}",
        reads_per_assembly = lambda w: get_reads_per_assembly(w.project, w.sample),
        sample_prefix = lambda w: get_sample_prefix(w.project, w.sample),
        pythonhashseed = lambda w: get_param(w.project, "pythonhashseed"),
        temp_dir = lambda w: get_param(w.project, "read_preparation", "temp_dir", sample_id=w.sample),
        assembly_min = lambda w: get_tarean_param_for_sample(w.project, w.sample, "assembly_min"),
        mincl_percent = lambda w: get_mincl_percent(w.project),
        min_lcov = lambda w: get_tarean_param_for_sample(w.project, w.sample, "min_lcov"),
        merge_threshold = lambda w: get_tarean_param_for_sample(w.project, w.sample, "merge_threshold"),
        r_value = lambda w: get_tarean_r_value(w.project),
        options = lambda w: get_tarean_param_for_sample(w.project, w.sample, "options"),
        paired = lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "paired")),
        automatic_filtering = lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "automatic_filtering")),
        tarean_mode = lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "tarean_mode")),
        keep_names = lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "keep_names")),
        cleanup = lambda w: bool(get_tarean_param_for_sample(w.project, w.sample, "cleanup")),
        domain_search = lambda w: get_tarean_param_for_sample(w.project, w.sample, "domain_search"),
        tarean_dir = "projects/{project}/samples/{sample}/tarean",
        rserv_port = lambda w: stable_rserv_port(42113, w.sample),  # Deterministic port per sample
        cleanup_after_prepare = lambda w: get_param(w.project, "cleanup_after_prepare")
    log:
        f'{LOG_DIR}/seqclust_{{project}}_{{sample}}.log'
    benchmark:
        "benchmarks/run_tarean_{project}_{sample}.tsv"
    threads:
        lambda w: get_seqclust_threads(
            w.project,
            get_tarean_param_for_sample(w.project, w.sample, "threads"),
        )
    resources:
        seqclust_slots = 1,
        # Keep Snakemake scheduling consistent with seqclust -r (kB).
        mem_mb = lambda w: max(1024, int(get_tarean_r_value(w.project) / 1024))
    conda:
        "../envs/reportr.yaml"
    script:
        "smk_scripts/run_tarean_step.py"

# =============================================================================
# CORE UTILITY RULES
# =============================================================================

rule validate_project:
    """Validate project structure and files - STABLE RULE"""
    output:
        touch("projects/{project}/validated.txt")
    params:
        project = "{project}"
    shell:
        """
        python3 -c "
from project_manager import ProjectManager
pm = ProjectManager()
if pm.validate_project('{params.project}'):
    print('Project {params.project} is valid')
else:
    print('Project {params.project} has issues')
"
        """ 