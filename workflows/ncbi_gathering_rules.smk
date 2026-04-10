# NCBI Data Gathering Rules for RepOrtR
# This file contains rules for gathering repeat sequences from NCBI
# based on taxonomic groups for annotation purposes

from pathlib import Path

# NCBI gathering configuration
ncbi_config = (config.get("global", {}) or {}).get("ncbi_gathering", {}) or {}
LOG_DIR = (config.get("global", {}) or {}).get("log_dir", "logs")
default_taxid = str(ncbi_config.get("default_taxonomy_id", "147389"))


def get_project_taxonomy_ids(project_name):
    """Get taxonomy IDs configured for a specific project.

    Resolution order:
      1) projects.<project_name>.taxonomy_ids (list of strings)
      2) projects.<project_name>.taxonomy_id (single string)
      3) global.ncbi_gathering.default_taxonomy_id
    """
    projects_cfg = config.get("projects", {})
    proj_cfg = projects_cfg.get(project_name, {})
    taxonomy_ids = proj_cfg.get("taxonomy_ids")
    if taxonomy_ids:
        return taxonomy_ids
    taxonomy_id = proj_cfg.get("taxonomy_id")
    if taxonomy_id:
        return [taxonomy_id]
    return [default_taxid]

# Rule to gather NCBI sequences for a project
rule gather_ncbi_sequences:
    """
    Gather repeat sequences from NCBI for a specific taxonomic group
    """
    output:
        fasta = "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta",
        blast_db = "projects/{project}/blast_db/ncbi_repeats_{taxid}.nhr",
        metadata = "projects/{project}/metadata/ncbi_metadata_{taxid}.csv"
    params:
        taxid = lambda wildcards: wildcards.taxid,
        project_id = lambda wildcards: wildcards.project,
        email = ncbi_config.get("email", "example@example.com")
    log:
        f'{LOG_DIR}/ncbi_gathering_{{project}}_{{taxid}}.log'
    threads: 2
    script:
        "smk_scripts/gather_ncbi_sequences.py"

# Rule to (re)build the on-disk BLAST database exactly where post_tarean
# expects it (next to the project-specific NCBI FASTA). This keeps the DB
# prefix aligned with post_tarean `post_tarean_params.databases.local_db_path`
# and `database_names.local` for projects like Triticeae.
rule rebuild_ncbi_blast_db:
    """
    Rebuild a BLAST database from the project-specific NCBI repeats FASTA.

    This is a lightweight wrapper around `makeblastdb` that ensures:
      - input  : projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta
      - db base: projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta

    so that BLAST can be invoked with:
      -db projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta
    and the index files (*.nhr, *.nin, *.nsq, etc.) live alongside the FASTA.
    """
    input:
        fasta = "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta",
    output:
        nhr = "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta.nhr",
    log:
        f'{LOG_DIR}/rebuild_ncbi_blast_db_{{project}}_{{taxid}}.log'
    threads: 2
    shell:
        "mkdir -p $(dirname {output.nhr}) && "
        "makeblastdb -in {input.fasta} -dbtype nucl "
        "-parse_seqids "
        "-out {input.fasta} "
        "-title NCBI_repeats_{wildcards.taxid}_{wildcards.project} "
        "> {log} 2>&1"

# Rule to update BLAST database with NCBI sequences
rule update_blast_database:
    """
    Update project BLAST database with NCBI sequences
    """
    input:
        fasta = "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta",
        existing_db = "projects/{project}/blast_db/existing_sequences.nhr"
    output:
        updated_db = "projects/{project}/blast_db/updated_sequences_{taxid}.nhr"
    log:
        f'{LOG_DIR}/blast_db_update_{{project}}_{{taxid}}.log'
    threads: 2
    script:
        "smk_scripts/update_blast_database.py"

# Rule to validate gathered sequences
rule validate_ncbi_sequences:
    """
    Validate gathered NCBI sequences for quality and completeness
    """
    input:
        fasta = "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta",
        metadata = "projects/{project}/metadata/ncbi_metadata_{taxid}.csv"
    output:
        validation_report = "projects/{project}/metadata/validation_report_{taxid}.txt"
    log:
        f'{LOG_DIR}/validation_{{project}}_{{taxid}}.log'
    threads: 1
    script:
        "smk_scripts/validate_ncbi_sequences.py"

# Rule to create taxonomic group summary
rule create_taxonomy_summary:
    """
    Create summary of gathered sequences by taxonomic group
    """
    input:
        metadata_files = lambda wildcards: expand(
            "projects/{project}/metadata/ncbi_metadata_{taxid}.csv",
            project=wildcards.project,
            taxid=get_project_taxonomy_ids(wildcards.project),
        )
    output:
        summary = "projects/{project}/metadata/taxonomy_summary.csv"
    log:
        f'{LOG_DIR}/taxonomy_summary_{{project}}.log'
    threads: 1
    script:
        "smk_scripts/create_taxonomy_summary.py"

# Rule to prepare sequences for ML training
rule prepare_ml_training_data:
    """
    Prepare gathered sequences for machine learning training
    """
    input:
        fasta_files = lambda wildcards: expand(
            "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta",
            project=wildcards.project,
            taxid=get_project_taxonomy_ids(wildcards.project),
        )
    output:
        training_data = "projects/{project}/ml/training_data.fasta",
        labels = "projects/{project}/ml/sequence_labels.csv"
    log:
        f'{LOG_DIR}/prepare_ml_data_{{project}}.log'
    threads: 2
    script:
        "smk_scripts/prepare_ml_training_data.py"

rule check_ncbi_freshness:
    """
    Check if NCBI sequences need to be updated based on age
    """
    input:
        metadata = "projects/{project}/metadata/ncbi_metadata_{taxid}.csv"
    output:
        freshness_report = "projects/{project}/metadata/freshness_report_{taxid}.txt"
    params:
        max_age_days = ncbi_config.get("max_age_days")  # Update if older than configured cadence
    log:
        f'{LOG_DIR}/freshness_check_{{project}}_{{taxid}}.log'
    threads: 1
    script:
        "smk_scripts/check_ncbi_freshness.py"


rule ensure_ncbi_fresh_for_project:
    """
    Ensure that NCBI-derived sequences for a project are fresh enough.
    This rule aggregates freshness checks for all configured taxonomy IDs;
    gathering itself is only re-triggered when explicitly requested.
    """
    input:
        lambda wildcards: expand(
            "projects/{project}/metadata/freshness_report_{taxid}.txt",
            project=wildcards.project,
            taxid=get_project_taxonomy_ids(wildcards.project),
        )
    output:
        "projects/{project}/metadata/ncbi_freshness_ok.txt"
    log:
        f'{LOG_DIR}/ncbi_freshness_ok_{{project}}.log'
    threads: 1
    shell:
        "echo 'NCBI data freshness checked for project {wildcards.project}' > {output}"

# Aggregated rule for complete NCBI data gathering workflow
rule gather_all_ncbi_data:
    """
    Complete NCBI data gathering workflow for a project
    """
    input:
        # Gather for all taxonomy IDs configured for this specific project.
        fasta_files = lambda wildcards: expand(
            "projects/{project}/blast_db/ncbi_repeats_{taxid}.fasta",
            project=wildcards.project,
            taxid=get_project_taxonomy_ids(wildcards.project),
        ),
        validation_reports = lambda wildcards: expand(
            "projects/{project}/metadata/validation_report_{taxid}.txt",
            project=wildcards.project,
            taxid=get_project_taxonomy_ids(wildcards.project),
        ),
        taxonomy_summary = "projects/{project}/metadata/taxonomy_summary.csv",
        ml_data = "projects/{project}/ml/training_data.fasta"
    output:
        completion_marker = "projects/{project}/ncbi_gathering_complete.txt"
    log:
        f'{LOG_DIR}/complete_ncbi_gathering_{{project}}.log'
    threads: 1
    script:
        "smk_scripts/gather_all_ncbi_data.py"

# Dynamic rule generation for multiple taxonomy IDs
def generate_taxonomy_rules():
    """Generate rules for all configured taxonomy IDs"""
    projects = config.get('projects', {})
    for project_name, project_config in projects.items():
        taxonomy_ids = project_config.get('taxonomy_ids', [default_taxid])
        for taxid in taxonomy_ids:
            # Generate specific rules for each taxonomy ID
            rule_name = f"gather_ncbi_{project_name}_{taxid}"
            # This would create specific rules for each project/taxonomy combination
            pass

# Export the rules
__all__ = [
    'gather_ncbi_sequences',
    'update_blast_database', 
    'validate_ncbi_sequences',
    'create_taxonomy_summary',
    'prepare_ml_training_data',
    'check_ncbi_freshness',
    'gather_all_ncbi_data'
]

