# Report Rules - Final Aggregation and Reporting
# HOLY WORKFLOW COMPLIANCE: All rules use the reportr environment exclusively
# This file contains rules for generating final reports and summaries

# =============================================================================
# AGGREGATE RULES
# =============================================================================

rule project_summary:
    """Generate project summary report"""
    input:
        post_tarean_reports = lambda w: expand("projects/{project}/samples/{sample}/post_tarean/analysis_complete.txt",
                                              project=w.project,
                                              sample=config["projects"][w.project]["samples"]),
        comparative_reports = lambda w: expand(
            "projects/{project}/comparative/{analysis}/post_tarean/analysis_complete.txt",
            project=w.project,
            analysis=config["projects"][w.project].get("comparative_analyses", {}),
        )
    output:
        summary = "projects/{project}/reports/project_summary.html"
    params:
        project = "{project}",
        sample_ids = lambda w: ", ".join(sorted((config["projects"][w.project].get("samples", {}) or {}).keys())),
        comparative_ids = lambda w: ", ".join(sorted((config["projects"][w.project].get("comparative_analyses", {}) or {}).keys())),
    shell:
        """
        mkdir -p projects/{params.project}/reports/
        SAMPLES="{params.sample_ids}"
        COMPS="{params.comparative_ids}"
        if [ -z "$SAMPLES" ]; then SAMPLES="(none)"; fi
        if [ -z "$COMPS" ]; then COMPS="(none)"; fi
        echo "<html><body><h1>Project Summary: {params.project}</h1><p>Analysis complete.</p><p>Post-TAREAN analysis completed for all samples.</p><p>Comparative analyses configured: $COMPS.</p><p>Samples: $SAMPLES.</p></body></html>" > {output.summary}
        """

# =============================================================================
# ALL RULES
# =============================================================================

rule all:
    """Main rule - run all analyses for a project"""
    input:
        # Dynamic sample completion tracking
        sample_completions = lambda w: expand("projects/{project}/samples/{sample}/tarean/tarean.done",
                                             project=w.project,
                                             sample=config["projects"][w.project]["samples"]),
        # Dynamic comparative analysis tracking
        comparative_completions = lambda w: expand("projects/{project}/comparative/{analysis}/comparative_reads.fasta",
                                                  project=w.project,
                                                  analysis=config["projects"][w.project].get("comparative_analyses", {}))
    output:
        touch("projects/{project,[^/]+}/analysis_complete.txt")
    params:
        project = "{project}",
        sample_ids = lambda w: ", ".join(sorted((config["projects"][w.project].get("samples", {}) or {}).keys())),
        comparative_ids = lambda w: ", ".join(sorted((config["projects"][w.project].get("comparative_analyses", {}) or {}).keys())),
    shell:
        """
        echo "All {params.project} analyses completed successfully!"
        echo "Individual samples: {params.sample_ids}"
        echo "Comparative analyses: {params.comparative_ids}"
        touch {output}
        """

# =============================================================================
# UTILITY RULES
# =============================================================================

rule clean_project:
    """Clean all outputs for a specific project"""
    output:
        touch("projects/{project}/cleaned.txt")
    params:
        project = "{project}"
    shell:
        """
        rm -rf projects/{params.project}/samples/*/cleaned_reads/
        rm -rf projects/{params.project}/samples/*/assembly/
        rm -rf projects/{params.project}/samples/*/tarean/
        rm -rf projects/{params.project}/samples/*/fastqc/
        rm -rf projects/{params.project}/comparative/
        rm -rf projects/{params.project}/reports/
        echo "Cleaned project {params.project}"
        """ 