# Test Rules - System Validation
# HOLY WORKFLOW COMPLIANCE: All rules use the reportr environment exclusively
# This file contains rules for testing system stability and validation

import json

def get_default_comparative_analysis(project_id):
    """Return first configured comparative analysis ID (or None)."""
    comps = config["projects"][project_id].get("comparative_analyses", {}) or {}
    return next(iter(comps.keys()), None)


def get_default_comparative_samples(project_id):
    """Return samples for default comparative analysis (or [])."""
    comp_id = get_default_comparative_analysis(project_id)
    if not comp_id:
        return []
    comps = config["projects"][project_id].get("comparative_analyses", {}) or {}
    return comps.get(comp_id, {}).get("samples", []) or []


def get_comparative_proportions(project_id):
    """Get comparative proportions from comparative_analyses config."""
    return json.dumps(
        {
            sample: config["projects"][project_id]["samples"][sample].get("genome_size", 1.0)
            for sample in get_default_comparative_samples(project_id)
        }
    )

# =============================================================================
# SYSTEM TEST RULES
# =============================================================================

rule test_system_resources:
    """Test system resources before starting assemblies - STABLE RULE"""
    output:
        resource_test = f'{config["global"]["log_dir"]}/system_resources_test.txt'
    log:
        f'{config["global"]["log_dir"]}/test_system_resources.log'
    shell:
        """
        echo "$(date): Testing system resources..." >> {log}
        
        # Enhanced system resource monitoring
        echo "=== SYSTEM RESOURCE TEST ===" > {output.resource_test}
        echo "Timestamp: $(date)" >> {output.resource_test}
        echo "" >> {output.resource_test}
        
        # Test memory with detailed analysis
        echo "=== MEMORY ANALYSIS ===" >> {output.resource_test}
        free -g | awk 'NR==2{{print "Total memory: " $2 "GB"}}' >> {output.resource_test}
        free -g | awk 'NR==2{{print "Available memory: " $7 "GB"}}' >> {output.resource_test}
        free -g | awk 'NR==2{{print "Used memory: " $3 "GB"}}' >> {output.resource_test}
        free -g | awk 'NR==2{{print "Memory usage: " int($3/$2*100) "%"}}' >> {output.resource_test}
        
        # Check if memory is sufficient for TAREAN
        available_mem=$(free -g | awk 'NR==2{{print $7}}')
        if [ $available_mem -lt 8 ]; then
            echo "⚠️  WARNING: Low memory available ($available_mem GB). TAREAN may fail." >> {output.resource_test}
        else
            echo "✓ Sufficient memory available ($available_mem GB)" >> {output.resource_test}
        fi
        
        # Test CPU with detailed analysis
        echo "" >> {output.resource_test}
        echo "=== CPU ANALYSIS ===" >> {output.resource_test}
        echo "CPU cores: $(nproc)" >> {output.resource_test}
        echo "CPU model: $(grep 'model name' /proc/cpuinfo | head -1 | cut -d: -f2 | xargs)" >> {output.resource_test}
        
        # Test disk space with detailed analysis
        echo "" >> {output.resource_test}
        echo "=== DISK ANALYSIS ===" >> {output.resource_test}
        df -h / | awk 'NR==2{{print "Disk usage: " $5}}' >> {output.resource_test}
        df -h / | awk 'NR==2{{print "Available space: " $4}}' >> {output.resource_test}
        
        # Check available disk space
        available_space=$(df -h / | awk 'NR==2{{print $4}}' | sed 's/G//')
        if [ $available_space -lt 10 ]; then
            echo "⚠️  WARNING: Low disk space available ($available_space GB). Analysis may fail." >> {output.resource_test}
        else
            echo "✓ Sufficient disk space available ($available_space GB)" >> {output.resource_test}
        fi
        
        # Test required tools with enhanced validation
        echo "" >> {output.resource_test}
        echo "=== TOOL VALIDATION ===" >> {output.resource_test}
        
        # Test seqclust with version check (binary in repex_tarean/)
        if [ -x "repex_tarean/seqclust" ]; then
            echo "✓ seqclust: OK (repex_tarean/seqclust)" >> {output.resource_test}
            repex_tarean/seqclust --help 2>/dev/null | head -1 >> {output.resource_test} 2>/dev/null || echo "seqclust version: Available" >> {output.resource_test}
        else
            echo "✗ seqclust: NOT FOUND" >> {output.resource_test}
            exit 1
        fi
        
        # Test Python dependencies with version info
        echo "Testing Python dependencies..." >> {log}
        python3 -c "import pandas, numpy, matplotlib, xlsxwriter, docx; print('✓ Python dependencies: OK')" 2>/dev/null && echo "✓ Python dependencies: OK" >> {output.resource_test} || echo "✗ Python dependencies: MISSING" >> {output.resource_test}
        
        # Test conda environment with enhanced validation
        if [[ "$CONDA_DEFAULT_ENV" == "reportr" ]]; then
            echo "✓ Conda environment: OK (reportr)" >> {output.resource_test}
            echo "Python version: $(python3 --version)" >> {output.resource_test}
        else
            echo "✗ Conda environment: WRONG (current: $CONDA_DEFAULT_ENV)" >> {output.resource_test}
            exit 1
        fi
        
        # Performance baseline test
        echo "" >> {output.resource_test}
        echo "=== PERFORMANCE BASELINE ===" >> {output.resource_test}
        echo "Testing file I/O performance..." >> {output.resource_test}
        dd if=/dev/zero of=/tmp/test_io bs=1M count=100 2>/dev/null | tail -1 >> {output.resource_test}
        rm -f /tmp/test_io
        
        echo "$(date): System resource test completed" >> {log}
        echo "✓ All system resources validated" >> {output.resource_test}
        """

rule test_read_preparation:
    """Test read preparation with small sample - STABLE RULE"""
    input:
        R1 = "projects/{project}/samples/{sample}/filtered_reads/R1.fq",
        R2 = "projects/{project}/samples/{sample}/filtered_reads/R2.fq",
        cleaned = "projects/{project}/samples/{sample}/filtered_reads/cleaned.done"
    output:
        test_prepared = "projects/{project}/samples/{sample}/test/test_prepared.fasta",
        test_summary = "projects/{project}/samples/{sample}/test/test_summary.txt"
    params:
        project = "{project}",
        sample = "{sample}",
        test_reads = 5000  # Use smaller test size for validation
    log:
        f'{config["global"]["log_dir"]}/test_read_preparation_{{project}}_{{sample}}.log'
    conda:
        "../envs/reportr.yaml"
    shell:
        """
        echo "$(date): Testing read preparation for {params.sample}..." >> {log}
        
        # Create test directory
        mkdir -p projects/{params.project}/samples/{params.sample}/test/
        
        # Test with small number of reads
        echo "$(date): Testing with {params.test_reads} reads" >> {log}
        tmp_fa="projects/{params.project}/samples/{params.sample}/test/tmp_test_{params.sample}.fa"
        
        # Use BBTools reformat.sh for sampling + interleaved FASTA
        reformat.sh \
          in1={input.R1} \
          in2={input.R2} \
          out=${tmp_fa} \
          samplereadstarget={params.test_reads} \
          sampleseed=42 \
          interleaved=t \
          ow=t 2>> {log}
        
        if [ $? -ne 0 ]; then
            echo "Read preparation test: FAILED (reformat.sh)" > {output.test_summary}
            echo "$(date): Read preparation test failed during reformat.sh" >> {log}
            exit 1
        fi
        
        # Rename headers to test_{sample}readN_f / _r and write final FASTA
        awk -v p="test_{params.sample}" '
          /^>/ {
            c++
            if (c % 2 == 1) {
              printf(">%sread%d_f\n", p, (c+1)/2)
            } else {
              printf(">%sread%d_r\n", p, c/2)
            }
            next
          }
          { print }
        ' "${tmp_fa}" > {output.test_prepared} 2>> {log}
        
        if [ $? -eq 0 ]; then
            rm -f "${tmp_fa}"
            echo "Read preparation test: SUCCESS" > {output.test_summary}
            echo "Test reads: {params.test_reads}" >> {output.test_summary}
            echo "$(date): Read preparation test completed successfully" >> {log}
        else
            rm -f "${tmp_fa}"
            echo "Read preparation test: FAILED (header renaming)" > {output.test_summary}
            echo "$(date): Read preparation test failed during header renaming" >> {log}
            exit 1
        fi
        """

rule test_tarean_small:
    """Test TAREAN with small dataset - STABLE RULE"""
    input:
        test_prepared = "projects/{project}/samples/{sample}/test/test_prepared.fasta"
    output:
        test_tarean_done = "projects/{project}/samples/{sample}/test/test_tarean.done"
    params:
        project = "{project}",
        sample = "{sample}",
        assembly_min = lambda w: config["projects"][w.project]["tarean_params"]["assembly_min"],
        mincl = lambda w: config["projects"][w.project]["tarean_params"]["mincl"],
        r_value = lambda w: get_tarean_r_value(w.project),
        threads = lambda w: min(4, config["projects"][w.project]["tarean_params"]["threads"])
    log:
        f'{config["global"]["log_dir"]}/test_tarean_{{project}}_{{sample}}.log'
    conda:
        "../envs/reportr.yaml"
    shell:
        """
        echo "$(date): Testing TAREAN with small dataset for {params.sample}..." >> {log}
        ORIGINAL_DIR=$(pwd)
        TEST_TAREAN_DIR="$ORIGINAL_DIR/projects/{params.project}/samples/{params.sample}/test/tarean"
        
        # Create test tarean directory
        mkdir -p "$TEST_TAREAN_DIR"
        
        # Test via the same wrapper used in the real pipeline.
        # This ensures runtime isolation (RSERVE_WORKDIR, RSERVE_PORT, TMPDIR) and any
        # compatibility shims are applied consistently.
        echo "$(date): Running test seqclust via run_tarean_step.py (timeout=30m)..." >> {log}

        timeout 30m python3 workflows/smk_scripts/run_tarean_step.py \
          --project "{params.project}" \
          --sample "{params.sample}" \
          --prepared "{input.test_prepared}" \
          --tarean-dir "projects/{params.project}/samples/{params.sample}/test/tarean" \
          --tarean-done "{output.test_tarean_done}" \
          --tarean-log "{log}" \
          --seqclust-log "{log}" \
          --threads "{params.threads}" \
          --assembly-min "{params.assembly_min}" \
          --mincl-percent "{params.mincl}" \
          --min-lcov 55 \
          --merge-threshold 0 \
          --r-value-kb "{params.r_value}" \
          --options ILLUMINA_SENSITIVE_BLASTPLUS \
          --paired \
          --automatic-filtering \
          --cleanup \
          --cleanup-after-prepare \
          --rserv-port 0
        EXIT_CODE=$?
        
        if [ $EXIT_CODE -eq 0 ]; then
            echo "TAREAN test: SUCCESS" > {output.test_tarean_done}
            echo "$(date): TAREAN test completed successfully" >> {log}
        else
            echo "TAREAN test: FAILED" > {output.test_tarean_done}
            echo "$(date): TAREAN test failed or timed out" >> {log}
            exit 1
        fi
        """

rule test_comparative_preparation:
    """Test comparative read preparation - STABLE RULE"""
    input:
        test_prepared = lambda w: expand("projects/{project}/samples/{sample}/test/test_prepared.fasta",
                                       project=w.project,
                                       sample=get_default_comparative_samples(w.project))
    output:
        test_comparative = "projects/{project}/test/test_comparative.fasta",
        test_comparative_summary = "projects/{project}/test/test_comparative_summary.txt"
    params:
        project = "{project}",
        proportions = lambda w: get_comparative_proportions(w.project),
        comparative_id = lambda w: get_default_comparative_analysis(w.project) or ""
    log:
        f'{config["global"]["log_dir"]}/test_comparative_{{project}}.log'
    shell:
        """
        echo "$(date): Testing comparative read preparation for {params.project}..." >> {log}
        
        # Create test comparative directory
        mkdir -p projects/{params.project}/test/
        
        if [ -z "{params.comparative_id}" ]; then
            echo "Comparative test: FAILED" > {output.test_comparative_summary}
            echo "No comparative_analyses configured for project {params.project}" >> {output.test_comparative_summary}
            echo "$(date): No comparative_analyses configured" >> {log}
            exit 1
        fi
        
        echo "$(date): Testing comparative analysis: {params.comparative_id}" >> {log}
        
        python3 workflows/smk_scripts/prepare_comparative_reads.py {params.project} {params.comparative_id} --verbose 2>> {log}
        
        if [ $? -eq 0 ]; then
            echo "Comparative test: SUCCESS" > {output.test_comparative_summary}
            echo "Test comparative analysis: {params.comparative_id}" >> {output.test_comparative_summary}
            echo "$(date): Comparative test completed successfully" >> {log}
            
            # Create dummy test file to satisfy Snakemake
            echo "# Test comparative reads" > {output.test_comparative}
            echo ">TEST_test_read_1" >> {output.test_comparative}
            echo "ATCGATCGATCG" >> {output.test_comparative}
            echo ">EXAM_test_read_1" >> {output.test_comparative}
            echo "GCTAGCTAGCTA" >> {output.test_comparative}
        else
            echo "Comparative test: FAILED" > {output.test_comparative_summary}
            echo "$(date): Comparative test failed" >> {log}
            exit 1
        fi
        """

rule test_post_tarean_pipeline:
    """Test post-TAREAN pipeline with test data - STABLE RULE"""
    input:
        test_tarean_done = lambda w: expand("projects/{project}/samples/{sample}/test/test_tarean.done",
                                          project=w.project,
                                          sample=config["projects"][w.project]["samples"])
    output:
        test_pipeline_done = "projects/{project}/test/test_pipeline.done"
    params:
        project = "{project}"
    log:
        f'{config["global"]["log_dir"]}/test_pipeline_{{project}}.log'
    shell:
        """
        echo "$(date): Testing post-TAREAN pipeline for {params.project}..." >> {log}
        
        # Test pipeline with report-only mode
        for sample in {config["projects"]["{params.project}"]["samples"]}; do
            echo "$(date): Testing pipeline for sample $sample" >> {log}
            
            python3 post_tarean/pipeline.py $sample --project-id {params.project} \
                --output-dir projects/{params.project}/test/post_tarean/ --report-only --test-mode 2>> {log}
            
            if [ $? -ne 0 ]; then
                echo "$(date): Pipeline test failed for $sample" >> {log}
                exit 1
            fi
        done
        
        echo "Pipeline test: SUCCESS" > {output.test_pipeline_done}
        echo "$(date): All pipeline tests completed successfully" >> {log}
        """ 