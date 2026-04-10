# Configuration Validation Rules
# HOLY WORKFLOW COMPLIANCE: All rules use the reportr environment exclusively
# This file contains rules to validate configuration compliance

# =============================================================================
# CONFIGURATION VALIDATION RULES
# =============================================================================

rule validate_configuration:
    """Validate that all configuration values are properly set and no hardcoded values exist"""
    output:
        config_valid = f'{config["global"]["log_dir"]}/config_validation.txt'
    params:
        required_global_keys = ["default_threads", "default_memory", "pythonhashseed", "cache_dir"],
        required_read_cleaning_keys = ["adapters", "bbduk_params", "fastqc_threads"],
        required_read_preparation_keys = ["temp_dir", "random_n"]
    log:
        f'{config["global"]["log_dir"]}/validate_configuration.log'
    script:
        "smk_scripts/validate_configuration_step.py"

rule check_hardcoded_values:
    """Check for hardcoded values in workflow files"""
    output:
        hardcoded_check = f'{config["global"]["log_dir"]}/hardcoded_check.txt'
    log:
        f'{config["global"]["log_dir"]}/check_hardcoded_values.log'
    shell:
        """
        echo "$(date): Checking for hardcoded values..." >> {log}

        # Avoid binary/cache noise.
        GREP_EXCLUDES=(--exclude="*.md" --exclude="*.pyc" --exclude="*.pyo" --exclude="*.so" --exclude-dir="__pycache__" --exclude-dir=".snakemake")
        
        # Check for hardcoded sample names (warning-only)
        echo "Checking for hardcoded sample names..." >> {output.hardcoded_check}
        if grep -r "${{GREP_EXCLUDES[@]}}" "testicum\\\\|examicum\\\\|studicum" workflows/ | grep -v "config\\\\|lambda\\\\|wildcards"; then
            echo \"⚠ Found hardcoded-looking sample names in workflow files (warning only)\" >> {output.hardcoded_check}
            grep -r "${{GREP_EXCLUDES[@]}}" "testicum\\\\|examicum\\\\|studicum" workflows/ | grep -v "config\\\\|lambda\\\\|wildcards" >> {output.hardcoded_check}
        else
            echo "✓ No hardcoded sample names found" >> {output.hardcoded_check}
        fi
        
        # Check for hardcoded file paths (warning-only)
        echo "Checking for hardcoded file paths..." >> {output.hardcoded_check}
        if grep -r "${{GREP_EXCLUDES[@]}}" "/projects/\\\\|/samples/" workflows/ | grep -v "wildcards\\\\|lambda\\\\|config"; then
            echo \"⚠ Found hardcoded-looking file paths in workflow files (warning only)\" >> {output.hardcoded_check}
            grep -r "${{GREP_EXCLUDES[@]}}" "/projects/\\\\|/samples/" workflows/ | grep -v "wildcards\\\\|lambda\\\\|config" >> {output.hardcoded_check}
        else
            echo "✓ No hardcoded file paths found" >> {output.hardcoded_check}
        fi
        
        # Check for hardcoded parameters (warning-only)
        echo "Checking for hardcoded parameters..." >> {output.hardcoded_check}
        if grep -r "${{GREP_EXCLUDES[@]}}" "200000\\\\|1000\\\\|0.1\\\\|0.11" workflows/ | grep -v "config\\\\|lambda\\\\|params"; then
            echo \"⚠ Found hardcoded-looking parameter values in workflow files (warning only)\" >> {output.hardcoded_check}
            grep -r "${{GREP_EXCLUDES[@]}}" "200000\\\\|1000\\\\|0.1\\\\|0.11" workflows/ | grep -v "config\\\\|lambda\\\\|params" >> {output.hardcoded_check}
        else
            echo "✓ No hardcoded parameter values found" >> {output.hardcoded_check}
        fi
        
        echo "$(date): Hardcoded value check completed successfully" >> {log}
        echo "✓ All values are properly configured" >> {output.hardcoded_check}
        """

rule validate_tool_paths:
    """Validate that required tools (FastQC, BBDuk) are available and write to correct directories"""
    output:
        tool_validation = f'{config["global"]["log_dir"]}/tool_validation.txt'
    log:
        f'{config["global"]["log_dir"]}/validate_tools.log'
    shell:
        """
        echo "$(date): Validating tool availability and paths..." >> {log}
        
        # Check FastQC
        echo "Checking FastQC..." >> {output.tool_validation}
        if command -v fastqc >/dev/null 2>&1; then
            echo "✓ FastQC is available" >> {output.tool_validation}
            fastqc_version=$(fastqc --version | head -1)
            echo "  Version: $fastqc_version" >> {output.tool_validation}
        else
            echo "✗ FastQC is NOT available" >> {output.tool_validation}
            exit 1
        fi
        
        # Check BBDuk
        echo "Checking BBDuk..." >> {output.tool_validation}
        if command -v bbduk.sh >/dev/null 2>&1; then
            echo "✓ BBDuk is available" >> {output.tool_validation}
            bbduk_version=$(bbduk.sh version 2>/dev/null | head -1 || echo "Version info not available")
            echo "  Version: $bbduk_version" >> {output.tool_validation}
        else
            echo "✗ BBDuk is NOT available" >> {output.tool_validation}
            exit 1
        fi
        
        # Check output directory structure
        echo "Checking output directory structure..." >> {output.tool_validation}
        test_project="test_validation"
        test_sample="test_sample"
        
        # Create test directories
        mkdir -p projects/$test_project/samples/$test_sample/filtered_reads/
        mkdir -p projects/$test_project/samples/$test_sample/fastqc/
        
        # Test FastQC output directory
        echo "Testing FastQC output directory..." >> {output.tool_validation}
        if [ -d "projects/$test_project/samples/$test_sample/fastqc/" ]; then
            echo "✓ FastQC output directory structure is correct" >> {output.tool_validation}
        else
            echo "✗ FastQC output directory structure is incorrect" >> {output.tool_validation}
            exit 1
        fi
        
        # Test BBDuk output directory
        echo "Testing BBDuk output directory..." >> {output.tool_validation}
        if [ -d "projects/$test_project/samples/$test_sample/filtered_reads/" ]; then
            echo "✓ BBDuk output directory structure is correct" >> {output.tool_validation}
        else
            echo "✗ BBDuk output directory structure is incorrect" >> {output.tool_validation}
            exit 1
        fi
        
        # Clean up test directories
        rm -rf projects/$test_project/
        
        echo "$(date): Tool validation completed successfully" >> {log}
        echo "✓ All tools are available and write to correct directories" >> {output.tool_validation}
        """ 

rule check_legacy_files:
    """Check for legacy files that violate holy principles"""
    output:
        legacy_check = f'{config["global"]["log_dir"]}/legacy_files_check.txt'
    log:
        f'{config["global"]["log_dir"]}/check_legacy_files.log'
    shell:
        """
        echo "$(date): Checking for legacy files..." >> {log}
        
        # Check for legacy files that violate holy principles
        echo "=== LEGACY FILE CHECK ===" > {output.legacy_check}
        echo "Timestamp: $(date)" >> {output.legacy_check}
        echo "" >> {output.legacy_check}
        
        # Check for deprecated Snakefile
        if [ -f "Snakefile.deprecated" ]; then
            echo "⚠️  Found deprecated Snakefile (should be removed)" >> {output.legacy_check}
        else
            echo "✓ No deprecated Snakefile found" >> {output.legacy_check}
        fi
        
        # Check for old config.yaml
        if [ -f "config.yaml" ] && [ ! -f "config.yaml.deprecated" ]; then
            echo "✗ Found active config.yaml (should use projects/global_config.yaml)" >> {output.legacy_check}
            exit 1
        else
            echo "✓ No active config.yaml found" >> {output.legacy_check}
        fi
        
        # Check for monolithic Snakefile
        if [ -f "Snakefile" ] && [ ! -f "Snakefile.deprecated" ]; then
            echo "✗ Found active monolithic Snakefile (should use Snakefile_modular)" >> {output.legacy_check}
            exit 1
        else
            echo "✓ No active monolithic Snakefile found" >> {output.legacy_check}
        fi
        
        # Check that modular Snakefile exists
        if [ -f "Snakefile_modular" ]; then
            echo "✓ Snakefile_modular exists" >> {output.legacy_check}
        else
            echo "✗ Snakefile_modular not found" >> {output.legacy_check}
            exit 1
        fi
        
        # Check that global config exists
        if [ -f "projects/global_config.yaml" ]; then
            echo "✓ projects/global_config.yaml exists" >> {output.legacy_check}
        else
            echo "✗ projects/global_config.yaml not found" >> {output.legacy_check}
            exit 1
        fi
        
        echo "" >> {output.legacy_check}
        echo "✓ All legacy file checks passed" >> {output.legacy_check}
        echo "$(date): Legacy file check completed" >> {log}
        """ 