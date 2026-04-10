# Progress Tracking Rules - seqclust Progress Monitoring
# HOLY WORKFLOW COMPLIANCE: All rules use the reportr environment exclusively

import json
import re
from datetime import datetime, timedelta

# =============================================================================
# PROGRESS TRACKING CONFIGURATION
# =============================================================================

# Stage mapping for seqclust stderr patterns
SEQCLUST_STAGES = {
    "initialization": {
        "patterns": [r"Initializing", r"Loading", r"Starting"],
        "description": "Initializing TAREAN analysis",
        "estimated_duration": 1800  # 30 minutes
    },
    "read_processing": {
        "patterns": [r"Processing reads", r"Read processing", r"Filtering"],
        "description": "Processing and filtering reads",
        "estimated_duration": 3600  # 1 hour
    },
    "clustering": {
        "patterns": [r"Clustering", r"Building clusters", r"Cluster analysis"],
        "description": "Building repeat clusters",
        "estimated_duration": 7200  # 2 hours
    },
    "assembly": {
        "patterns": [r"Assembly", r"Contig building", r"Assembling"],
        "description": "Assembling repeat contigs",
        "estimated_duration": 5400  # 1.5 hours
    },
    "annotation": {
        "patterns": [r"Annotation", r"BLAST", r"Database search"],
        "description": "Annotating repeats",
        "estimated_duration": 3600  # 1 hour
    },
    "finalization": {
        "patterns": [r"Finalizing", r"Writing results", r"Completing"],
        "description": "Finalizing analysis",
        "estimated_duration": 1800  # 30 minutes
    }
}

# =============================================================================
# PROGRESS TRACKING RULES
# =============================================================================

rule track_seqclust_progress:
    """Track seqclust progress with beautiful output - HOLY RULE"""
    input:
        log_file = f'{config["global"]["log_dir"]}/seqclust_{{project}}_{{sample}}.log'
    output:
        progress_file = f'{config["global"]["log_dir"]}/progress_{{project}}_{{sample}}.json',
        progress_display = f'{config["global"]["log_dir"]}/progress_{{project}}_{{sample}}.txt'
    params:
        project = "{project}",
        sample = "{sample}",
        stages = lambda w: json.dumps(SEQCLUST_STAGES),
        tracking_interval = lambda w: int(
            (config.get("global", {}) or {})
            .get("progress_tracking", {})
            .get("tracking_interval", 30)
        ),
    log:
        f'{config["global"]["log_dir"]}/progress_tracking_{{project}}_{{sample}}.log'
    shell:
        """
        python3 workflows/smk_scripts/track_seqclust_progress.py \
          --project "{params.project}" \
          --sample "{params.sample}" \
          --seqclust-log "{input.log_file}" \
          --progress-json "{output.progress_file}" \
          --progress-txt "{output.progress_display}" \
          --runner-log "{log}" \
          --stages-json '{params.stages}' \
          --interval "{params.tracking_interval}"
        """

rule display_progress:
    """Display beautiful progress for all running samples"""
    output:
        touch(f'{config["global"]["log_dir"]}/progress_display_updated.txt')
    params:
        display_interval = lambda w: int(
            (config.get("global", {}) or {})
            .get("progress_tracking", {})
            .get("display_interval", 10)
        ),
        log_dir = config["global"]["log_dir"]
    shell:
        """
        # Display progress for all samples
        echo "🎯 RepOrtR Progress Dashboard"
        echo "=" * 60
        
        for progress_file in {params.log_dir}/progress_*.json; do
            if [ -f "$progress_file" ]; then
                sample=$(basename "$progress_file" .json | sed 's/progress_//')
                if [ -f "{params.log_dir}/progress_${sample}.txt" ]; then
                    cat "{params.log_dir}/progress_${sample}.txt"
                    echo ""
                fi
            fi
        done
        
        echo "Press Ctrl+C to exit"
        sleep {params.display_interval}
        """

# =============================================================================
# PROGRESS VALIDATION RULES
# =============================================================================

rule validate_progress_tracking:
    """Validate progress tracking configuration - HOLY RULE"""
    output:
        progress_validation = f'{config["global"]["log_dir"]}/progress_tracking_validation.txt'
    log:
        f'{config["global"]["log_dir"]}/validate_progress.log'
    shell:
        """
        echo "$(date): Validating progress tracking configuration..." >> {log}
        
        # Check if progress tracking is enabled
        if [ "$(python3 -c 'import yaml; cfg=yaml.safe_load(open("projects/global_config.yaml")) or {}; print(1 if (cfg.get("global", {}) or {}).get("progress_tracking", {}).get("enabled", False) else 0)' 2>/dev/null)" -eq 1 ]; then
            echo "✅ Progress tracking is enabled" >> {output.progress_validation}
        else
            echo "⚠️  Progress tracking disabled (global.progress_tracking.enabled = false)" >> {output.progress_validation}
        fi
        
        # Validate stage patterns
        python3 -c "
import json
stages = {json.dumps(SEQCLUST_STAGES)}
for stage_name, stage_info in stages.items():
    if not stage_info.get('patterns'):
        print(f'✗ Stage {stage_name} has no patterns')
        exit(1)
    if not stage_info.get('estimated_duration'):
        print(f'✗ Stage {stage_name} has no estimated duration')
        exit(1)
print('✅ All stage configurations are valid')
"
        
        echo "$(date): Progress tracking validation completed" >> {log}
        echo "✅ Progress tracking validation passed" >> {output.progress_validation}
        """
