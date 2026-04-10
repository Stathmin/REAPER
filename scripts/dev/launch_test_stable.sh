#!/bin/bash
# Launch script for probename project with stability features
# Handles server disconnections and ensures assembly stability
# 🏛️ HOLY PRINCIPLE: Comprehensive logging for development analysis

set -e

PROJECT_ID="probename_project"
CORES=28
MEMORY="180G"
MAX_MEMORY_PERCENT=80
CHECKPOINT_INTERVAL=3600

# 🏛️ HOLY PRINCIPLE: Development logging setup (config-driven)
LOG_DIR=$(python3 -c 'import yaml; cfg=yaml.safe_load(open("projects/global_config.yaml")) or {}; print((cfg.get("global") or {}).get("log_dir","logs"))')
DEV_LOG_DIR="${LOG_DIR}/performance_logs"
RUN_ID=$(date +%Y%m%d_%H%M%S)
DEV_LOG_FILE="${DEV_LOG_DIR}/pipeline_run_${RUN_ID}.log"
PERFORMANCE_LOG="${DEV_LOG_DIR}/performance_${RUN_ID}.json"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - $1" >> "$DEV_LOG_FILE"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - ERROR: $1" >> "$DEV_LOG_FILE"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - WARNING: $1" >> "$DEV_LOG_FILE"
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
    echo "$(date +'%Y-%m-%d %H:%M:%S') - SUCCESS: $1" >> "$DEV_LOG_FILE"
}

# 🏛️ HOLY PRINCIPLE: Input characteristics logging
log_input_characteristics() {
    log "🏛️ HOLY PRINCIPLE: Logging input characteristics for development analysis..."
    
    # Create performance log structure
    cat > "$PERFORMANCE_LOG" << EOF
{
    "run_id": "$RUN_ID",
    "timestamp": "$(date -Iseconds)",
    "project_id": "$PROJECT_ID",
    "input_characteristics": {
        "samples": {}
    },
    "system_resources": {
        "total_memory_gb": $TOTAL_MEMORY,
        "available_memory_gb": $AVAILABLE_MEMORY,
        "total_cores": $TOTAL_CORES,
        "requested_cores": $CORES,
        "requested_memory": "$MEMORY"
    },
    "pipeline_config": {
        "max_memory_percent": $MAX_MEMORY_PERCENT,
        "checkpoint_interval": $CHECKPOINT_INTERVAL
    },
    "runtime_data": {
        "start_time": "$(date -Iseconds)",
        "end_time": null,
        "duration_seconds": null,
        "resource_usage": {}
    }
}
EOF
    
    # Log input characteristics for each sample
    for sample in testicum examicum studicum; do
        r1_file="data/probename_reads/${sample}_R1.fq"
        r2_file="data/probename_reads/${sample}_R2.fq"
        
        if [[ -f "$r1_file" && -f "$r2_file" ]]; then
            # Get file sizes
            r1_size=$(stat -c%s "$r1_file")
            r2_size=$(stat -c%s "$r2_file")
            
            # Count reads (approximate)
            r1_reads=$(echo "$r1_size / 100" | bc -l | cut -d. -f1)  # Rough estimate
            
            # Get genome size from config
            genome_size=$(python3 -c "
import yaml
with open('projects/global_config.yaml', 'r') as f:
    config = yaml.safe_load(f)
print(config['projects']['$PROJECT_ID']['samples']['$sample']['genome_size'])
" 2>/dev/null || echo "1.0")
            
            # Update performance log
            python3 -c "
import json
with open('$PERFORMANCE_LOG', 'r') as f:
    data = json.load(f)
data['input_characteristics']['samples']['$sample'] = {
    'r1_file': '$r1_file',
    'r2_file': '$r2_file',
    'r1_size_bytes': $r1_size,
    'r2_size_bytes': $r2_size,
    'estimated_reads': $r1_reads,
    'genome_size': $genome_size,
    'total_size_mb': round(($r1_size + $r2_size) / 1024 / 1024, 2)
}
with open('$PERFORMANCE_LOG', 'w') as f:
    json.dump(data, f, indent=2)
"
            
            log "📊 Sample $sample: ${r1_reads} reads, ${genome_size} genome size, $((($r1_size + $r2_size) / 1024 / 1024))MB total"
        else
            warning "Sample $sample files not found: $r1_file, $r2_file"
        fi
    done
    
    success "Input characteristics logged to $PERFORMANCE_LOG"
}

# 🏛️ HOLY PRINCIPLE: Runtime monitoring
monitor_runtime() {
    # Start time will be set when actual assembly begins
    local start_time=""
    
    # Function to start timing (called when assembly actually starts)
    start_timing() {
        start_time=$(date +%s)
        log "🏛️ HOLY PRINCIPLE: Assembly timing started at $(date)"
    }
    
    # Function to stop monitoring
    stop_monitoring() {
        # Only calculate duration if timing was started
        if [[ -n "$start_time" ]]; then
            local end_time=$(date +%s)
            local duration=$((end_time - start_time))
        else
            # If timing wasn't started, use a reasonable default
            duration=0
            log "🏛️ HOLY PRINCIPLE: No assembly timing data available"
        fi
        
        python3 -c "
import json
import time
try:
    with open('$PERFORMANCE_LOG', 'r') as f:
        data = json.load(f)
    
    data['runtime_data']['end_time'] = time.strftime('%Y-%m-%dT%H:%M:%S')
    data['runtime_data']['duration_seconds'] = $duration
    
    with open('$PERFORMANCE_LOG', 'w') as f:
        json.dump(data, f, indent=2)
except Exception as e:
    pass
"
        
        # Format duration nicely with proper handling (cap display, not actual value)
        local display_duration=$duration
        if [[ $display_duration -gt 86400 ]]; then
            # For display purposes only, show max 24 hours but log actual value
            display_duration=86400
            log "🏛️ HOLY PRINCIPLE: Display duration capped at 24 hours, actual duration: ${duration}s"
        fi
        
        if [[ $display_duration -lt 60 ]]; then
            duration_str="${display_duration}s"
        elif [[ $display_duration -lt 3600 ]]; then
            minutes=$((display_duration / 60))
            seconds=$((display_duration % 60))
            duration_str="${minutes}m ${seconds}s"
        elif [[ $display_duration -lt 86400 ]]; then
            hours=$((display_duration / 3600))
            minutes=$(((display_duration % 3600) / 60))
            seconds=$((display_duration % 60))
            duration_str="${hours}h ${minutes}m ${seconds}s"
        else
            days=$((display_duration / 86400))
            hours=$(((display_duration % 86400) / 3600))
            minutes=$(((display_duration % 3600) / 60))
            duration_str="${days}d ${hours}h ${minutes}m"
        fi
        
        log "🏛️ HOLY PRINCIPLE: Runtime monitoring completed. Duration: ${duration_str} (actual: ${duration}s)"
    }
    
    # Set up cleanup
    trap stop_monitoring EXIT
}

# Check if we're in the right environment
if [[ "$CONDA_DEFAULT_ENV" != "reportr" ]]; then
    error "This script must be run from the 'reportr' conda environment"
    echo "   Activate it with: conda activate reportr"
    exit 1
fi

# Check if required files exist
if [[ ! -f "Snakefile_modular" ]]; then
    error "Snakefile_modular not found in current directory"
    exit 1
fi

# Check if modular workflow files exist
for module in core_rules.smk test_rules.smk analysis_rules.smk report_rules.smk config_validation.smk; do
    if [[ ! -f "workflows/$module" ]]; then
        error "Modular workflow file workflows/$module not found"
        exit 1
    fi
done

if [[ ! -f "monitoring/stability_monitor.py" ]]; then
    error "monitoring/stability_monitor.py not found"
    exit 1
fi

# 🏛️ HOLY PRINCIPLE: Create development logging directories
log "🏛️ HOLY PRINCIPLE: Setting up development logging infrastructure..."
mkdir -p "$DEV_LOG_DIR"
mkdir -p "${LOG_DIR}"
mkdir -p tmp
mkdir -p .snakemake

# Check system resources
log "Checking system resources..."
TOTAL_MEMORY=$(free -g | awk 'NR==2{print $2}')
AVAILABLE_MEMORY=$(free -g | awk 'NR==2{print $7}')
TOTAL_CORES=$(nproc)

log "System resources:"
log "  Total memory: ${TOTAL_MEMORY}GB"
log "  Available memory: ${AVAILABLE_MEMORY}GB"
log "  Total cores: ${TOTAL_CORES}"

if [[ $AVAILABLE_MEMORY -lt 64 ]]; then
    warning "Available memory (${AVAILABLE_MEMORY}GB) is less than recommended (64GB)"
fi

if [[ $TOTAL_CORES -lt $CORES ]]; then
    warning "Available cores (${TOTAL_CORES}) is less than requested (${CORES})"
    CORES=$TOTAL_CORES
fi

# 🏛️ HOLY PRINCIPLE: Log input characteristics
log_input_characteristics

# Function to handle graceful shutdown
cleanup() {
    log "Received shutdown signal, cleaning up..."
    # Kill any background processes
    jobs -p | xargs -r kill
    success "Cleanup completed"
    exit 0
}

# Set up signal handlers
trap cleanup SIGTERM SIGINT

# 🏛️ HOLY PRINCIPLE: Check and clean up lock files before execution
log "🏛️ HOLY PRINCIPLE: Checking for existing lock files and cleaning up..."
if [[ -d ".snakemake/locks" ]]; then
    warning "Found existing Snakemake locks, cleaning up..."
    rm -rf .snakemake/locks
    log "Cleaned up existing locks"
fi

if [[ -d ".snakemake/incomplete" ]]; then
    warning "Found incomplete Snakemake state, cleaning up..."
    rm -rf .snakemake/incomplete
    log "Cleaned up incomplete state"
fi

# Check for existing checkpoint
if [[ -f ".snakemake/checkpoint_${PROJECT_ID}.yaml" ]]; then
    log "Found existing checkpoint, will resume from previous state"
fi

# Run with stability monitoring
run_with_stability() {
    local target="$1"
    
    log "Starting probename project with stability features..."
    log "Project: ${PROJECT_ID}"
    log "Cores: ${CORES}"
    log "Memory: ${MEMORY}"
    log "Max memory: ${MAX_MEMORY_PERCENT}%"
    log "Checkpoint interval: ${CHECKPOINT_INTERVAL}s"
    
    if [[ -n "$target" ]]; then
        log "Target: ${target}"
    fi
    
    # 🏛️ HOLY PRINCIPLE: Ensure complete pipeline execution with all dependencies
    log "🏛️ HOLY PRINCIPLE: Building complete pipeline with all preprocessing steps..."
    
    # First run system validation tests
    log "Running system validation tests..."
    snakemake -s Snakefile_modular --configfile projects/global_config.yaml --cores 4 "${LOG_DIR}/system_resources_test.txt" --quiet > /dev/null 2>&1
    
    if [[ $? -ne 0 ]]; then
        error "System resource test failed"
        exit 1
    fi
    
    # 🏛️ HOLY PRINCIPLE: Build preprocessing dependencies for target only
    log "🏛️ HOLY PRINCIPLE: Building preprocessing dependencies for target..."
    
    # Determine which samples need preprocessing based on target
    if [[ -n "$target" ]]; then
        case "$target" in
            all)
                # Build for all samples
                target_samples="testicum examicum studicum"
                ;;
            comparative)
                # Build for comparative analysis samples
                target_samples="testicum examicum"
                ;;
            testicum|examicum|studicum)
                # Build for specific sample only
                target_samples="$target"
                ;;
            *)
                # Default to all samples
                target_samples="testicum examicum studicum"
                ;;
        esac
    else
        # Default to all samples
        target_samples="testicum examicum studicum"
    fi
    
    log "Building dependencies for samples: $target_samples"
    
    # Build cleaned reads for target samples only
    for sample in $target_samples; do
        log "Building cleaned reads for ${sample}..."
        snakemake -s Snakefile_modular --configfile projects/global_config.yaml --cores 4 "projects/${PROJECT_ID}/samples/${sample}/filtered_reads/cleaned.done" --quiet > /dev/null 2>&1
        
        if [[ $? -ne 0 ]]; then
            error "Failed to build cleaned reads for ${sample}"
            exit 1
        fi
    done
    
    # Build prepared reads for target samples only
    for sample in $target_samples; do
        log "Building prepared reads for ${sample}..."
        snakemake -s Snakefile_modular --configfile projects/global_config.yaml --cores 4 "projects/${PROJECT_ID}/samples/${sample}/tarean/prepared_forRE.fasta" --quiet > /dev/null 2>&1
        
        if [[ $? -ne 0 ]]; then
            error "Failed to build prepared reads for ${sample}"
            exit 1
        fi
    done
    
    success "Preprocessing dependencies built for target samples! Starting assembly..."
    
    # Run with stability monitor - specify actual assembly targets
    if [[ -n "$target" ]]; then
        case "$target" in
            all)
                # Run all individual samples and comparative analysis
                assembly_target="projects/${PROJECT_ID}/samples/testicum/tarean/tarean.done projects/${PROJECT_ID}/samples/examicum/tarean/tarean.done projects/${PROJECT_ID}/samples/studicum/tarean/tarean.done projects/${PROJECT_ID}/comparative/comp_testicum_examicum/comparative_reads.fasta"
                ;;
            comparative)
                # Run only comparative analysis
                assembly_target="projects/${PROJECT_ID}/comparative/comp_testicum_examicum/comparative_reads.fasta"
                ;;
            testicum|examicum|studicum)
                # Run specific sample
                assembly_target="projects/${PROJECT_ID}/samples/${target}/tarean/tarean.done"
                ;;
            *)
                # Use target as-is
                assembly_target="$target"
                ;;
        esac
    else
        # Default to running TAREAN analysis for all samples and comparative analysis
        assembly_target="projects/${PROJECT_ID}/samples/testicum/tarean/tarean.done projects/${PROJECT_ID}/samples/examicum/tarean/tarean.done projects/${PROJECT_ID}/samples/studicum/tarean/tarean.done projects/${PROJECT_ID}/comparative/comp_testicum_examicum/comparative_reads.fasta"
    fi
    
    log "🏛️ HOLY PRINCIPLE: Launching complete pipeline with stability monitoring..."
    
    # Start timing when assembly actually begins
    start_timing
    
    # Start stability monitor and wait for completion
    log "Starting stability monitor with target: $assembly_target"
    python3 monitoring/stability_monitor.py \
        "${PROJECT_ID}" \
        --cores "${CORES}" \
        --memory "${MEMORY}" \
        --max-memory "${MAX_MEMORY_PERCENT}" \
        --checkpoint-interval "${CHECKPOINT_INTERVAL}" \
        --target "$assembly_target"
    
    # Check if stability monitor completed successfully
    if [[ $? -eq 0 ]]; then
        success "Stability monitor completed successfully"
    else
        error "Stability monitor failed"
        exit 1
    fi
}

# Main execution
main() {
    local target=""
    
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --target)
                target="$2"
                shift 2
                ;;
            --cores)
                CORES="$2"
                shift 2
                ;;
            --memory)
                MEMORY="$2"
                shift 2
                ;;
            --help|-h)
                echo "Usage: $0 [OPTIONS]"
                echo ""
                echo "Options:"
                echo "  --target TARGET    Specific target to run"
                echo "  --cores N          Number of cores (default: 16)"
                echo "  --memory SIZE      Memory limit (default: 80G)"
                echo "  --help, -h         Show this help message"
                echo ""
                echo "Examples:"
                echo "  $0                                    # Run all targets (individual + comparative)"
                echo "  $0 --target all                      # Run all targets (individual + comparative)"
                echo "  $0 --target comparative              # Run only comparative analysis"
                echo "  $0 --target testicum                 # Run only testicum"
                echo "  $0 --target examicum                 # Run only examicum"
                echo "  $0 --target studicum                 # Run only studicum"
                echo "  $0 --cores 8 --memory 40G           # Use 8 cores, 40GB memory"
                exit 0
                ;;
            *)
                error "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done
    
    # Validate target if specified
    if [[ -n "$target" ]]; then
                    case "$target" in
                testicum|examicum|studicum|all|comparative)
                    log "Target validated: $target"
                    ;;
                *)
                    error "Invalid target: $target"
                    echo "Valid targets: testicum, examicum, studicum, all, comparative"
                    echo "  - all: Run all individual samples + comparative analysis"
                    echo "  - comparative: Run only comparative analysis (testicum vs examicum)"
                    echo "  - testicum/examicum/studicum: Run specific sample"
                    exit 1
                    ;;
            esac
    fi
    
    # 🏛️ HOLY PRINCIPLE: Start runtime monitoring
    log "🏛️ HOLY PRINCIPLE: Starting comprehensive runtime monitoring for development analysis..."
    monitor_runtime
    
    # Start the assembly
    success "Starting probename assembly with stability features and comprehensive logging"
    run_with_stability "$target"
}

# Run main function
main "$@" 