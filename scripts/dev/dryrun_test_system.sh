#!/bin/bash
# Test script for probename project - validates system without running full assemblies

set -e

PROJECT_ID="probename_project"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging function
log() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
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

# Create necessary directories
log "Creating necessary directories..."
LOG_DIR=$(python3 -c 'import yaml; cfg=yaml.safe_load(open("projects/global_config.yaml")) or {}; print((cfg.get("global") or {}).get("log_dir","logs"))')
mkdir -p "${LOG_DIR}"
mkdir -p tmp
mkdir -p .snakemake

# Function to run tests
run_tests() {
    local test_name=$1
    local target=$2
    
    log "Running ${test_name}..."
    snakemake -s Snakefile_modular --configfile projects/global_config.yaml --cores 4 "$target" --quiet
    
    if [[ $? -eq 0 ]]; then
        success "${test_name} passed"
        return 0
    else
        error "${test_name} failed"
        return 1
    fi
}

# Main test execution
main() {
    log "Starting probename system tests..."
    
    # Test 1: System resources
    if ! run_tests "System resources test" "${LOG_DIR}/system_resources_test.txt"; then
        exit 1
    fi
    
    # Test 2: Read preparation for each sample
    for sample in testicum examicum studicum; do
        if ! run_tests "Read preparation test (${sample})" "projects/${PROJECT_ID}/samples/${sample}/test/test_summary.txt"; then
            exit 1
        fi
    done
    
    # Test 3: TAREAN with small datasets
    for sample in testicum examicum studicum; do
        if ! run_tests "TAREAN test (${sample})" "projects/${PROJECT_ID}/samples/${sample}/test/test_tarean.done"; then
            exit 1
        fi
    done
    
    # Test 4: Comparative preparation
    if ! run_tests "Comparative preparation test" "projects/${PROJECT_ID}/test/test_comparative_summary.txt"; then
        exit 1
    fi
    
    # Test 5: Post-TAREAN pipeline
    if ! run_tests "Post-TAREAN pipeline test" "projects/${PROJECT_ID}/test/test_pipeline.done"; then
        exit 1
    fi
    
    success "All system tests passed! The system is ready for full assembly."
    log "You can now run: bash scripts/dev/launch_test_stable.sh"
}

# Run main function
main "$@"
