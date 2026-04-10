#!/bin/bash
# RepOrtR Cleanup Script
# This script provides easy cleanup commands for runtime files and installation artifacts

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to show usage
show_usage() {
    echo "RepOrtR Cleanup Script"
    echo ""
    echo "Usage: $0 [OPTION]"
    echo ""
    echo "Options:"
    echo "  runtime     Clean runtime files (outputs, temp files, cache)"
    echo "  install     Clean installation artifacts (repex_tarean, scripts)"
    echo "  all         Clean both runtime and installation files"
    echo "  help        Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 runtime    # Clean only runtime files"
    echo "  $0 install    # Clean only installation artifacts"
    echo "  $0 all        # Clean everything (full reset)"
    echo ""
    echo "Runtime files include:"
    echo "  - cleaned_reads/"
    echo "  - fastqc_cleaned_reads/"
    echo "  - interlaced/"
    echo "  - tmp/"
    echo "  - __pycache__/"
    echo "  - *.log files"
    echo ""
    echo "Installation artifacts include:"
    echo "  - repex_tarean/"
    echo "  - activate_reportr.sh"
    echo "  - environment_*.yml files"
}

# Function to clean runtime files
clean_runtime() {
    print_status "Cleaning runtime files..."
    
    # Directories to clean
    runtime_dirs=(
        "cleaned_reads"
        "fastqc_cleaned_reads"
        "interlaced"
        "tmp"
        "__pycache__"
        "tests/cleaned_reads"
        "tests/interlaced"
        "tests/tmp"
    )
    
    # Files to clean
    runtime_files=(
        "*.log"
        "*.tmp"
        "*.temp"
        "*.bak"
        "*.backup"
    )
    
    # Clean directories
    for dir in "${runtime_dirs[@]}"; do
        if [[ -d "$dir" ]]; then
            print_status "Removing directory: $dir"
            rm -rf "$dir"
            print_success "Removed: $dir"
        else
            print_warning "Directory not found: $dir"
        fi
    done
    
    # Clean files
    for pattern in "${runtime_files[@]}"; do
        if ls $pattern 1> /dev/null 2>&1; then
            print_status "Removing files matching: $pattern"
            rm -f $pattern
            print_success "Removed files matching: $pattern"
        fi
    done
    
    print_success "Runtime cleanup completed!"
}

# Function to clean installation artifacts
clean_install() {
    print_status "Cleaning installation artifacts..."
    
    # Installation artifacts to clean
    install_items=(
        "repex_tarean"
        "activate_reportr.sh"
        "environment_repeatexplorer.yml"
        "environment_reportr.yml"
    )
    
    for item in "${install_items[@]}"; do
        if [[ -e "$item" ]]; then
            print_status "Removing: $item"
            rm -rf "$item"
            print_success "Removed: $item"
        else
            print_warning "Not found: $item"
        fi
    done
    
    print_success "Installation cleanup completed!"
}

# Function to clean everything
clean_all() {
    print_warning "This will remove ALL runtime files and installation artifacts!"
    print_warning "This is equivalent to a complete reset of the RepOrtR installation."
    echo ""
    read -p "Are you sure you want to continue? (y/N): " -n 1 -r
    echo ""
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        clean_runtime
        clean_install
        print_success "Complete cleanup finished!"
        print_status "To reinstall, run: python3 install_reportr.py"
    else
        print_status "Cleanup cancelled."
    fi
}

# Main script logic
case "${1:-help}" in
    "runtime")
        clean_runtime
        ;;
    "install")
        clean_install
        ;;
    "all")
        clean_all
        ;;
    "help"|"-h"|"--help")
        show_usage
        ;;
    *)
        print_error "Unknown option: $1"
        echo ""
        show_usage
        exit 1
        ;;
esac
