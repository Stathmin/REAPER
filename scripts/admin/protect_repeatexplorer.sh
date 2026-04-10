#!/bin/bash
# Protect Repeatexplorer Environment Script
# This script prevents accidental removal of the critical repeatexplorer environment

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

# Check if repeatexplorer environment exists
if conda env list | grep -q "repeatexplorer"; then
    print_success "Repeatexplorer environment found and protected"
    
    # Create a backup of the environment
    BACKUP_FILE="repeatexplorer_backup_$(date +%Y%m%d_%H%M%S).yml"
    print_warning "Creating backup of repeatexplorer environment..."
    
    if conda env export -n repeatexplorer > "$BACKUP_FILE" 2>/dev/null; then
        print_success "Backup created: $BACKUP_FILE"
    else
        print_warning "Could not create backup, but environment exists"
    fi
    
    # Check if seqclust is working (in-repo binary under repeatexplorer env)
    if conda run -n repeatexplorer ./repex_tarean/seqclust --help >/dev/null 2>&1; then
        print_success "seqclust is working correctly (repex_tarean/seqclust)"
    else
        print_error "seqclust is not working - environment or build may be corrupted"
        print_warning "Run: python3 install_reportr.py to fix"
    fi
else
    print_error "Repeatexplorer environment NOT FOUND!"
    print_error "This will break all seqclust functionality!"
    print_warning "Run: python3 install_reportr.py to reinstall"
fi

# Check for any recent removal attempts
if [ -f "install.log" ]; then
    if grep -q "repeatexplorer.*remove" install.log; then
        print_error "Recent removal attempt detected in install.log!"
        print_warning "Check the log for details"
    fi
fi

echo ""
echo "🏛️ HOLY PRINCIPLE: Repeatexplorer environment protection active"
echo "   - Environment: $(conda env list | grep repeatexplorer || echo 'NOT FOUND')"
echo "   - Seqclust: $(conda run -n repeatexplorer ./repex_tarean/seqclust --help >/dev/null 2>&1 && echo 'WORKING' || echo 'BROKEN')"
echo "   - Backup: $BACKUP_FILE"
