#!/usr/bin/env python3
"""
High-level test runner for RepOrtR pipeline
Runs consolidated tests with consistent formatting

Reference: See tests/TEST_STYLE_GUIDE.md for formatting standards
"""

import os
import sys
import time

# Consistent formatting constants
SUCCESS = "✅"
ERROR = "❌"
WARNING = "⚠️"
INFO = "ℹ️"

def print_header(title):
    """Print consistent header"""
    print(f"\n{'='*80}")
    print(f"{title}")
    print(f"{'='*80}")

def check_dependencies():
    """Check if required dependencies are available"""
    print("Checking dependencies...")
    
    dependencies = {
        "Biopython": "import Bio",
        "PyYAML": "import yaml",
        "FastQC": "fastqc --version",
        "BBDuk": "bbduk.sh version",
        "Snakemake": "snakemake --version"
    }
    
    available = {}
    
    # Check Python dependencies
    for name, import_cmd in [("Biopython", "import Bio"), ("PyYAML", "import yaml")]:
        try:
            exec(import_cmd)
            available[name] = True
            print(f"{SUCCESS} {name}")
        except ImportError:
            available[name] = False
            print(f"{ERROR} {name}")
    
    # Check external tools
    for name in ["FastQC", "BBDuk", "Snakemake"]:
        try:
            import subprocess
            if name == "FastQC":
                result = subprocess.run(["fastqc", "--version"], capture_output=True, text=True, timeout=10)
            elif name == "BBDuk":
                result = subprocess.run(["bbduk.sh", "version"], capture_output=True, text=True, timeout=10)
            elif name == "Snakemake":
                result = subprocess.run(["snakemake", "--version"], capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0:
                available[name] = True
                print(f"{SUCCESS} {name}")
            else:
                available[name] = False
                print(f"{ERROR} {name}")
        except:
            available[name] = False
            print(f"{ERROR} {name}")
    
    return available

def run_consolidated_test():
    """Run the consolidated test suite"""
    print_header("Running Consolidated Test Suite")
    
    start_time = time.time()
    
    try:
        # Import and run the consolidated test
        sys.path.insert(0, 'tests')
        import test_consolidated
        
        success = test_consolidated.main()
        duration = time.time() - start_time
        
        if success:
            print(f"\n{SUCCESS} Consolidated test PASSED ({duration:.2f}s)")
            return True, duration
        else:
            print(f"\n{ERROR} Consolidated test FAILED ({duration:.2f}s)")
            return False, duration
            
    except Exception as e:
        duration = time.time() - start_time
        print(f"\n{ERROR} Consolidated test ERROR: {e}")
        return False, duration

def main():
    """Run the consolidated test suite"""
    print_header("RepOrtR Pipeline - Test Suite")
    
    # Check dependencies
    deps = check_dependencies()
    
    # Run consolidated test
    passed, duration = run_consolidated_test()
    
    # Print summary
    print_header("Final Summary")
    print(f"Consolidated test: {'PASSED' if passed else 'FAILED'} ({duration:.2f}s)")
    
    print(f"\nDependency Status:")
    for dep, available in deps.items():
        status = f"{SUCCESS} Available" if available else f"{ERROR} Missing"
        print(f"  {dep:<15} {status}")
    
    # Recommendations
    print(f"\nRecommendations:")
    if not deps.get("Biopython", False):
        print(f"  - Install Biopython: conda install -c conda-forge biopython")
    if not deps.get("PyYAML", False):
        print(f"  - Install PyYAML: conda install -c conda-forge pyyaml")
    if not deps.get("Snakemake", False):
        print(f"  - Install Snakemake: conda install -c conda-forge snakemake")
    if not deps.get("BBDuk", False):
        print(f"  - Install BBDuk: conda install -c bioconda bbmap")
    if not deps.get("FastQC", False):
        print(f"  - Install FastQC: conda install -c bioconda fastqc")
    
    if passed:
        print(f"\n{SUCCESS} All tests passed! Pipeline is ready for use.")
        return True
    else:
        print(f"\n{WARNING} Some tests failed. Check the output above for details.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 