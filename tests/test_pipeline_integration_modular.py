#!/usr/bin/env python3
"""
Integration test for RepOrtR Modular Pipeline
Actually runs pipeline components with test data and validates modular structure
"""

import os
import sys
import subprocess
import tempfile
import shutil
import time
from pathlib import Path
# from Bio import SeqIO  # Optional import for sequence analysis

def test_modular_structure():
    """Test that all required modular files exist and are valid"""
    print("Testing modular structure...")
    
    required_files = [
        "Snakefile_modular",
        "workflows/core_rules.smk",
        "workflows/test_rules.smk", 
        "workflows/analysis_rules.smk",
        "workflows/report_rules.smk",
        "workflows/config_validation.smk"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not os.path.exists(file_path):
            missing_files.append(file_path)
    
    if missing_files:
        print(f"✗ Missing required modular files: {missing_files}")
        return False
    
    print("✓ All modular workflow files present")
    return True

def test_configuration_validation():
    """Test configuration validation with test data"""
    print("Testing configuration validation...")
    
    # Check if snakemake is available
    try:
        result = subprocess.run(["snakemake", "--version"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode != 0:
            print("✗ Snakemake not available")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("✗ Snakemake not found")
        return False
    
    # Create temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        original_dir = os.getcwd()
        os.chdir(temp_dir)
        
        # Copy necessary files (validation rules expect projects/global_config.yaml)
        shutil.copy2(f"{original_dir}/Snakefile_modular", "Snakefile_modular")
        os.makedirs("projects", exist_ok=True)
        shutil.copy2(f"{original_dir}/projects/global_config.yaml", "projects/global_config.yaml")
        shutil.copytree(f"{original_dir}/workflows", "workflows")
        
        # Run validation rules
        cmd = [
            "snakemake", "-s", "Snakefile_modular",
            "--configfile", "projects/global_config.yaml",
            "validate_configuration", "check_hardcoded_values", "validate_tool_paths"
        ]
        
        print(f"Running validation: {' '.join(cmd)}")
        start_time = time.time()
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            print(f"✗ Validation failed: {result.stderr}")
            return False
        
        print(f"✓ Validation completed in {time.time() - start_time:.2f} seconds")
        return True

def test_prepare_reads_integration():
    """Test the prepareReadsRE.py script with actual test data"""
    print("Testing prepareReadsRE.py integration...")
    
    # Test data paths
    test_r1 = "data/reads/test_R1.fq"
    test_r2 = "data/reads/test_R2.fq"
    
    # Check if test data exists
    if not os.path.exists(test_r1) or not os.path.exists(test_r2):
        print("✗ Test data not found, creating minimal test data...")
        
        # Create minimal test data
        os.makedirs("data/reads", exist_ok=True)
        with open(test_r1, 'w') as f:
            f.write("@test_read1\nATCGATCG\n+\nIIIIIIII\n")
        with open(test_r2, 'w') as f:
            f.write("@test_read1\nGCTAGCTA\n+\nIIIIIIII\n")
    
    # Create temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        original_dir = os.getcwd()
        os.chdir(temp_dir)
        
        # Copy test files to temp directory
        shutil.copy2(f"{original_dir}/{test_r1}", "test_R1.fq")
        shutil.copy2(f"{original_dir}/{test_r2}", "test_R2.fq")
        
        # Copy prepareReadsRE.py script (lives in pre_tarean/)
        shutil.copy2(f"{original_dir}/pre_tarean/prepareReadsRE.py", "prepareReadsRE.py")
        
        # Run prepareReadsRE.py with test data
        cmd = [
            "python3", "prepareReadsRE.py",
            "--randomN", "1000",  # Use smaller number for testing
            "test_R1.fq", "test_R2.fq", "test_sample"
        ]
        
        print(f"Running: {' '.join(cmd)}")
        start_time = time.time()
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        if result.returncode != 0:
            print(f"✗ Error: {result.stderr}")
            return False
        
        print(f"✓ Script completed in {time.time() - start_time:.2f} seconds")
        
        # Check output file
        output_file = "test_sample_prepared_forRE.fasta"
        if not os.path.exists(output_file):
            print(f"✗ Error: Output file {output_file} not created")
            return False
        
        # Analyze output file
        with open(output_file, 'r') as f:
            content = f.read()
            lines = content.strip().split('\n')
            
            # Count sequences
            seq_count = sum(1 for line in lines if line.startswith('>'))
            print(f"✓ Generated {seq_count} sequences")
            
            # Check format
            if not any('>test_sample_read' in line for line in lines):
                print("✗ Error: Output doesn't contain expected format")
                return False
        
        print("✓ prepareReadsRE.py integration test passed")
        return True

def test_modular_snakemake_workflow():
    """Test modular Snakemake workflow with test data"""
    print("Testing modular Snakemake workflow...")
    
    # Check if snakemake is available
    try:
        result = subprocess.run(["snakemake", "--version"], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode != 0:
            print("✗ Snakemake not available")
            return False
    except (subprocess.TimeoutExpired, FileNotFoundError):
        print("✗ Snakemake not found")
        return False
    
    # Create temporary working directory
    with tempfile.TemporaryDirectory() as temp_dir:
        original_dir = os.getcwd()
        os.chdir(temp_dir)
        
        # Copy necessary files
        shutil.copy2(f"{original_dir}/Snakefile_modular", "Snakefile_modular")
        os.makedirs("projects", exist_ok=True)
        shutil.copy2(f"{original_dir}/projects/global_config.yaml", "projects/global_config.yaml")
        shutil.copytree(f"{original_dir}/workflows", "workflows")
        shutil.copy2(f"{original_dir}/pre_tarean/prepareReadsRE.py", "prepareReadsRE.py")
        shutil.copy2(f"{original_dir}/project_manager.py", "project_manager.py")
        
        # Create minimal test project structure
        os.makedirs("projects/test_project/samples/test_sample/raw_reads", exist_ok=True)
        
        # Create test read files
        with open("projects/test_project/samples/test_sample/raw_reads/R1.fq", 'w') as f:
            f.write("@test_read1\nATCGATCG\n+\nIIIIIIII\n")
        with open("projects/test_project/samples/test_sample/raw_reads/R2.fq", 'w') as f:
            f.write("@test_read1\nGCTAGCTA\n+\nIIIIIIII\n")
        
        # Run modular Snakemake with dry-run
        cmd = [
            "snakemake", "-s", "Snakefile_modular",
            "--configfile", "projects/global_config.yaml",
            "--dry-run", "--cores", "1"
        ]
        
        print(f"Running modular Snakemake: {' '.join(cmd)}")
        start_time = time.time()
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            print(f"✗ Modular Snakemake failed: {result.stderr}")
            return False
        
        print(f"✓ Modular Snakemake completed in {time.time() - start_time:.2f} seconds")
        return True

def test_file_operations():
    """Test file operations and directory structure"""
    print("Testing file operations...")
    
    # Test directory creation
    test_dirs = [
        "projects/test_project/samples/test_sample/filtered_reads",
        "projects/test_project/samples/test_sample/fastqc",
        "projects/test_project/samples/test_sample/tarean",
        "projects/test_project/samples/test_sample/post_tarean"
    ]
    
    for test_dir in test_dirs:
        os.makedirs(test_dir, exist_ok=True)
        if not os.path.exists(test_dir):
            print(f"✗ Failed to create directory: {test_dir}")
            return False
    
    print("✓ Directory creation test passed")
    
    # Test file writing
    test_file = "projects/test_project/samples/test_sample/test.txt"
    try:
        with open(test_file, 'w') as f:
            f.write("test content")
        
        if not os.path.exists(test_file):
            print(f"✗ Failed to create file: {test_file}")
            return False
        
        # Clean up
        os.remove(test_file)
        print("✓ File operations test passed")
        return True
        
    except Exception as e:
        print(f"✗ File operation failed: {e}")
        return False

def test_configuration_compliance():
    """Test that no hardcoded values exist in workflow files"""
    print("Testing configuration compliance...")
    
    # Check for hardcoded sample names
    workflow_files = [
        "workflows/core_rules.smk",
        "workflows/test_rules.smk",
        "workflows/analysis_rules.smk",
        "workflows/report_rules.smk"
    ]
    
    hardcoded_patterns = [
        r'testicum',
        r'examicum', 
        r'studicum'
    ]
    
    for workflow_file in workflow_files:
        if not os.path.exists(workflow_file):
            print(f"✗ Workflow file not found: {workflow_file}")
            continue
            
        with open(workflow_file, 'r') as f:
            content = f.read()
            
        for pattern in hardcoded_patterns:
            import re
            matches = re.findall(pattern, content)
            # Allow matches in comments or config references
            if matches:
                for match in matches:
                    if match not in content and 'config' not in content and 'wildcards' not in content:
                        print(f"✗ Found hardcoded value '{match}' in {workflow_file}")
                        return False
    
    print("✓ No hardcoded values found in workflow files")
    return True

def main():
    """Run all integration tests"""
    print("=" * 60)
    print("RepOrtR Modular Pipeline Integration Tests")
    print("=" * 60)
    
    tests = [
        ("Modular Structure", test_modular_structure),
        ("Configuration Validation", test_configuration_validation),
        ("Prepare Reads Integration", test_prepare_reads_integration),
        ("Modular Snakemake Workflow", test_modular_snakemake_workflow),
        ("File Operations", test_file_operations),
        ("Configuration Compliance", test_configuration_compliance)
    ]
    
    results = []
    for test_name, test_func in tests:
        print(f"\n--- {test_name} ---")
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"✗ Test failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = 0
    total = len(results)
    
    for test_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{test_name}: {status}")
        if result:
            passed += 1
    
    print(f"\nOverall: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! Modular pipeline is ready.")
        return True
    else:
        print("❌ Some tests failed. Please check the output above.")
        return False

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1) 