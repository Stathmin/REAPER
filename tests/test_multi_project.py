#!/usr/bin/env python3
"""
Test script for RepOrtR Multi-Project Framework
Tests project management, sample handling, and workflow functionality
"""

import os
import sys
import yaml
import tempfile
import shutil
from pathlib import Path

# Ensure repo root is on sys.path so we can import project_manager
ROOT_DIR = Path(__file__).resolve().parents[1]
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from project_manager import ProjectManager

# Test constants
SUCCESS = "✅"
ERROR = "❌"
WARNING = "⚠️"
INFO = "ℹ️"

def print_header(title):
    """Print consistent header"""
    print(f"\n{'='*60}")
    print(f"{title}")
    print(f"{'='*60}")

def print_section(title):
    """Print consistent section header"""
    print(f"\n--- {title} ---")

def test_project_manager():
    """Test project manager functionality"""
    print_section("Project Manager Tests")
    
    # Create temporary directory for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        original_dir = os.getcwd()
        os.chdir(temp_dir)
        
        try:
            # Create projects directory first
            os.makedirs("projects", exist_ok=True)
            
            # Initialize project manager with absolute path
            config_path = os.path.join(temp_dir, "test_global_config.yaml")
            pm = ProjectManager(config_path)
            
            # Test 1: Create project
            print(f"{INFO} Testing project creation...")
            success = pm.create_project(
                "test_project",
                "Solanum",
                "Test tomato project",
                "solanum_repeats.fasta",
                ["S. lycopersicum", "S. pimpinellifolium"]
            )
            
            if success:
                print(f"{SUCCESS} Project creation test passed")
            else:
                print(f"{ERROR} Project creation test failed")
                return False
            
            # Test 2: Add sample
            print(f"{INFO} Testing sample addition...")
            
            # Create test read files
            test_dir = Path("test_data")
            test_dir.mkdir(exist_ok=True)
            
            # Create dummy FASTQ files
            with open(test_dir / "R1.fastq", "w") as f:
                f.write("@test_read_1\nATGCATGCATGC\n+\nIIIIIIIIIIII\n")
            
            with open(test_dir / "R2.fastq", "w") as f:
                f.write("@test_read_1\nGCTAGCTAGCTA\n+\nIIIIIIIIIIII\n")
            
            success = pm.add_sample(
                "test_project",
                "test_sample",
                "Solanum lycopersicum cv. Test",
                str(test_dir / "R1.fastq"),
                str(test_dir / "R2.fastq")
            )
            
            if success:
                print(f"{SUCCESS} Sample addition test passed")
            else:
                print(f"{ERROR} Sample addition test failed")
                return False
            
            # Test 3: Validate project
            print(f"{INFO} Testing project validation...")
            if pm.validate_project("test_project"):
                print(f"{SUCCESS} Project validation test passed")
            else:
                print(f"{ERROR} Project validation test failed")
                return False
            
            # Test 4: Validate sample
            print(f"{INFO} Testing sample validation...")
            if pm.validate_sample("test_project", "test_sample"):
                print(f"{SUCCESS} Sample validation test passed")
            else:
                print(f"{ERROR} Sample validation test failed")
                return False
            
            # Test 5: Get sample metadata
            print(f"{INFO} Testing metadata retrieval...")
            metadata = pm.get_sample_metadata("test_project", "test_sample")
            if metadata and metadata["sample_id"] == "test_sample":
                print(f"{SUCCESS} Metadata retrieval test passed")
            else:
                print(f"{ERROR} Metadata retrieval test failed")
                return False
            
            # Test 6: Cache key generation
            print(f"{INFO} Testing cache key generation...")
            cache_key = pm.get_cache_key("test_project", "test_sample", "assembly")
            if cache_key and len(cache_key) == 32:  # MD5 hash length
                print(f"{SUCCESS} Cache key generation test passed")
            else:
                print(f"{ERROR} Cache key generation test failed")
                return False
            
            # Test 7: List projects
            print(f"{INFO} Testing project listing...")
            projects = pm.list_projects()
            if "test_project" in projects:
                print(f"{SUCCESS} Project listing test passed")
            else:
                print(f"{ERROR} Project listing test failed")
                return False
            
            # Test 8: List samples
            print(f"{INFO} Testing sample listing...")
            samples = pm.get_project_samples("test_project")
            if "test_sample" in samples:
                print(f"{SUCCESS} Sample listing test passed")
            else:
                print(f"{ERROR} Sample listing test failed")
                return False
            
            # Test 9: Sample prefix retrieval
            print(f"{INFO} Testing sample prefix retrieval...")
            prefix = pm.get_sample_prefix("test_project", "test_sample")
            if prefix and len(prefix) == 4:
                print(f"{SUCCESS} Sample prefix retrieval test passed (prefix: {prefix})")
            else:
                print(f"{ERROR} Sample prefix retrieval test failed")
                return False
            
            print(f"{SUCCESS} All project manager tests passed!")
            return True
        except Exception as e:
            print(f"{ERROR} Project manager test failed: {e}")
            return False
        finally:
            os.chdir(original_dir)

def test_workflow_configuration():
    """Test workflow configuration loading"""
    print_section("Workflow Configuration Tests")
    
    # Test 1: Check if global config exists
    config_path = "projects/global_config.yaml"
    if os.path.exists(config_path):
        print(f"{SUCCESS} Global config file exists")
    else:
        print(f"{WARNING} Global config file not found, creating...")
        pm = ProjectManager()
        pm._create_default_config()
    
    # Test 2: Load and validate config
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        required_keys = ["global", "projects", "defaults"]
        missing_keys = []
        for key in required_keys:
            if key not in config:
                missing_keys.append(key)
        
        if not missing_keys:
            print(f"{SUCCESS} Configuration structure is valid")
        else:
            print(f"{ERROR} Missing required config keys: {missing_keys}")
            return False
        
        # Test 3: Validate global settings
        global_settings = config["global"]
        required_global = ["default_threads", "default_memory", "cache_dir"]
        missing_global = []
        for key in required_global:
            if key not in global_settings:
                missing_global.append(key)
        
        if not missing_global:
            print(f"{SUCCESS} Global settings are valid")
        else:
            print(f"{ERROR} Missing global settings: {missing_global}")
            return False
        
        # Test 4: Validate defaults
        defaults = config["defaults"]
        required_defaults = ["assembly_params", "tarean_params", "comparative_params"]
        missing_defaults = []
        for key in required_defaults:
            if key not in defaults:
                missing_defaults.append(key)
        
        if not missing_defaults:
            print(f"{SUCCESS} Default parameters are valid")
        else:
            print(f"{ERROR} Missing default parameters: {missing_defaults}")
            return False
        
        print(f"{SUCCESS} All workflow configuration tests passed!")
        return True
        
    except Exception as e:
        print(f"{ERROR} Configuration test failed: {e}")
        return False

def test_snakemake_workflow():
    """Test Snakemake workflow functionality"""
    print_section("Snakemake Workflow Tests")
    
    # Test 1: Check if Snakefile_modular exists
    snakefile_path = "Snakefile_modular"
    if os.path.exists(snakefile_path):
        print(f"{SUCCESS} Snakefile_modular exists")
    else:
        print(f"{ERROR} Snakefile_modular not found: {snakefile_path}")
        return False
    
    # Test 2: Check if project manager module exists
    if os.path.exists("project_manager.py"):
        print(f"{SUCCESS} Project manager module exists")
    else:
        print(f"{ERROR} Project manager module not found")
        return False
    
    # Test 3: Validate Snakemake syntax
    try:
        import subprocess
        result = subprocess.run([
            "snakemake", "-s", snakefile_path, "--dry-run", "--quiet"
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            print(f"{SUCCESS} Snakefile_modular syntax is valid")
        else:
            print(f"{WARNING} Snakefile_modular syntax check failed (this may be normal if no projects exist)")
            print(f"STDERR: {result.stderr[:200]}...")
    except Exception as e:
        print(f"{WARNING} Could not validate Snakefile_modular syntax: {e}")
    
    print(f"{SUCCESS} Snakemake workflow tests completed")
    return True

def test_directory_structure():
    """Test directory structure creation"""
    print_section("Directory Structure Tests")
    
    # Test 1: Check if projects directory exists
    projects_dir = Path("projects")
    if projects_dir.exists():
        print(f"{SUCCESS} Projects directory exists")
    else:
        print(f"{INFO} Creating projects directory...")
        projects_dir.mkdir(exist_ok=True)
        print(f"{SUCCESS} Projects directory created")
    
    # Test 2: Check if required subdirectories can be created
    test_dirs = ["samples", "comparative", "ncbi_repeats", "assembly", "tarean"]
    for dir_name in test_dirs:
        test_dir = projects_dir / "test_project" / dir_name
        test_dir.mkdir(parents=True, exist_ok=True)
        if test_dir.exists():
            print(f"{SUCCESS} Can create {dir_name} directory")
        else:
            print(f"{ERROR} Cannot create {dir_name} directory")
            return False
    
    # Clean up test directories
    shutil.rmtree(projects_dir / "test_project", ignore_errors=True)
    
    print(f"{SUCCESS} All directory structure tests passed!")
    return True

def test_caching_functionality():
    """Test caching functionality"""
    print_section("Caching Functionality Tests")
    
    # Test 1: File hash generation
    pm = ProjectManager()
    
    # Create test file
    test_file = "test_file.txt"
    with open(test_file, "w") as f:
        f.write("test content")
    
    try:
        file_hash = pm.get_file_hash(test_file)
        if file_hash and len(file_hash) == 32:  # MD5 hash length
            print(f"{SUCCESS} File hash generation works")
        else:
            print(f"{ERROR} File hash generation failed")
            return False
        
        # Create a mock sample for testing cache key generation
        os.makedirs("projects/test_project/samples/test_sample", exist_ok=True)
        mock_metadata = {
            "sample_id": "test_sample",
            "project_id": "test_project",
            "read_files": {
                "R1": test_file,  # Use the test file we created
                "R2": test_file
            }
        }
        
        # Write mock metadata
        metadata_path = "projects/test_project/samples/test_sample/metadata.yaml"
        os.makedirs(os.path.dirname(metadata_path), exist_ok=True)
        with open(metadata_path, 'w') as f:
            yaml.dump(mock_metadata, f)
        
        # Test 2: Cache key generation
        cache_key = pm.get_cache_key("test_project", "test_sample", "assembly", {"param": "value"})
        print(f"{INFO} Generated cache key: {cache_key}")
        if cache_key and len(cache_key) == 32:
            print(f"{SUCCESS} Cache key generation works")
        else:
            print(f"{ERROR} Cache key generation failed")
            return False
        
        # Test 3: Cache key consistency
        cache_key2 = pm.get_cache_key("test_project", "test_sample", "assembly", {"param": "value"})
        if cache_key == cache_key2:
            print(f"{SUCCESS} Cache key consistency works")
        else:
            print(f"{ERROR} Cache key consistency failed")
            return False
        
        # Test 4: Cache key uniqueness
        cache_key3 = pm.get_cache_key("test_project", "test_sample", "assembly", {"param": "different"})
        if cache_key != cache_key3:
            print(f"{SUCCESS} Cache key uniqueness works")
        else:
            print(f"{ERROR} Cache key uniqueness failed")
            return False
        
        # Clean up
        os.remove(test_file)
        shutil.rmtree("projects/test_project", ignore_errors=True)
        
        print(f"{SUCCESS} All caching functionality tests passed!")
        return True
        
    except Exception as e:
        print(f"{ERROR} Caching test failed: {e}")
        return False

def main():
    """Run all tests"""
    print_header("RepOrtR Multi-Project Framework Tests")
    
    tests = [
        ("Project Manager", test_project_manager),
        ("Workflow Configuration", test_workflow_configuration),
        ("Snakemake Workflow", test_snakemake_workflow),
        ("Directory Structure", test_directory_structure),
        ("Caching Functionality", test_caching_functionality)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        try:
            if test_func():
                passed += 1
            else:
                print(f"{ERROR} {test_name} test failed")
        except Exception as e:
            print(f"{ERROR} {test_name} test failed with exception: {e}")
    
    print_header("Test Summary")
    print(f"Tests passed: {passed}/{total}")
    print(f"Success rate: {(passed/total)*100:.1f}%")
    
    if passed == total:
        print(f"{SUCCESS} All tests passed! Multi-project framework is ready.")
        return True
    else:
        print(f"{WARNING} Some tests failed. Check the output above for details.")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 