#!/usr/bin/env python3
"""
RepOrtR Modular Snakemake Workflow Runner
This script runs the modular Snakemake workflow with configuration validation
"""

import argparse
import subprocess
import sys
import os
import yaml
from pathlib import Path
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_config(config_file="projects/global_config.yaml"):
    """Load configuration from YAML file"""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def validate_configuration(config):
    """Validate that configuration has required keys"""
    required_global_keys = ["default_threads", "default_memory", "pythonhashseed", "cache_dir"]
    required_read_cleaning_keys = ["adapters", "bbduk_params", "fastqc_threads"]
    required_read_preparation_keys = ["temp_dir", "random_n"]
    
    # Check global configuration
    for key in required_global_keys:
        if key not in config.get("global", {}):
            raise ValueError(f"Missing required global config key: {key}")
    
    # Check read cleaning configuration
    for key in required_read_cleaning_keys:
        if key not in config.get("global", {}).get("read_cleaning", {}):
            raise ValueError(f"Missing required read_cleaning config key: {key}")
    
    # Check read preparation configuration
    for key in required_read_preparation_keys:
        if key not in config.get("global", {}).get("read_preparation", {}):
            raise ValueError(f"Missing required read_preparation config key: {key}")
    
    # Check project configurations
    for project_id, project_config in config.get("projects", {}).items():
        if project_id in ["defaults", "comparative_params", "post_tarean_params"]:
            continue  # Skip special config sections
        
        required_keys = ["samples", "tarean_params", "comparative_species"]
        for key in required_keys:
            if key not in project_config:
                raise ValueError(f"Missing required key '{key}' in project '{project_id}'")
        
        # Check TAREAN params
        tarean_keys = ["assembly_min", "mincl", "threads"]
        for key in tarean_keys:
            if key not in project_config.get("tarean_params", {}):
                raise ValueError(f"Missing required TAREAN param '{key}' in project '{project_id}'")
    
    logger.info("Configuration validation passed")
    return True

def run_snakemake_modular(config_file, snakefile="Snakefile_modular", dry_run=False, cores=None):
    """Run modular Snakemake workflow with validation"""
    cmd = ["snakemake", "-s", snakefile, "--configfile", config_file]
    
    if dry_run:
        cmd.append("--dry-run")
    
    if cores:
        cmd.extend(["--cores", str(cores)])
    else:
        cmd.extend(["--cores", str(1)])
    
    # Add other useful Snakemake options
    cmd.extend([
        "--rerun-incomplete",
        "--keep-going",
        "--printshellcmds"
    ])
    
    logger.info(f"Running modular Snakemake command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Modular Snakemake workflow completed successfully")
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Modular Snakemake workflow failed: {e}")
        logger.error(f"STDOUT: {e.stdout}")
        logger.error(f"STDERR: {e.stderr}")
        print(e.stdout)
        print(e.stderr, file=sys.stderr)
        return False

def run_validation_only(config_file, snakefile="Snakefile_modular"):
    """Run only validation rules"""
    cmd = ["snakemake", "-s", snakefile, "--configfile", config_file, 
           "--cores", "1", "validate_workflow"]
    
    logger.info(f"Running validation: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Validation completed successfully")
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Validation failed: {e}")
        print(e.stdout)
        print(e.stderr, file=sys.stderr)
        return False

def run_single_project(project_id, config_file):
    """Run workflow for a single project"""
    logger.info(f"Processing project: {project_id}")
    
    # Create project-specific config
    project_config = load_config(config_file)
    
    # Validate configuration
    try:
        validate_configuration(project_config)
    except ValueError as e:
        logger.error(f"Configuration validation failed: {e}")
        return False
    
    # Run modular Snakemake for this project
    success = run_snakemake_modular(config_file, cores=1)
    
    return success

def run_parallel_projects(projects, config_file, max_workers):
    """Run workflow for multiple projects in parallel"""
    logger.info(f"Running parallel processing for {len(projects)} projects with {max_workers} workers")
    
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_single_project, project, config_file) for project in projects]
        
        results = []
        for future in futures:
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                logger.error(f"Project processing failed: {e}")
                results.append(False)
    
    return results

def check_modular_structure():
    """Check that all required modular files exist"""
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
        logger.error(f"Missing required modular files: {missing_files}")
        return False
    
    logger.info("All modular workflow files present")
    return True

def main():
    parser = argparse.ArgumentParser(description="RepOrtR Modular Snakemake Workflow Runner")
    parser.add_argument('--config', default='projects/global_config.yaml',
                       help='Configuration file')
    parser.add_argument('--snakefile', default='Snakefile_modular',
                       help='Snakemake file to use')
    parser.add_argument('--dry-run', action='store_true',
                       help='Dry run - show what would be executed')
    parser.add_argument('--cores', type=int, default=1,
                       help='Number of cores to use')
    parser.add_argument('--validate-only', action='store_true',
                       help='Run only validation rules')
    parser.add_argument('--projects', nargs='+',
                       help='Specific projects to process (default: all)')
    
    args = parser.parse_args()
    
    # Check modular structure
    if not check_modular_structure():
        sys.exit(1)
    
    # Load configuration
    try:
        config = load_config(args.config)
    except FileNotFoundError:
        logger.error(f"Configuration file {args.config} not found")
        sys.exit(1)
    
    # Validate configuration
    try:
        validate_configuration(config)
    except ValueError as e:
        logger.error(f"Configuration validation failed: {e}")
        sys.exit(1)
    
    # Determine projects to process
    if args.projects:
        projects_to_process = args.projects
    else:
        # Get all projects except special config sections
        projects_to_process = [p for p in config["projects"].keys() 
                             if p not in ["defaults", "comparative_params", "post_tarean_params"]]
    
    logger.info(f"Processing {len(projects_to_process)} projects: {projects_to_process}")
    
    if args.validate_only:
        # Run only validation
        success = run_validation_only(args.config, args.snakefile)
    else:
        # Run full workflow
        if len(projects_to_process) == 1:
            # Single project
            success = run_single_project(projects_to_process[0], args.config)
        else:
            # Multiple projects - run in parallel
            max_workers = min(args.cores, len(projects_to_process))
            results = run_parallel_projects(projects_to_process, args.config, max_workers)
            success = all(results)
    
    if success:
        logger.info("All projects processed successfully")
        sys.exit(0)
    else:
        logger.error("Some projects failed to process")
        sys.exit(1)

if __name__ == '__main__':
    main() 