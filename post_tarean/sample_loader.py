#!/usr/bin/env python3
"""
Dynamic Sample Loader for RepOrtR

This module provides dynamic sample loading functionality to replace hardcoded test samples.

Author: RepOrtR Team
Date: 2025
"""

import os
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass

@dataclass
class SampleConfig:
    """Represents a sample configuration"""
    sample_id: str
    genome_size: float
    r1_path: str
    r2_path: str
    metadata: Dict = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

class DynamicSampleLoader:
    """
    Dynamic sample loader that replaces hardcoded test samples
    """
    
    def __init__(self, project_id: str):
        """
        Initialize the sample loader
        
        Args:
            project_id: Project identifier
        """
        self.project_id = project_id
        self.config_path = f"projects/{project_id}/project_metadata.yaml"
        self.global_config_path = "projects/global_config.yaml"
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def load_project_samples(self) -> Dict[str, SampleConfig]:
        """
        Load samples dynamically from project configuration
        
        Returns:
            Dictionary of sample configurations
        """
        self.logger.info(f"Loading samples for project: {self.project_id}")
        
        # Try to load from project-specific config first
        if os.path.exists(self.config_path):
            return self._load_from_project_config()
        
        # Fall back to global config
        return self._load_from_global_config()
    
    def _load_from_project_config(self) -> Dict[str, SampleConfig]:
        """Load samples from project-specific configuration"""
        try:
            with open(self.config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            samples = {}
            for sample_id, sample_data in config.get('samples', {}).items():
                samples[sample_id] = SampleConfig(
                    sample_id=sample_id,
                    genome_size=sample_data.get('genome_size', 1.0),
                    r1_path=sample_data.get('r1_path', ''),
                    r2_path=sample_data.get('r2_path', ''),
                    metadata=sample_data.get('metadata', {})
                )
            
            self.logger.info(f"Loaded {len(samples)} samples from project config")
            return samples
            
        except Exception as e:
            self.logger.error(f"Failed to load project config: {e}")
            return {}
    
    def _load_from_global_config(self) -> Dict[str, SampleConfig]:
        """Load samples from global configuration"""
        try:
            with open(self.global_config_path, 'r') as f:
                config = yaml.safe_load(f)
            
            project_config = config.get('projects', {}).get(self.project_id, {})
            samples = {}
            
            for sample_id, sample_data in project_config.get('samples', {}).items():
                samples[sample_id] = SampleConfig(
                    sample_id=sample_id,
                    genome_size=sample_data.get('genome_size', 1.0),
                    r1_path=sample_data.get('r1_path', ''),
                    r2_path=sample_data.get('r2_path', ''),
                    metadata=sample_data.get('metadata', {})
                )
            
            self.logger.info(f"Loaded {len(samples)} samples from global config")
            return samples
            
        except Exception as e:
            self.logger.error(f"Failed to load global config: {e}")
            return {}
    
    def validate_samples(self, samples: Dict[str, SampleConfig]) -> List[str]:
        """
        Validate that samples have required files
        
        Args:
            samples: Dictionary of sample configurations
            
        Returns:
            List of validation errors
        """
        errors = []
        
        for sample_id, sample_config in samples.items():
            # Check if R1 file exists
            if not os.path.exists(sample_config.r1_path):
                errors.append(f"Sample {sample_id}: R1 file not found: {sample_config.r1_path}")
            
            # Check if R2 file exists
            if not os.path.exists(sample_config.r2_path):
                errors.append(f"Sample {sample_id}: R2 file not found: {sample_config.r2_path}")
            
            # Check if genome size is reasonable
            if sample_config.genome_size <= 0 or sample_config.genome_size > 100:
                errors.append(f"Sample {sample_id}: Invalid genome size: {sample_config.genome_size}")
        
        return errors
    
    def get_sample_list(self) -> List[str]:
        """Get list of sample IDs"""
        samples = self.load_project_samples()
        return list(samples.keys())
    
    def get_sample_paths(self) -> Dict[str, Dict[str, str]]:
        """Get sample paths for Snakemake integration"""
        samples = self.load_project_samples()
        paths = {}
        
        for sample_id, sample_config in samples.items():
            paths[sample_id] = {
                'r1_path': sample_config.r1_path,
                'r2_path': sample_config.r2_path,
                'genome_size': str(sample_config.genome_size)
            }
        
        return paths

def load_samples_for_project(project_id: str) -> Dict[str, SampleConfig]:
    """
    Convenience function to load samples for a project
    
    Args:
        project_id: Project identifier
        
    Returns:
        Dictionary of sample configurations
    """
    loader = DynamicSampleLoader(project_id)
    return loader.load_project_samples()

def validate_project_samples(project_id: str) -> bool:
    """
    Validate all samples for a project
    
    Args:
        project_id: Project identifier
        
    Returns:
        True if all samples are valid, False otherwise
    """
    loader = DynamicSampleLoader(project_id)
    samples = loader.load_project_samples()
    errors = loader.validate_samples(samples)
    
    if errors:
        for error in errors:
            print(f"❌ {error}")
        return False
    
    print(f"✅ All {len(samples)} samples validated successfully")
    return True

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python sample_loader.py <project_id>")
        sys.exit(1)
    
    project_id = sys.argv[1]
    
    # Load and validate samples
    samples = load_samples_for_project(project_id)
    
    print(f"\nSamples for project '{project_id}':")
    for sample_id, sample_config in samples.items():
        print(f"  {sample_id}:")
        print(f"    Genome size: {sample_config.genome_size}")
        print(f"    R1: {sample_config.r1_path}")
        print(f"    R2: {sample_config.r2_path}")
    
    # Validate samples
    if validate_project_samples(project_id):
        print("\n✅ All samples are valid")
    else:
        print("\n❌ Some samples have issues")
        sys.exit(1)

