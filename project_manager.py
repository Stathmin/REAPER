#!/usr/bin/env python3
"""
REAPER project manager
Handles multi-project configuration and sample management

HOLY WORKFLOW COMPLIANCE: This module manages project configuration for the modular workflow.
All project setup is automated - users never create directories manually.
All Snakemake rules use the reportr environment exclusively.
"""

import os
import yaml  # type: ignore[import-not-found]
import shutil
from pathlib import Path
from typing import Dict, List, Optional
import hashlib
import json

class ProjectManager:
    """Manages multiple projects and their configurations"""
    
    def __init__(self, global_config_path: str = "projects/global_config.yaml"):
        self.global_config_path = global_config_path
        self.config = self._load_global_config()
        self.projects_dir = Path("projects")
        self.projects_dir.mkdir(exist_ok=True)
    
    def _load_global_config(self) -> Dict:
        """Load global configuration"""
        if os.path.exists(self.global_config_path):
            with open(self.global_config_path, 'r') as f:
                return yaml.safe_load(f)
        else:
            return self._create_default_config()
    
    def _create_default_config(self) -> Dict:
        """Create default global configuration"""
        config = {
            "global": {
                "default_threads": 4,
                "default_memory": "8G",
                "cache_dir": ".snakemake/cache",
                "log_dir": "logs",
                "temp_dir": "tmp"
            },
            "projects": {},
            "defaults": {
                "assembly_params": {
                    "kmer_sizes": [21, 33, 55],
                    "min_contig_length": 200,
                    "threads": 4
                },
                "tarean_params": {
                    "assembly_min": 4,
                    "mincl": 0.001,
                    "r_value": 178067846,
                    "threads": 4
                },
                "comparative_params": {
                    "blast_evalue": 1e-5,
                    "min_identity": 80,
                    "min_coverage": 80
                }
            }
        }
        self._save_global_config(config)
        return config
    
    def _save_global_config(self, config: Dict):
        """Save global configuration"""
        os.makedirs(os.path.dirname(self.global_config_path), exist_ok=True)
        with open(self.global_config_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
    
    def create_project(
        self,
        project_id: str,
        taxonomy: str,
        description: str,
        ncbi_repeats: str,
        comparative_species: List[str],
        total_reads_per_assembly: Optional[int] = None,
        taxonomy_id: Optional[str] = None,
        tarean_overrides: Optional[Dict] = None,
    ) -> bool:
        """Create a new project"""
        try:
            # Create project directory
            project_dir = self.projects_dir / project_id
            project_dir.mkdir(exist_ok=True)
            
            # Create subdirectories
            (project_dir / "samples").mkdir(exist_ok=True)
            (project_dir / "comparative").mkdir(exist_ok=True)
            (project_dir / "ncbi_repeats").mkdir(exist_ok=True)
            
            # Create project metadata
            project_metadata = {
                "project_id": project_id,
                "taxonomy": taxonomy,
                "description": description,
                "ncbi_repeats": ncbi_repeats,
                # Legacy field: comparative analyses should be configured via
                # projects.<id>.comparative_analyses (or `add-comparative` CLI).
                "comparative_species": comparative_species,
                "total_reads_per_assembly": total_reads_per_assembly,
                "created_date": self._get_timestamp(),
                "samples": {},
            }
            if taxonomy_id:
                project_metadata["taxonomy_id"] = taxonomy_id
            
            with open(project_dir / "project_metadata.yaml", 'w') as f:
                yaml.dump(project_metadata, f, default_flow_style=False, sort_keys=False)
            
            # Update global config
            tarean_params = self.config["defaults"]["tarean_params"].copy()
            if tarean_overrides:
                # Only update known keys; keep other defaults intact.
                tarean_params.update({k: v for k, v in tarean_overrides.items() if v is not None})

            self.config["projects"][project_id] = {
                "taxonomy": taxonomy,
                "description": description,
                "ncbi_repeats": ncbi_repeats,
                "tarean_params": tarean_params,
                "samples": {},
                "comparative_analyses": {},
            }
            if comparative_species:
                # Back-compat shorthand: treat `--comparative-species` as a request
                # to create a default comparative analysis. New code should use
                # `add-comparative` / `comparative_analyses` directly.
                self.config["projects"][project_id]["comparative_analyses"]["default"] = {
                    "samples": comparative_species,
                    "description": f"Default comparative analysis of {', '.join(comparative_species)}",
                }
                # Keep legacy key for older configs/tools that still read it.
                self.config["projects"][project_id]["comparative_species"] = comparative_species
            if taxonomy_id:
                # Store both a single taxonomy_id and a list taxonomy_ids
                # for compatibility with ncbi_gathering_rules.smk.
                self.config["projects"][project_id]["taxonomy_id"] = taxonomy_id
                self.config["projects"][project_id]["taxonomy_ids"] = [taxonomy_id]
            if total_reads_per_assembly is not None:
                self.config["projects"][project_id]["total_reads_per_assembly"] = total_reads_per_assembly
            
            self._save_global_config(self.config)
            print(f"✅ Project '{project_id}' created successfully")
            return True
            
        except Exception as e:
            print(f"❌ Failed to create project '{project_id}': {e}")
            return False
    
    def _generate_unique_prefix(self, project_id: str, sample_id: str) -> str:
        """Generate a unique prefix for a sample (no fixed slicing)."""
        # Get existing prefixes in this project
        existing_prefixes = set()
        if project_id in self.config["projects"]:
            samples = self.config["projects"][project_id].get("samples", {})
            for existing_sample_id, sample_config in samples.items():
                if "prefix" in sample_config:
                    existing_prefixes.add(sample_config["prefix"])
        
        # Primary scheme: full sample_id uppercased (no truncation).
        base = str(sample_id).upper()
        if base not in existing_prefixes:
            return base

        # If user already has a colliding prefix, deterministically disambiguate.
        suffix = hashlib.md5(str(sample_id).encode("utf-8")).hexdigest()[:6].upper()
        candidate = f"{base}_{suffix}"
        if candidate not in existing_prefixes:
            return candidate

        # Last resort: keep extending the hash (extremely unlikely to be needed).
        full = hashlib.md5(str(sample_id).encode("utf-8")).hexdigest().upper()
        for n in range(8, len(full) + 1, 2):
            candidate = f"{base}_{full[:n]}"
            if candidate not in existing_prefixes:
                return candidate
        return f"{base}_{full}"

    def add_sample(self, project_id: str, sample_id: str, taxonomy: str,
                   r1_path: str, r2_path: str, genome_size: float = 1.0) -> bool:
        """Add a sample to a project with a unique seqclust read prefix."""
        try:
            if project_id not in self.config["projects"]:
                print(f"❌ Project '{project_id}' not found")
                return False
            
            # Generate unique prefix (uppercased sample_id, disambiguated on collision)
            prefix = self._generate_unique_prefix(project_id, sample_id)
            
            # Create sample directory with correct structure
            sample_dir = self.projects_dir / project_id / "samples" / sample_id
            sample_dir.mkdir(exist_ok=True)
            raw_reads_dir = sample_dir / "raw_reads"
            raw_reads_dir.mkdir(exist_ok=True)
            (sample_dir / "filtered_reads").mkdir(exist_ok=True)
            (sample_dir / "fastqc").mkdir(exist_ok=True)
            (sample_dir / "tarean").mkdir(exist_ok=True)
            (sample_dir / "post_tarean").mkdir(exist_ok=True)
            (sample_dir / "test").mkdir(exist_ok=True)
            
            # Create symbolic links to read files
            # Use .fq.gz links when the underlying files are gzipped so tools
            # can auto-detect compression based on extension.
            r1_ext = ".fq.gz" if r1_path.endswith(".gz") else ".fq"
            r2_ext = ".fq.gz" if r2_path.endswith(".gz") else ".fq"
            r1_link = raw_reads_dir / f"R1{r1_ext}"
            r2_link = raw_reads_dir / f"R2{r2_ext}"
            
            # Create absolute symbolic links.
            # Relative links produce long ../../.. paths depending on where the repo
            # lives vs the input files (e.g. /storage/...), and they become harder
            # to debug. Absolute targets are clearer and more robust across moves
            # of the project directory.
            r1_abs_path = str(Path(r1_path).expanduser().resolve())
            r2_abs_path = str(Path(r2_path).expanduser().resolve())
            
            if r1_link.exists():
                r1_link.unlink()
            if r2_link.exists():
                r2_link.unlink()
            
            r1_link.symlink_to(r1_abs_path)
            r2_link.symlink_to(r2_abs_path)
            
            # Create sample metadata (tarean_status only; no assembly_status)
            sample_metadata = {
                "sample_id": sample_id,
                "project_id": project_id,
                "taxonomy": taxonomy,
                "genome_size": genome_size,
                "prefix": prefix,  # Store the unique prefix
                "read_files": {
                    "R1": str(r1_link),
                    "R2": str(r2_link)
                },
                "tarean_status": "pending",
                "added_date": self._get_timestamp()
            }
            
            with open(sample_dir / "metadata.yaml", 'w') as f:
                yaml.dump(sample_metadata, f, default_flow_style=False, sort_keys=False)
            
            # Update global config with new sample structure
            if "samples" not in self.config["projects"][project_id]:
                self.config["projects"][project_id]["samples"] = {}
            
            self.config["projects"][project_id]["samples"][sample_id] = {
                "genome_size": genome_size,
                "r1_path": r1_path,
                "r2_path": r2_path,
                "prefix": prefix  # Store the unique prefix in global config
            }
            
            self._save_global_config(self.config)
            
            print(f"✅ Sample '{sample_id}' added to project '{project_id}' (genome_size: {genome_size}, prefix: {prefix})")
            print(f"   R1 link: {r1_link} -> {r1_abs_path}")
            print(f"   R2 link: {r2_link} -> {r2_abs_path}")
            return True
            
        except Exception as e:
            print(f"❌ Failed to add sample '{sample_id}': {e}")
            return False
    
    def get_project_samples(self, project_id: str) -> List[str]:
        """Get list of samples in a project"""
        if project_id in self.config["projects"]:
            samples = self.config["projects"][project_id].get("samples", {})
            if isinstance(samples, dict):
                return list(samples.keys())
            elif isinstance(samples, list):
                return samples
        return []
    
    def get_project_config(self, project_id: str) -> Optional[Dict]:
        """Get project configuration"""
        return self.config["projects"].get(project_id)
    
    def list_projects(self) -> List[str]:
        """List all projects"""
        return list(self.config["projects"].keys())
    
    def get_sample_metadata(self, project_id: str, sample_id: str) -> Optional[Dict]:
        """Get sample metadata"""
        metadata_path = self.projects_dir / project_id / "samples" / sample_id / "metadata.yaml"
        if metadata_path.exists():
            with open(metadata_path, 'r') as f:
                return yaml.safe_load(f)
        return None
    
    def update_sample_tarean_status(self, project_id: str, sample_id: str, tarean_status: str):
        """Update sample metadata.yaml with current tarean_status."""
        metadata = self.get_sample_metadata(project_id, sample_id)
        if metadata:
            metadata["tarean_status"] = tarean_status
            metadata_path = self.projects_dir / project_id / "samples" / sample_id / "metadata.yaml"
            with open(metadata_path, 'w') as f:
                yaml.dump(metadata, f, default_flow_style=False, sort_keys=False)
    
    def get_file_hash(self, file_path: str) -> str:
        """Calculate file hash for caching"""
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def get_cache_key(self, project_id: str, sample_id: str, step: str, 
                     params: Dict = None) -> str:
        """Generate cache key for a processing step"""
        # Get input file hashes
        sample_metadata = self.get_sample_metadata(project_id, sample_id)
        if not sample_metadata:
            return None
        
        try:
            r1_hash = self.get_file_hash(sample_metadata["read_files"]["R1"])
            r2_hash = self.get_file_hash(sample_metadata["read_files"]["R2"])
        except Exception as e:
            # If files don't exist, use placeholder hashes
            r1_hash = "placeholder_r1_hash"
            r2_hash = "placeholder_r2_hash"
        
        # Combine hashes with parameters
        cache_data = {
            "project_id": project_id,
            "sample_id": sample_id,
            "step": step,
            "r1_hash": r1_hash,
            "r2_hash": r2_hash,
            "params": params or {}
        }
        
        return hashlib.md5(json.dumps(cache_data, sort_keys=True).encode()).hexdigest()
    
    def _get_timestamp(self) -> str:
        """Get current timestamp"""
        from datetime import datetime
        return datetime.now().isoformat()
    
    def validate_project(self, project_id: str) -> bool:
        """Validate project structure and basic config/filesystem consistency"""
        project_dir = self.projects_dir / project_id
        if not project_dir.exists():
            print(f"❌ Project directory not found: {project_dir}")
            return False
        
        required_dirs = ["samples", "comparative", "ncbi_repeats"]
        for dir_name in required_dirs:
            if not (project_dir / dir_name).exists():
                print(f"❌ Required directory missing: {project_dir / dir_name}")
                return False

        # Optional: cross-check config vs filesystem
        project_config = self.get_project_config(project_id)
        if project_config:
            config_samples = project_config.get("samples", {}) or {}
            config_sample_ids = set(config_samples.keys())

            samples_dir = project_dir / "samples"
            fs_sample_ids = set()
            if samples_dir.exists():
                for entry in samples_dir.iterdir():
                    if entry.is_dir():
                        fs_sample_ids.add(entry.name)

            # Warn about mismatches without failing the whole validation
            missing_on_disk = config_sample_ids - fs_sample_ids
            missing_in_config = fs_sample_ids - config_sample_ids
            if missing_on_disk:
                print(f"⚠️ Samples in config but missing on disk: {', '.join(sorted(missing_on_disk))}")
            if missing_in_config:
                print(f"⚠️ Sample directories missing from config: {', '.join(sorted(missing_in_config))}")

        return True
    
    def validate_sample(self, project_id: str, sample_id: str) -> bool:
        """Validate sample structure and raw read links"""
        sample_dir = self.projects_dir / project_id / "samples" / sample_id
        if not sample_dir.exists():
            print(f"❌ Sample directory not found: {sample_dir}")
            return False
        
        # Check read files (support fq/fastq, plain and gzipped symlinks)
        raw_reads_dir = sample_dir / "raw_reads"
        r1_candidates = [
            raw_reads_dir / "R1.fq.gz",
            raw_reads_dir / "R1.fastq.gz",
            raw_reads_dir / "R1.fq",
            raw_reads_dir / "R1.fastq",
        ]
        r2_candidates = [
            raw_reads_dir / "R2.fq.gz",
            raw_reads_dir / "R2.fastq.gz",
            raw_reads_dir / "R2.fq",
            raw_reads_dir / "R2.fastq",
        ]

        # `Path.exists()` follows symlinks; a broken symlink returns False.
        # For UX, treat "broken symlink" differently from "missing file".
        r1_path = next((p for p in r1_candidates if p.exists() or p.is_symlink()), None)
        r2_path = next((p for p in r2_candidates if p.exists() or p.is_symlink()), None)

        if r1_path is None:
            print(f"❌ R1 file not found in {raw_reads_dir} (expected R1.fq[.gz] or R1.fastq[.gz])")
            return False

        if r2_path is None:
            print(f"❌ R2 file not found in {raw_reads_dir} (expected R2.fq[.gz] or R2.fastq[.gz])")
            return False
        
        # Check that symlinks (if used) are not broken
        for label, path in [("R1", r1_path), ("R2", r2_path)]:
            if path.is_symlink():
                target = os.readlink(path)
                target_abs = (path.parent / target).resolve()
                if not target_abs.exists():
                    print(f"❌ {label} symlink is broken: {path} -> {target}")
                    return False

            # Non-symlink but missing target (e.g. file deleted) should also fail.
            if not path.exists():
                print(f"❌ {label} file missing: {path}")
                return False

        # Warn if metadata is missing
        metadata = self.get_sample_metadata(project_id, sample_id)
        if not metadata:
            print(f"⚠️ metadata.yaml not found for sample '{sample_id}' in project '{project_id}'")
        
        return True

    def validate_comparative_analysis(self, project_id: str, analysis_id: str) -> bool:
        """Validate comparative analysis config and directory layout (no assemblies required)."""
        project_cfg = self.get_project_config(project_id) or {}
        analyses = project_cfg.get("comparative_analyses", {}) or {}
        analysis = analyses.get(analysis_id) or {}

        samples = analysis.get("samples") or []
        if not isinstance(samples, list) or not samples:
            print(f"❌ Comparative '{analysis_id}' has no samples configured")
            return False

        # Referenced samples must exist in project config.
        project_samples = set((project_cfg.get("samples", {}) or {}).keys())
        missing = [s for s in samples if s not in project_samples]
        if missing:
            print(f"❌ Comparative '{analysis_id}' references missing samples: {', '.join(missing)}")
            return False

        comp_dir = self.projects_dir / project_id / "comparative" / analysis_id
        if not comp_dir.exists():
            print(f"❌ Comparative directory not found: {comp_dir}")
            return False

        return True

    def drop_project(self, project_id: str, delete_files: bool = True) -> bool:
        """Drop a project from configuration, optionally removing its directory."""
        try:
            if project_id not in self.config["projects"]:
                print(f"❌ Project '{project_id}' not found")
                return False
            
            project_dir = self.projects_dir / project_id
            if delete_files and project_dir.exists():
                shutil.rmtree(project_dir)
            
            del self.config["projects"][project_id]
            self._save_global_config(self.config)
            
            print(f"✅ Project '{project_id}' dropped successfully")
            if not delete_files:
                print(f"   Note: project directory '{project_dir}' was kept on disk")
            return True
        
        except Exception as e:
            print(f"❌ Failed to drop project '{project_id}': {e}")
            return False

    def add_comparative_analysis(self, project_id: str, analysis_id: str, 
                                sample_names: List[str], description: str = "") -> bool:
        """Add a user-defined comparative analysis"""
        try:
            if project_id not in self.config["projects"]:
                print(f"❌ Project '{project_id}' not found")
                return False
            
            # Validate that all samples exist
            project_samples = self.get_project_samples(project_id)
            for sample_name in sample_names:
                if sample_name not in project_samples:
                    print(f"❌ Sample '{sample_name}' not found in project '{project_id}'")
                    return False
            
            # Initialize comparative_analyses if it doesn't exist
            if "comparative_analyses" not in self.config["projects"][project_id]:
                self.config["projects"][project_id]["comparative_analyses"] = {}
            
            # Add the comparative analysis
            self.config["projects"][project_id]["comparative_analyses"][analysis_id] = {
                "samples": sample_names,
                "description": description or f"Comparative analysis of {', '.join(sample_names)}"
            }
            
            # Create comparative analysis directory
            comparative_dir = self.projects_dir / project_id / "comparative" / analysis_id
            comparative_dir.mkdir(parents=True, exist_ok=True)
            
            self._save_global_config(self.config)
            print(f"✅ Comparative analysis '{analysis_id}' added to project '{project_id}'")
            print(f"   Samples: {', '.join(sample_names)}")
            print(f"   Description: {self.config['projects'][project_id]['comparative_analyses'][analysis_id]['description']}")
            return True
            
        except Exception as e:
            print(f"❌ Error adding comparative analysis: {e}")
            return False
    
    def list_comparative_analyses(self, project_id: str) -> List[str]:
        """List all comparative analyses for a project"""
        if project_id not in self.config["projects"]:
            return []
        
        project_config = self.config["projects"][project_id]
        return list(project_config.get("comparative_analyses", {}).keys())
    
    def get_comparative_analysis(self, project_id: str, analysis_id: str) -> Optional[Dict]:
        """Get details of a specific comparative analysis"""
        if project_id not in self.config["projects"]:
            return None
        
        project_config = self.config["projects"][project_id]
        return project_config.get("comparative_analyses", {}).get(analysis_id)
    
    def remove_comparative_analysis(self, project_id: str, analysis_id: str) -> bool:
        """Remove a comparative analysis"""
        try:
            if project_id not in self.config["projects"]:
                print(f"❌ Project '{project_id}' not found")
                return False
            
            project_config = self.config["projects"][project_id]
            if "comparative_analyses" not in project_config:
                print(f"❌ No comparative analyses found in project '{project_id}'")
                return False
            
            if analysis_id not in project_config["comparative_analyses"]:
                print(f"❌ Comparative analysis '{analysis_id}' not found in project '{project_id}'")
                return False
            
            # Remove from config
            del project_config["comparative_analyses"][analysis_id]
            
            # Remove directory if it exists
            comparative_dir = self.projects_dir / project_id / "comparative" / analysis_id
            if comparative_dir.exists():
                shutil.rmtree(comparative_dir)
            
            self._save_global_config(self.config)
            print(f"✅ Comparative analysis '{analysis_id}' removed from project '{project_id}'")
            return True
            
        except Exception as e:
            print(f"❌ Error removing comparative analysis: {e}")
            return False
    
    def update_comparative_analysis(self, project_id: str, analysis_id: str, 
                                  sample_names: List[str] = None, description: str = None) -> bool:
        """Update an existing comparative analysis"""
        try:
            if project_id not in self.config["projects"]:
                print(f"❌ Project '{project_id}' not found")
                return False
            
            project_config = self.config["projects"][project_id]
            if "comparative_analyses" not in project_config or analysis_id not in project_config["comparative_analyses"]:
                print(f"❌ Comparative analysis '{analysis_id}' not found in project '{project_id}'")
                return False
            
            analysis_config = project_config["comparative_analyses"][analysis_id]
            
            # Update samples if provided
            if sample_names is not None:
                # Validate that all samples exist
                project_samples = self.get_project_samples(project_id)
                for sample_name in sample_names:
                    if sample_name not in project_samples:
                        print(f"❌ Sample '{sample_name}' not found in project '{project_id}'")
                        return False
                analysis_config["samples"] = sample_names
            
            # Update description if provided
            if description is not None:
                analysis_config["description"] = description
            
            self._save_global_config(self.config)
            print(f"✅ Comparative analysis '{analysis_id}' updated in project '{project_id}'")
            return True
            
        except Exception as e:
            print(f"❌ Error updating comparative analysis: {e}")
            return False

    def get_sample_prefix(self, project_id: str, sample_id: str) -> Optional[str]:
        """Get the assigned prefix for a sample"""
        if project_id in self.config["projects"]:
            samples = self.config["projects"][project_id].get("samples", {})
            if sample_id in samples:
                return samples[sample_id].get("prefix")
        return None
    
    def update_existing_samples_with_prefixes(self, project_id: str) -> bool:
        """Update existing samples in a project to have unique prefixes"""
        try:
            if project_id not in self.config["projects"]:
                print(f"❌ Project '{project_id}' not found")
                return False
            
            samples = self.config["projects"][project_id].get("samples", {})
            updated = False
            
            for sample_id in samples:
                if "prefix" not in samples[sample_id]:
                    # Generate unique prefix for this sample
                    prefix = self._generate_unique_prefix(project_id, sample_id)
                    samples[sample_id]["prefix"] = prefix
                    updated = True
                    print(f"✅ Added prefix '{prefix}' to sample '{sample_id}'")
            
            if updated:
                self._save_global_config(self.config)
                print(f"✅ Updated {project_id} with prefixes for {len(samples)} samples")
            
            return True
            
        except Exception as e:
            print(f"❌ Failed to update prefixes for project '{project_id}': {e}")
            return False

def main():
    """Command-line interface for project management
    
    HOLY WORKFLOW COMPLIANCE: This CLI manages project setup for the modular workflow.
    All directory creation is automated - users never create directories manually.
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="REAPER project manager")
    parser.add_argument(
        "action",
        choices=[
            "create-project",
            "update-project",
            "add-sample",
            "list-projects",
            "list-samples",
            "show-project",
            "show-sample",
            "validate",
            "add-comparative",
            "list-comparatives",
            "remove-comparative",
            "update-comparative",
            "update-prefixes",
            "drop-project",
        ],
    )
    parser.add_argument("--project-id", help="Project ID")
    parser.add_argument("--taxonomy", help="Taxonomy")
    parser.add_argument(
        "--taxonomy-id",
        help="Top-level NCBI taxonomy ID for this project (enables automatic NCBI gathering)",
    )
    parser.add_argument("--description", help="Project description")
    parser.add_argument("--ncbi-repeats", help="NCBI repeats file")
    parser.add_argument("--comparative-species", nargs="+", help="Comparative species")
    parser.add_argument("--sample-id", help="Sample ID")
    parser.add_argument("--r1-path", help="R1 FASTQ file path")
    parser.add_argument("--r2-path", help="R2 FASTQ file path")
    parser.add_argument("--genome-size", type=float, default=1.0, help="Genome size (default: 1.0)")
    parser.add_argument(
        "--total-reads-per-assembly",
        type=int,
        help="Project-specific total_reads_per_assembly override",
    )
    parser.add_argument(
        "--iterative-depth",
        type=int,
        help="Per-project iterative_assembly.depth override (e.g. 3)",
    )
    parser.add_argument(
        "--iterative-enabled",
        choices=["true", "false"],
        help="Per-project iterative_assembly.enabled override",
    )
    # TAREAN / seqclust per-project overrides (optional)
    parser.add_argument("--tarean-options", help="seqclust --options preset (e.g. ILLUMINA_SENSITIVE_BLASTPLUS)")
    parser.add_argument(
        "--tarean-mincl",
        type=float,
        help="Cluster threshold as fraction-of-reads (e.g. 0.0001 means 0.01%); converted internally for seqclust",
    )
    parser.add_argument("--tarean-assembly-min", type=int, help="seqclust --assembly_min")
    parser.add_argument(
        "--tarean-domain-search",
        choices=["BLASTX_W2", "BLASTX_W3", "DIAMOND"],
        help="seqclust --domain_search mode",
    )
    parser.add_argument("--tarean-prefix-length", type=int, help="seqclust --prefix_length (comparatives usually auto-set)")
    parser.add_argument(
        "--tarean-cleanup",
        dest="tarean_cleanup",
        action="store_true",
        help="Enable seqclust --cleanup",
    )
    parser.add_argument(
        "--no-tarean-cleanup",
        dest="tarean_cleanup",
        action="store_false",
        help="Disable seqclust --cleanup",
    )
    parser.set_defaults(tarean_cleanup=None)
    parser.add_argument(
        "--tarean-automatic-filtering",
        dest="tarean_automatic_filtering",
        action="store_true",
        help="Enable seqclust --automatic_filtering",
    )
    parser.add_argument(
        "--no-tarean-automatic-filtering",
        dest="tarean_automatic_filtering",
        action="store_false",
        help="Disable seqclust --automatic_filtering",
    )
    parser.set_defaults(tarean_automatic_filtering=None)
    parser.add_argument("--analysis-id", help="Comparative analysis ID")
    parser.add_argument("--samples", nargs="+", help="Sample names for comparative analysis")
    parser.add_argument("--analysis-description", help="Description for comparative analysis")
    parser.add_argument(
        "--keep-files",
        action="store_true",
        help="When used with drop-project, keep project files on disk",
    )
    
    args = parser.parse_args()
    
    pm = ProjectManager()
    
    if args.action == "create-project":
        if not all([args.project_id, args.taxonomy, args.description, args.ncbi_repeats]):
            print("❌ Missing required arguments for project creation")
            return

        tarean_overrides = {
            "options": args.tarean_options,
            "mincl": args.tarean_mincl,
            "assembly_min": args.tarean_assembly_min,
            "domain_search": args.tarean_domain_search,
            "prefix_length": args.tarean_prefix_length,
            "cleanup": args.tarean_cleanup,
            "automatic_filtering": args.tarean_automatic_filtering,
        }

        pm.create_project(
            args.project_id,
            args.taxonomy,
            args.description,
            args.ncbi_repeats,
            args.comparative_species or [],
            total_reads_per_assembly=args.total_reads_per_assembly,
            taxonomy_id=args.taxonomy_id,
            tarean_overrides=tarean_overrides,
        )
    
    elif args.action == "update-project":
        if not args.project_id:
            print("❌ Missing required argument --project-id")
            return
        if args.project_id not in pm.config.get("projects", {}):
            print(f"❌ Project '{args.project_id}' not found in config")
            return

        updated = False
        proj = pm.config["projects"][args.project_id]

        if args.total_reads_per_assembly is not None:
            proj["total_reads_per_assembly"] = int(args.total_reads_per_assembly)
            print(f"✅ Set total_reads_per_assembly={args.total_reads_per_assembly} for {args.project_id}")
            updated = True

        if args.iterative_depth is not None:
            proj.setdefault("iterative_assembly", {})
            proj["iterative_assembly"]["depth"] = int(args.iterative_depth)
            print(f"✅ Set iterative_assembly.depth={args.iterative_depth} for {args.project_id}")
            updated = True

        if args.iterative_enabled is not None:
            proj.setdefault("iterative_assembly", {})
            proj["iterative_assembly"]["enabled"] = (args.iterative_enabled == "true")
            print(f"✅ Set iterative_assembly.enabled={args.iterative_enabled} for {args.project_id}")
            updated = True

        if updated:
            pm._save_global_config(pm.config)
        else:
            print("ℹ️ No changes requested (provide flags to update).")

    elif args.action == "add-sample":
        if not all([args.project_id, args.sample_id, args.taxonomy, args.r1_path, args.r2_path]):
            print("❌ Missing required arguments for sample addition")
            return
        
        pm.add_sample(args.project_id, args.sample_id, args.taxonomy, args.r1_path, args.r2_path, args.genome_size)
    
    elif args.action == "list-projects":
        projects = pm.list_projects()
        if projects:
            print("📁 Projects:")
            for project_id in projects:
                config = pm.get_project_config(project_id)
                sample_count = len(config.get('samples', {}))
                print(f"  - {project_id}: {config['taxonomy']} ({sample_count} samples)")
        else:
            print("📁 No projects found")
    
    elif args.action == "list-samples":
        if not args.project_id:
            print("❌ Please specify --project-id to list samples")
            return
        samples = pm.get_project_samples(args.project_id)
        if samples:
            print(f"📄 Samples in project '{args.project_id}':")
            for sample_id in samples:
                prefix = pm.get_sample_prefix(args.project_id, sample_id)
                if prefix:
                    print(f"  - {sample_id} (prefix: {prefix})")
                else:
                    print(f"  - {sample_id}")
        else:
            print(f"📄 No samples found in project '{args.project_id}'")
    
    elif args.action == "show-project":
        if not args.project_id:
            print("❌ Please specify --project-id to show project details")
            return
        config = pm.get_project_config(args.project_id)
        if not config:
            print(f"❌ Project '{args.project_id}' not found")
            return
        samples = config.get("samples", {}) or {}
        comps = config.get("comparative_analyses", {}) or {}
        print(f"📁 Project '{args.project_id}':")
        print(f"   Taxonomy: {config.get('taxonomy', 'N/A')}")
        print(f"   Description: {config.get('description', 'N/A')}")
        print(f"   Samples: {len(samples)}")
        if samples:
            sample_list = ", ".join(sorted(samples.keys()))
            print(f"   Sample IDs: {sample_list}")
        print(f"   Comparative analyses: {len(comps)}")
        if comps:
            comp_list = ", ".join(sorted(comps.keys()))
            print(f"   Analyses: {comp_list}")
    
    elif args.action == "show-sample":
        if not (args.project_id and args.sample_id):
            print("❌ Please specify --project-id and --sample-id to show sample details")
            return
        metadata = pm.get_sample_metadata(args.project_id, args.sample_id)
        if not metadata:
            print(f"❌ Sample '{args.sample_id}' not found in project '{args.project_id}'")
            return
        prefix = pm.get_sample_prefix(args.project_id, args.sample_id)
        print(f"📄 Sample '{args.sample_id}' in project '{args.project_id}':")
        print(f"   Taxonomy: {metadata.get('taxonomy', 'N/A')}")
        print(f"   Genome size: {metadata.get('genome_size', 'N/A')}")
        if prefix:
            print(f"   Prefix: {prefix}")
        reads = metadata.get("read_files", {})
        if reads:
            print(f"   R1: {reads.get('R1', 'N/A')}")
            print(f"   R2: {reads.get('R2', 'N/A')}")
        print(f"   TAREAN status: {metadata.get('tarean_status', 'N/A')}")
    
    elif args.action == "validate":
        if args.project_id:
            if pm.validate_project(args.project_id):
                print(f"✅ Project '{args.project_id}' is valid")
                samples = pm.get_project_samples(args.project_id)
                for sample_id in samples:
                    if pm.validate_sample(args.project_id, sample_id):
                        print(f"  ✅ Sample '{sample_id}' is valid")
                    else:
                        print(f"  ❌ Sample '{sample_id}' has issues")

                analyses = pm.list_comparative_analyses(args.project_id)
                for analysis_id in analyses:
                    if pm.validate_comparative_analysis(args.project_id, analysis_id):
                        print(f"  ✅ Comparative '{analysis_id}' is valid")
                    else:
                        print(f"  ❌ Comparative '{analysis_id}' has issues")
            else:
                print(f"❌ Project '{args.project_id}' has issues")
        else:
            print("❌ Please specify --project-id for validation")
    
    elif args.action == "add-comparative":
        if not all([args.project_id, args.analysis_id, args.samples]):
            print("❌ Missing required arguments for comparative analysis addition")
            print("   Required: --project-id, --analysis-id, --samples")
            return
        
        pm.add_comparative_analysis(args.project_id, args.analysis_id, args.samples, args.analysis_description)
    
    elif args.action == "list-comparatives":
        if not args.project_id:
            print("❌ Please specify --project-id to list comparative analyses")
            return
        
        analyses = pm.list_comparative_analyses(args.project_id)
        if analyses:
            print(f"📊 Comparative analyses in project '{args.project_id}':")
            for analysis_id in analyses:
                analysis = pm.get_comparative_analysis(args.project_id, analysis_id)
                if analysis:
                    print(f"  - {analysis_id}: {analysis['description']}")
                    print(f"    Samples: {', '.join(analysis['samples'])}")
        else:
            print(f"📊 No comparative analyses found in project '{args.project_id}'")
    
    elif args.action == "remove-comparative":
        if not all([args.project_id, args.analysis_id]):
            print("❌ Missing required arguments for comparative analysis removal")
            print("   Required: --project-id, --analysis-id")
            return
        
        pm.remove_comparative_analysis(args.project_id, args.analysis_id)
    
    elif args.action == "update-comparative":
        if not all([args.project_id, args.analysis_id]):
            print("❌ Missing required arguments for comparative analysis update")
            print("   Required: --project-id, --analysis-id")
            print("   Optional: --samples, --analysis-description")
            return
        
        pm.update_comparative_analysis(args.project_id, args.analysis_id, args.samples, args.analysis_description)
    
    elif args.action == "update-prefixes":
        if not args.project_id:
            print("❌ Please specify --project-id to update sample prefixes")
            return
        
        pm.update_existing_samples_with_prefixes(args.project_id)
    
    elif args.action == "drop-project":
        if not args.project_id:
            print("❌ Please specify --project-id to drop a project")
            return
        
        pm.drop_project(args.project_id, delete_files=not args.keep_files)

if __name__ == "__main__":
    main() 