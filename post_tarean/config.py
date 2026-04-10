#!/usr/bin/env python3
"""
Configuration settings for RepOrtR Post-TAREAN Pipeline
Centralized configuration management

HOLY WORKFLOW COMPLIANCE: This module supports the modular workflow's post-processing.
All configuration integrates with projects/global_config.yaml.
"""

import os
from typing import Dict, List, Any
from dataclasses import dataclass, field
import logging

logger = logging.getLogger(__name__)

GLOBAL_CONFIG_PATH = os.environ.get("REPORTR_GLOBAL_CONFIG", "projects/global_config.yaml")


@dataclass
class BLASTConfig:
    """BLAST analysis configuration"""
    # Default BLAST tasks
    default_tasks: List[str] = field(default_factory=lambda: ['megablast', 'dc-megablast', 'blastn'])
    
    # Default databases
    default_dbs: List[str] = field(default_factory=lambda: ['local'])
    
    # Threading
    default_threads: int = 20
    
    # E-value thresholds
    evalue_threshold: float = 0.1
    evalue_max_olig: float = 0.001
    
    # Percent identity threshold
    pident_threshold: float = 60.0
    
    # Coverage thresholds
    coverage_thresholds: Dict[str, float] = field(default_factory=lambda: {
        'near-full': 0.95,
        'composite': 0.8,
        'partial': 0.6,
        'weak': 0.0
    })
    
    # BLAST executable path
    blastn_path: str = 'blastn'
    
    # Verbose output
    verbose: bool = True


@dataclass
class DatabaseConfig:
    """Database configuration"""
    # Database paths
    ncbi_db_path: str = './ncbi_repeats_db'
    local_db_path: str = './local_db_solo'
    comparative_db_path: str = './comparatives_db'
    reference_db_path: str = './important_db'
    
    # Database file patterns
    ncbi_patterns: List[str] = field(default_factory=lambda: ['*.nhr', '*.nin', '*.nsq'])
    fasta_patterns: List[str] = field(default_factory=lambda: ['*.fasta', '*.fa', '*.fas'])
    
    # Database names
    database_names: Dict[str, str] = field(default_factory=lambda: {
        'ncbi': 'ncbi_repeats.fasta',
        'ncbi_x3': 'ncbi_repeats_x3.fasta',
        'local': 'multifasta.fasta',
        'local_x3': 'multifasta_x3.fasta',
        'comp': 'COMPBASE.fasta',
        'comp_x3': 'COMPBASE_x3.fasta',
        'ref': 'reference.fasta'
    })


@dataclass
class AnalysisConfig:
    """Analysis configuration"""
    # Output settings
    output_dir: str = 'reports'
    create_excel: bool = True
    create_word: bool = True
    create_csv: bool = True
    
    # File naming
    output_prefix: str = 'repeat_analysis'
    
    # Sequence formatting
    sequence_width: int = 80
    
    # Cluster naming
    cluster_prefix: str = 'CL'
    cluster_padding: int = 4
    
    # Image settings
    default_image: str = None  # Let user specify or use actual images
    image_pattern: str = 'graph_layout.png'

    # Step orchestration
    # Ordered list of logical post-TAREAN steps to run for each subject.
    # Valid values currently: "blast", "summary".
    enabled_steps: List[str] = field(default_factory=lambda: ["blast", "summary", "quality_gating"])


@dataclass
class FilteringConfig:
    """Filtering configuration"""
    # Categories to exclude from normalization
    misc_categories: List[str] = field(default_factory=lambda: [
        'organelle',
        '|--plastid',
        '|--mitochondria',
        'Unclassified repeat (No evidence)',
        'contamination'
    ])
    
    # Levenshtein similarity threshold
    levenshtein_threshold: float = 90.0
    
    # Self-blast filtering
    remove_self_blasts: bool = True
    
    # E-value filtering for weak/partial alignments
    filter_weak_by_evalue: bool = True
    weak_evalue_threshold: float = 1.0


@dataclass
class ReportConfig:
    """Report generation configuration"""
    # Report types
    generate_summary: bool = True
    generate_detailed: bool = True
    generate_comparative: bool = True
    
    # Language settings
    language: str = 'en'  # 'en' for English, 'ru' for Russian
    
    # Formatting
    decimal_places: int = 2
    scientific_notation_precision: int = 2
    
    # Word document settings
    word_template: str = None
    word_style: str = 'default'
    
    # Excel settings
    excel_sheet_name: str = 'Repeat Analysis'
    excel_auto_filter: bool = True


@dataclass
class QualityGatingConfig:
    """Heuristic confidence scoring / quality gating configuration."""

    enabled: bool = True

    # Weights for composite confidence score (normalized internally).
    w_annotation: float = 0.35
    w_abundance: float = 0.25
    w_blast: float = 0.40

    # BLAST-derived thresholds (used to normalize scores).
    min_pident: float = 60.0
    max_evalue: float = 0.001

    # Recommendation cutoffs on composite confidence.
    publish: float = 0.85
    validate: float = 0.65
    refine: float = 0.40


@dataclass
class PipelineConfig:
    """Main pipeline configuration"""
    # Core components
    blast: BLASTConfig = field(default_factory=BLASTConfig)
    databases: DatabaseConfig = field(default_factory=DatabaseConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    filtering: FilteringConfig = field(default_factory=FilteringConfig)
    reports: ReportConfig = field(default_factory=ReportConfig)
    quality_gating: QualityGatingConfig = field(default_factory=QualityGatingConfig)
    
    # Logging
    log_level: str = 'INFO'
    log_format: str = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    log_file: str = None
    
    # Performance
    chunk_size: int = 1000
    max_memory_mb: int = 2048
    
    # Error handling
    continue_on_error: bool = True
    max_retries: int = 3
    
    # Validation
    validate_inputs: bool = True
    validate_outputs: bool = True


class ConfigManager:
    """Configuration manager for the pipeline"""
    
    def __init__(self, config_file: str = None, project_id: str = None):
        """
        If project_id is provided, load post-TAREAN parameters for that project
        from projects/global_config.yaml first, then optional config_file
        overrides.
        """
        self.config = PipelineConfig()
        self.project_id = project_id
        
        # Try to load from global config first
        if project_id:
            self._load_from_global_config(project_id)
            # Canonical in-workflow LocalDB location (project-scoped).
            # This is required for legacy post_tarean BLAST ("local" DB).
            if getattr(self.config.databases, "local_db_path", "") in ("", "./local_db_solo"):
                self.config.databases.local_db_path = f"projects/{project_id}/blast_db"
        
        # Then load from specific config file if provided
        if config_file and os.path.exists(config_file):
            self.load_config(config_file)
        
        self._setup_logging()
    
    def _load_from_global_config(self, project_id: str):
        """Load configuration from global config file"""
        try:
            import yaml
            global_config_path = GLOBAL_CONFIG_PATH
            if os.path.exists(global_config_path):
                with open(global_config_path, 'r') as f:
                    global_config = yaml.safe_load(f)
                
                if project_id in global_config.get("projects", {}):
                    project_config = global_config["projects"][project_id]
                    if "post_tarean_params" in project_config:
                        self._update_from_global_config(project_config["post_tarean_params"])
                        logger.info(f"Loaded post-TAREAN config for project {project_id} from global config")
                    else:
                        # Use defaults if no project-specific config
                        defaults = global_config.get("defaults", {}).get("post_tarean_params", {})
                        if defaults:
                            self._update_from_global_config(defaults)
                            logger.info(f"Loaded default post-TAREAN config for project {project_id}")
        except Exception as e:
            logger.warning(f"Failed to load from global config: {e}")
    
    def _update_from_global_config(self, config_dict: Dict[str, Any]):
        """Update configuration from global config dictionary"""
        # Update BLAST config
        if "blast" in config_dict:
            for key, value in config_dict["blast"].items():
                if hasattr(self.config.blast, key):
                    setattr(self.config.blast, key, value)
        
        # Update database config
        if "databases" in config_dict:
            for key, value in config_dict["databases"].items():
                if hasattr(self.config.databases, key):
                    setattr(self.config.databases, key, value)
        
        # Update analysis config
        if "analysis" in config_dict:
            for key, value in config_dict["analysis"].items():
                if hasattr(self.config.analysis, key):
                    setattr(self.config.analysis, key, value)
        
        # Update filtering config
        if "filtering" in config_dict:
            for key, value in config_dict["filtering"].items():
                if hasattr(self.config.filtering, key):
                    setattr(self.config.filtering, key, value)
        
        # Update reports config
        if "reports" in config_dict:
            for key, value in config_dict["reports"].items():
                if hasattr(self.config.reports, key):
                    setattr(self.config.reports, key, value)

        # Update quality gating config
        if "quality_gating" in config_dict:
            for key, value in config_dict["quality_gating"].items():
                if hasattr(self.config.quality_gating, key):
                    setattr(self.config.quality_gating, key, value)
        
        # Update performance config
        if "performance" in config_dict:
            for key, value in config_dict["performance"].items():
                if hasattr(self.config, key):
                    setattr(self.config, key, value)
    
    def load_config(self, config_file: str):
        """Load configuration from file"""
        try:
            import yaml
            with open(config_file, 'r') as f:
                config_data = yaml.safe_load(f)
            self._update_config_from_dict(config_data)
            logger.info(f"Configuration loaded from {config_file}")
        except ImportError:
            logger.warning("PyYAML not available, using default configuration")
        except Exception as e:
            logger.error(f"Failed to load configuration: {e}")
    
    def _update_config_from_dict(self, config_dict: Dict[str, Any]):
        """Update configuration from dictionary"""
        for section, values in config_dict.items():
            if hasattr(self.config, section):
                section_config = getattr(self.config, section)
                for key, value in values.items():
                    if hasattr(section_config, key):
                        setattr(section_config, key, value)
    
    def _setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=getattr(logging, self.config.log_level),
            format=self.config.log_format,
            filename=self.config.log_file
        )
    
    def get_blast_config(self) -> BLASTConfig:
        """Get BLAST configuration"""
        return self.config.blast
    
    def get_database_config(self) -> DatabaseConfig:
        """Get database configuration"""
        return self.config.databases
    
    def get_analysis_config(self) -> AnalysisConfig:
        """Get analysis configuration"""
        return self.config.analysis
    
    def get_filtering_config(self) -> FilteringConfig:
        """Get filtering configuration"""
        return self.config.filtering
    
    def get_report_config(self) -> ReportConfig:
        """Get report configuration"""
        return self.config.reports
    
    def update_config(self, **kwargs):
        """Update configuration parameters"""
        for key, value in kwargs.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)
            else:
                logger.warning(f"Unknown configuration parameter: {key}")
    
    def save_config(self, config_file: str):
        """Save current configuration to file"""
        try:
            import yaml
            config_dict = self._config_to_dict()
            with open(config_file, 'w') as f:
                yaml.dump(config_dict, f, default_flow_style=False)
            logger.info(f"Configuration saved to {config_file}")
        except ImportError:
            logger.warning("PyYAML not available, cannot save configuration")
        except Exception as e:
            logger.error(f"Failed to save configuration: {e}")
    
    def _config_to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary"""
        config_dict = {}
        for section_name in ['blast', 'databases', 'analysis', 'filtering', 'reports', 'quality_gating']:
            section = getattr(self.config, section_name)
            config_dict[section_name] = {
                key: getattr(section, key) 
                for key in section.__dataclass_fields__.keys()
            }
        
        # Add top-level config
        config_dict.update({
            'log_level': self.config.log_level,
            'log_format': self.config.log_format,
            'log_file': self.config.log_file,
            'chunk_size': self.config.chunk_size,
            'max_memory_mb': self.config.max_memory_mb,
            'continue_on_error': self.config.continue_on_error,
            'max_retries': self.config.max_retries,
            'validate_inputs': self.config.validate_inputs,
            'validate_outputs': self.config.validate_outputs,
        })
        
        return config_dict
    
    def validate_config(self) -> bool:
        """Validate configuration settings"""
        errors = []
        
        # Validate BLAST settings
        if self.config.blast.evalue_threshold <= 0:
            errors.append("E-value threshold must be positive")
        
        if self.config.blast.pident_threshold < 0 or self.config.blast.pident_threshold > 100:
            errors.append("Percent identity threshold must be between 0 and 100")
        
        # Validate database paths
        for path_name, path in [
            ('NCBI DB', self.config.databases.ncbi_db_path),
            ('Local DB', self.config.databases.local_db_path),
            ('Comparative DB', self.config.databases.comparative_db_path),
            ('Reference DB', self.config.databases.reference_db_path)
        ]:
            if not os.path.exists(path):
                logger.warning(f"{path_name} path does not exist: {path}")
        
        # Validate output directory
        if not os.path.exists(self.config.analysis.output_dir):
            try:
                os.makedirs(self.config.analysis.output_dir, exist_ok=True)
            except Exception as e:
                errors.append(f"Cannot create output directory: {e}")
        
        if errors:
            for error in errors:
                logger.error(f"Configuration error: {error}")
            return False
        
        logger.info("Configuration validation passed")
        return True

    # ------------------------------------------------------------------
    # Backwards-compatible helper for pipeline.py
    # ------------------------------------------------------------------

    def load_post_tarean_config(self, project_id: str = None) -> PipelineConfig:
        """
        Legacy-style helper used by post_tarean.pipeline.

        Returns a PipelineConfig instance with project-specific post_tarean
        settings applied (if available).

        If project_id is provided, a new ConfigManager scoped to that project
        is created and its config returned. Otherwise, this instance's config
        is returned.
        """
        if project_id:
            mgr = ConfigManager(project_id=project_id)
            return mgr.config
        return self.config


# Default configuration instance
default_config = ConfigManager()


def get_config() -> ConfigManager:
    """Get default configuration instance"""
    return default_config


def create_config_from_file(config_file: str) -> ConfigManager:
    """Create configuration from file"""
    return ConfigManager(config_file)


def create_config_for_project(project_id: str, config_file: str = None) -> ConfigManager:
    """Create configuration for a specific project"""
    return ConfigManager(config_file, project_id)


def create_custom_config(**kwargs) -> ConfigManager:
    """Create custom configuration with overrides"""
    config = ConfigManager()
    config.update_config(**kwargs)
    return config 