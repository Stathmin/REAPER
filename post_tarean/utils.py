#!/usr/bin/env python3
"""
Utility functions for RepOrtR Post-TAREAN Pipeline
Common functions used across multiple modules
"""

import os
import re
import glob
import pandas as pd
import textwrap
from typing import Dict, List, Optional, Tuple
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def get_cyphers() -> Dict[str, List[str]]:
    """Return per-sample metadata mapping from projects/global_config.yaml.

    Legacy code used to read `cypher.csv` from the current working directory.
    That implicit dependency is forbidden in the modern codebase; the workflow's
    source of truth is `projects/global_config.yaml`.

    Returns a dict: sample_id -> [organism, accession_or_empty, genome_size_or_empty]
    """
    import yaml

    global_cfg_path = os.environ.get("REPORTR_GLOBAL_CONFIG", "projects/global_config.yaml")
    with open(global_cfg_path, "r") as f:
        cfg = yaml.safe_load(f) or {}

    out: Dict[str, List[str]] = {}
    projects = (cfg.get("projects") or {})
    if not isinstance(projects, dict):
        raise TypeError("global config: projects must be a mapping")

    for proj_id, proj_cfg in projects.items():
        if not isinstance(proj_cfg, dict):
            continue
        samples = (proj_cfg.get("samples") or {})
        if not isinstance(samples, dict):
            continue
        for sample_id, sm in samples.items():
            if not isinstance(sm, dict):
                continue
            organism = sm.get("organism")
            genome_size = sm.get("genome_size")
            if organism is None or str(organism).strip() == "":
                raise KeyError(f"Missing required projects.{proj_id}.samples.{sample_id}.organism in {global_cfg_path}")
            if genome_size is None:
                raise KeyError(f"Missing required projects.{proj_id}.samples.{sample_id}.genome_size in {global_cfg_path}")
            out[str(sample_id)] = [str(organism), "", str(genome_size)]

    return out


def get_ncbi_naming() -> pd.DataFrame:
    """Load NCBI naming data"""
    try:
        with open('./ncbi_naming.csv', 'r') as file:
            return pd.read_csv(file).drop('Unnamed: 0', axis=1, errors='ignore')
    except FileNotFoundError:
        logger.warning("ncbi_naming.csv not found")
        return pd.DataFrame()


def get_tareans(path, index):
    """Parse TAREAN output data with proper image mapping"""
    repeats_paths = glob.glob(f'{path}/*TAREAN*')
    names, seqs, clusters = [], [], []
    
    # Create image mapping like the original script
    img_map = {}
    for img_path in glob.glob(f'{path}/seqclust/clustering/clusters/dir_*/graph_layout.png'):
        dir_name = os.path.basename(os.path.dirname(img_path))
        cluster_id = dir_name.replace('dir_', '')
        img_map[cluster_id] = img_path
    
    for repeat_path in repeats_paths:
        with open(repeat_path) as f:
            content = f.read().split('>')
            for fasta in filter(None, content):
                lines = list(filter(None, fasta.splitlines()))
                header = lines[0]
                seq = ''.join(lines[1:])
                names.append(f"{index}_{header}")
                clusters.append(re.search(r'CL\d+', header).group())
                seqs.append(textwrap.fill(seq, width=80))

    # Map cluster IDs to image paths
    cluster_ids = [f"CL{int(cl[2:]):04d}" for cl in clusters]
    pics = [img_map.get(cid, None) for cid in cluster_ids]  # Handle missing images properly
    
    return pd.DataFrame({
        'cluster': names,
        'seq': seqs,
        'pic_path': pics,
        'Cluster': clusters
    })


def parse_copy_data(path):
    """Parse copy number data from CLUSTER_TABLE.csv"""
    cluster_table = Path(str(path)) / "CLUSTER_TABLE.csv"
    if not cluster_table.exists():
        raise FileNotFoundError(f"Missing required CLUSTER_TABLE.csv at {cluster_table}")
    # Strict: this table is required for downstream copy/abundance reporting.
    # RE2 emits a short 2-column metadata preamble, then a real header row.
    # We must skip the preamble to avoid pandas ParserError on ragged rows.
    skiprows = 0
    analyzed_reads = None
    with cluster_table.open("r") as f:
        for i, line in enumerate(f):
            # Parse metadata preamble like: "Number_of_analyzed_reads"\t100000
            if analyzed_reads is None and line.startswith('"Number_of_analyzed_reads"'):
                try:
                    analyzed_reads = int(line.strip().split("\t", 1)[1].strip().strip('"'))
                except Exception:
                    analyzed_reads = None
            if line.lstrip().startswith('"Cluster"\t') or line.lstrip().startswith("Cluster\t"):
                skiprows = i
                break
    try:
        df = pd.read_csv(cluster_table, sep="\t", engine="python", skiprows=skiprows)
        # Align merge key with TAREAN parser (e.g. "CL2", "CL10", ...).
        if "Cluster" in df.columns:
            df["Cluster"] = df["Cluster"].apply(lambda x: f"CL{int(x)}" if str(x).strip() != "" else x)
        # Legacy reports expect `size, %` like satMiner tables.
        if "size, %" not in df.columns:
            base_col = "Size_adjusted" if "Size_adjusted" in df.columns else ("Size" if "Size" in df.columns else None)
            if base_col is not None and analyzed_reads:
                df["size, %"] = (100.0 * pd.to_numeric(df[base_col], errors="coerce") / float(analyzed_reads)).round(2)
        return df
    except Exception as e:
        raise RuntimeError(f"Failed to parse {cluster_table}: {e}") from e


def parse_comparative_data(index, path):
    """Parse comparative analysis data"""
    try:
        cluster_table_path = os.path.join(str(path), "COMPARATIVE_CLUSTER_TABLE.csv")
        counts_table_path = os.path.join(str(path), "COMPARATIVE_ANALYSIS_COUNTS.csv")

        if os.path.exists(cluster_table_path):
            # seqclust comparative mode output
            df = pd.read_csv(cluster_table_path, skiprows=2, sep="\t", nrows=500)
        else:
            df = pd.read_csv(counts_table_path, skiprows=2, sep="\t", nrows=500)
        
        # Handle different column name formats
        cluster_col = 'cluster' if 'cluster' in df.columns else 'Cluster'
        supercluster_col = 'supercluster' if 'supercluster' in df.columns else 'Supercluster'
        
        # Species/sample count columns are everything except cluster identifiers.
        index_cols = [c for c in df.columns if c not in {cluster_col, supercluster_col}]
        if not index_cols:
            return pd.DataFrame()
        total_per_row = df[index_cols].sum(axis=1)
        total_per_row = total_per_row.replace(0, pd.NA)
        
        for col in index_cols:
            df[f'{col}, %'] = (100 * df[col] / total_per_row).map('{:,.2f}'.format)
        
        df['merger'] = index + '_CL' + df[cluster_col].astype(str)
        # Align merge key with TAREAN parser (e.g. "CL2", "CL10", ...).
        df["Cluster"] = df[cluster_col].apply(lambda x: f"CL{int(x)}" if str(x).strip() != "" else x)

        # Drop original identifier columns except canonical `Cluster`.
        if cluster_col != "Cluster" and cluster_col in df.columns:
            df = df.drop(columns=[cluster_col])
        if supercluster_col in df.columns:
            df = df.drop(columns=[supercluster_col])

        return df
    except FileNotFoundError:
        logger.warning(
            f"Comparative table not found in {path} "
            f"(expected COMPARATIVE_CLUSTER_TABLE.csv or COMPARATIVE_ANALYSIS_COUNTS.csv)"
        )
        return pd.DataFrame()


def get_coverage_description(coverage: float) -> str:
    """Get human-readable coverage description"""
    if 0.6 < coverage < 0.8:
        return 'по части длины'
    elif 0.8 <= coverage < 0.9:
        return 'по большей части длины'
    elif 0.9 <= coverage < 0.95:
        return 'практически по всей длине'
    elif coverage >= 0.95:
        return 'по всей длине'
    else:
        return 'слабое покрытие'


def get_task_description(task: str) -> str:
    """Get human-readable task description"""
    task_descriptions = {
        'blastn': 'некую гомологию',
        'dc-megablast': 'умеренную гомологию',
        'megablast': 'высокую гомологию'
    }
    return task_descriptions.get(task, 'гомологию')


def get_copy_number_description(copy_number: float) -> str:
    """Get human-readable copy number description"""
    if copy_number < 0.1:
        return 'низкокопийный'
    elif 0.1 <= copy_number < 0.5:
        return 'среднекопийный'
    elif copy_number >= 0.5:
        return 'высококопийный'
    else:
        return 'неизвестная копийность'


def get_better_name(name: str) -> str:
    """Extract better name from sequence identifier"""
    return name.split('_TR_')[0] if '_TR_' in name else name


def wrap_sequence(seq: str, width: int = 80) -> str:
    """Wrap sequence to specified width"""
    return textwrap.fill(seq, width=width)


def extract_cluster_number(cluster_name: str) -> int:
    """Extract cluster number from cluster name"""
    match = re.search(r'CL(\d+)', cluster_name)
    return int(match.group(1)) if match else 0


def format_percentage(value: float, decimal_places: int = 2) -> str:
    """Format percentage value"""
    return f"{value:.{decimal_places}f}%"


def safe_merge_dataframes(df1: pd.DataFrame, df2: pd.DataFrame, 
                         on: str, how: str = 'left') -> pd.DataFrame:
    """Safely merge two dataframes with error handling"""
    try:
        return pd.merge(df1, df2, on=on, how=how)
    except Exception as e:
        logger.warning(f"Merge failed: {e}")
        return df1


def validate_file_path(path: str) -> bool:
    """Validate if file path exists"""
    return os.path.exists(path)


def create_output_directory(directory: str) -> bool:
    """Create output directory if it doesn't exist"""
    try:
        os.makedirs(directory, exist_ok=True)
        return True
    except Exception as e:
        logger.error(f"Failed to create directory {directory}: {e}")
        return False


def get_file_size_mb(file_path: str) -> float:
    """Get file size in megabytes"""
    try:
        return os.path.getsize(file_path) / (1024 * 1024)
    except OSError:
        return 0.0


def clean_filename(filename: str) -> str:
    """Clean filename for safe file operations"""
    # Remove or replace problematic characters
    cleaned = re.sub(r'[<>:"/\\|?*]', '_', filename)
    # Remove leading/trailing spaces and dots
    cleaned = cleaned.strip('. ')
    return cleaned


def get_unique_filename(base_path: str, extension: str = '') -> str:
    """Generate unique filename to avoid overwrites"""
    counter = 1
    original_path = base_path
    while os.path.exists(base_path + extension):
        base_path = f"{original_path}_{counter}"
        counter += 1
    return base_path + extension


def parse_summarized_annotation(path: str) -> pd.DataFrame:
    """Parse summarized annotation file"""
    try:
        df = pd.read_fwf(path, skiprows=13, widths=[49, 16, 17, 12, 30])
        df.columns = ['type', 'percent', 'scl', 'cl', 'reads']
        
        # Drop problematic rows
        df = df.drop([74, 75, 79, 80, 82, 83], errors='ignore')
        
        # Clean numeric columns
        for column in ['percent', 'scl', 'cl', 'reads']:
            df[column] = df[column].apply(
                lambda x: float(str(x).replace(' ', '').replace('|', '').replace('</pre>', ''))
            )
        
        # Clean type column
        df['type'] = df['type'].apply(
            lambda x: re.sub(r'[ ]+$', '', ''.join(x[0:-1])).replace("'", "|")
        )
        
        return df[['type', 'percent']]
    except Exception as e:
        logger.error(f"Failed to parse summarized annotation: {e}")
        return pd.DataFrame()


def normalize_repeat_percentages(df: pd.DataFrame, sample_name: str, 
                               misc_categories: List[str]) -> pd.DataFrame:
    """Normalize repeat percentages by excluding misc categories"""
    df = df.copy()
    df.columns = ['type', sample_name]
    df = df.set_index('type')
    
    # Calculate normalization factor
    not_clusters = df.loc[misc_categories][sample_name].sum()
    factor = (100 - not_clusters) / 100
    
    # Apply normalization
    df[f'{sample_name} normalized'] = df[sample_name].apply(lambda x: round(x * factor, 2))
    df.loc[misc_categories, f'{sample_name} normalized'] = ''
    
    return df


def merge_multiple_samples(dfs: List[pd.DataFrame]) -> pd.DataFrame:
    """Merge multiple sample dataframes"""
    if not dfs:
        return pd.DataFrame()
    
    result = dfs[0]
    for other_df in dfs[1:]:
        result = pd.merge(left=result, right=other_df, left_on='type', right_on='type', how='outer')
    
    return result


def calculate_depth_statistics(df: pd.DataFrame, depth_dict: Dict[str, int]) -> Dict[str, float]:
    """Calculate statistics by repeat depth"""
    stats = {}
    
    for repeat_type, depth in depth_dict.items():
        if repeat_type in df.index:
            stats[f'depth_{depth}'] = df.loc[repeat_type].sum()
    
    return stats


def generate_summary_statistics(data: pd.DataFrame) -> Dict[str, any]:
    """Generate summary statistics for repeat data"""
    stats = {
        'total_clusters': len(data),
        'total_coverage': data.get('size, %', pd.Series()).sum(),
        'mean_cluster_size': data.get('size, %', pd.Series()).mean(),
        'max_cluster_size': data.get('size, %', pd.Series()).max(),
        'min_cluster_size': data.get('size, %', pd.Series()).min(),
    }
    
    # Add BLAST statistics if available
    if 'best_evalue' in data.columns:
        stats.update({
            'clusters_with_blast_hits': data['best_evalue'].notna().sum(),
            'mean_best_evalue': data['best_evalue'].mean(),
            'min_best_evalue': data['best_evalue'].min(),
        })
    
    return stats


def format_scientific_notation(value: float, precision: int = 2) -> str:
    """Format number in scientific notation"""
    return f"{value:.{precision}e}"


def extract_genera_from_blast_hits(blast_data: pd.DataFrame) -> List[str]:
    """Extract genera from BLAST hit descriptions"""
    if 'full' not in blast_data.columns:
        return []
    
    genera = set()
    for description in blast_data['full'].dropna():
        parts = str(description).split()
        if len(parts) > 1:
            potential_genus = parts[1]
            if re.match(r"^[A-Z][a-z]+$", potential_genus):
                genera.add(potential_genus)
    
    return sorted(list(genera))


def create_progress_tracker(total_items: int, description: str = "Processing"):
    """Create a simple progress tracker"""
    class ProgressTracker:
        def __init__(self, total, desc):
            self.total = total
            self.current = 0
            self.desc = desc
        
        def update(self, increment=1):
            self.current += increment
            percentage = (self.current / self.total) * 100
            logger.info(f"{self.desc}: {self.current}/{self.total} ({percentage:.1f}%)")
        
        def finish(self):
            logger.info(f"{self.desc}: Completed {self.current}/{self.total}")
    
    return ProgressTracker(total_items, description) 