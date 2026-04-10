#!/usr/bin/env python3
"""
Flow Cytometry Correlation Module for RepOrtR

This module correlates genome size data from flow cytometry with repeat content
from TAREAN analysis to validate computational predictions.

Author: RepOrtR Team
Date: 2025
"""

import os
import sys
import argparse
import json
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# BioPython imports
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

@dataclass
class CytometryData:
    """Represents flow cytometry data for a sample"""
    sample_id: str
    genome_size_gb: float
    peak_fluorescence: float
    coefficient_variation: float
    cell_count: int
    standard_sample: str
    measurement_date: str
    
    def __post_init__(self):
        if not hasattr(self, 'measurement_date'):
            from datetime import datetime
            self.measurement_date = datetime.now().isoformat()

@dataclass
class RepeatContentData:
    """Represents repeat content data from TAREAN analysis"""
    sample_id: str
    total_repeat_percentage: float
    satellite_percentage: float
    ltr_percentage: float
    line_percentage: float
    sine_percentage: float
    dna_percentage: float
    unknown_percentage: float
    total_consensus_count: int
    total_repeat_length: int
    
    def __post_init__(self):
        # Calculate total if not provided
        if not hasattr(self, 'total_repeat_percentage'):
            self.total_repeat_percentage = (
                self.satellite_percentage + 
                self.ltr_percentage + 
                self.line_percentage + 
                self.sine_percentage + 
                self.dna_percentage + 
                self.unknown_percentage
            )

class FlowCytometryCorrelator:
    """
    Flow Cytometry Correlation Module
    
    This class correlates genome size data from flow cytometry with repeat content
    from TAREAN analysis to validate computational predictions.
    """
    
    def __init__(self, project_id: str = None):
        """
        Initialize the flow cytometry correlator
        
        Args:
            project_id: Project identifier for data organization
        """
        self.project_id = project_id or "default"
        self.cytometry_data = []
        self.repeat_data = []
        self.correlation_results = {}
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Statistical parameters
        self.correlation_methods = ['pearson', 'spearman']
        self.significance_level = 0.05
        
        # Visualization settings
        plt.style.use('seaborn-v0_8')
        self.figure_size = (12, 8)
        self.dpi = 300
    
    def load_cytometry_data(self, data_file: str) -> List[CytometryData]:
        """
        Load flow cytometry data from file
        
        Args:
            data_file: Path to cytometry data file (CSV format)
            
        Returns:
            List of CytometryData objects
        """
        self.logger.info(f"Loading cytometry data from {data_file}")
        
        try:
            df = pd.read_csv(data_file)
            cytometry_data = []
            
            for _, row in df.iterrows():
                cytometry = CytometryData(
                    sample_id=row['sample_id'],
                    genome_size_gb=float(row['genome_size_gb']),
                    peak_fluorescence=float(row['peak_fluorescence']),
                    coefficient_variation=float(row['coefficient_variation']),
                    cell_count=int(row['cell_count']),
                    standard_sample=row['standard_sample'],
                    measurement_date=row.get('measurement_date', '')
                )
                cytometry_data.append(cytometry)
            
            self.cytometry_data = cytometry_data
            self.logger.info(f"Loaded {len(cytometry_data)} cytometry measurements")
            return cytometry_data
            
        except Exception as e:
            self.logger.error(f"Failed to load cytometry data: {e}")
            return []
    
    def load_repeat_data_from_tarean(self, tarean_dirs: List[str]) -> List[RepeatContentData]:
        """
        Load repeat content data from TAREAN analysis directories
        
        Args:
            tarean_dirs: List of TAREAN analysis directory paths
            
        Returns:
            List of RepeatContentData objects
        """
        self.logger.info(f"Loading repeat data from {len(tarean_dirs)} TAREAN directories")
        
        repeat_data = []
        
        for tarean_dir in tarean_dirs:
            sample_path = Path(tarean_dir)
            sample_id = sample_path.name
            
            try:
                # Parse TAREAN results
                repeat_content = self._parse_tarean_results(sample_path)
                
                if repeat_content:
                    repeat_data.append(repeat_content)
                    self.logger.debug(f"Loaded repeat data for {sample_id}")
                
            except Exception as e:
                self.logger.error(f"Failed to load repeat data for {sample_id}: {e}")
        
        self.repeat_data = repeat_data
        self.logger.info(f"Loaded repeat data for {len(repeat_data)} samples")
        return repeat_data
    
    def _parse_tarean_results(self, tarean_dir: Path) -> Optional[RepeatContentData]:
        """
        Parse TAREAN analysis results to extract repeat content
        
        Args:
            tarean_dir: Path to TAREAN analysis directory
            
        Returns:
            RepeatContentData object or None if parsing fails
        """
        sample_id = tarean_dir.name
        
        # Look for TAREAN summary files
        summary_files = [
            "TAREAN_summary.txt",
            "analysis_summary.txt",
            "repeat_summary.csv"
        ]
        
        repeat_content = {
            'satellite_percentage': 0.0,
            'ltr_percentage': 0.0,
            'line_percentage': 0.0,
            'sine_percentage': 0.0,
            'dna_percentage': 0.0,
            'unknown_percentage': 0.0,
            'total_consensus_count': 0,
            'total_repeat_length': 0
        }
        
        # Try to find and parse summary file
        for summary_file in summary_files:
            file_path = tarean_dir / summary_file
            if file_path.exists():
                try:
                    repeat_content = self._parse_summary_file(file_path)
                    break
                except Exception as e:
                    self.logger.debug(f"Failed to parse {summary_file}: {e}")
        
        # If no summary file, estimate from consensus sequences
        if repeat_content['total_consensus_count'] == 0:
            repeat_content = self._estimate_from_consensus(tarean_dir)
        
        if repeat_content['total_consensus_count'] > 0:
            # Normalize percentages to sum to 100
            total_percentage = repeat_content['satellite_percentage'] + repeat_content['ltr_percentage'] + repeat_content['line_percentage'] + repeat_content['sine_percentage'] + repeat_content['dna_percentage'] + repeat_content['unknown_percentage']
            
            if total_percentage > 0:
                normalized_satellite = (repeat_content['satellite_percentage'] / total_percentage) * 100
                normalized_ltr = (repeat_content['ltr_percentage'] / total_percentage) * 100
                normalized_line = (repeat_content['line_percentage'] / total_percentage) * 100
                normalized_sine = (repeat_content['sine_percentage'] / total_percentage) * 100
                normalized_dna = (repeat_content['dna_percentage'] / total_percentage) * 100
                normalized_unknown = (repeat_content['unknown_percentage'] / total_percentage) * 100
            else:
                normalized_satellite = normalized_ltr = normalized_line = normalized_sine = normalized_dna = normalized_unknown = 0
            
            return RepeatContentData(
                sample_id=sample_id,
                total_repeat_percentage=normalized_satellite + normalized_ltr + normalized_line + normalized_sine + normalized_dna + normalized_unknown,
                satellite_percentage=normalized_satellite,
                ltr_percentage=normalized_ltr,
                line_percentage=normalized_line,
                sine_percentage=normalized_sine,
                dna_percentage=normalized_dna,
                unknown_percentage=normalized_unknown,
                total_consensus_count=repeat_content['total_consensus_count'],
                total_repeat_length=repeat_content['total_repeat_length']
            )
        
        return None
    
    def _parse_summary_file(self, file_path: Path) -> Dict:
        """Parse TAREAN summary file to extract repeat content"""
        repeat_content = {
            'satellite_percentage': 0.0,
            'ltr_percentage': 0.0,
            'line_percentage': 0.0,
            'sine_percentage': 0.0,
            'dna_percentage': 0.0,
            'unknown_percentage': 0.0,
            'total_consensus_count': 0,
            'total_repeat_length': 0
        }
        
        with open(file_path, 'r') as f:
            content = f.read()
            
            # Extract percentages (basic parsing)
            lines = content.split('\n')
            for line in lines:
                if 'satellite' in line.lower() and '%' in line:
                    try:
                        repeat_content['satellite_percentage'] = float(line.split('%')[0].split()[-1])
                    except:
                        pass
                elif 'ltr' in line.lower() and '%' in line:
                    try:
                        repeat_content['ltr_percentage'] = float(line.split('%')[0].split()[-1])
                    except:
                        pass
                elif 'line' in line.lower() and '%' in line:
                    try:
                        repeat_content['line_percentage'] = float(line.split('%')[0].split()[-1])
                    except:
                        pass
                elif 'sine' in line.lower() and '%' in line:
                    try:
                        repeat_content['sine_percentage'] = float(line.split('%')[0].split()[-1])
                    except:
                        pass
                elif 'dna' in line.lower() and '%' in line:
                    try:
                        repeat_content['dna_percentage'] = float(line.split('%')[0].split()[-1])
                    except:
                        pass
                elif 'total' in line.lower() and 'consensus' in line.lower():
                    try:
                        repeat_content['total_consensus_count'] = int(line.split()[-1])
                    except:
                        pass
        
        return repeat_content
    
    def _estimate_from_consensus(self, tarean_dir: Path) -> Dict:
        """Estimate repeat content from consensus sequences"""
        repeat_content = {
            'satellite_percentage': 0.0,
            'ltr_percentage': 0.0,
            'line_percentage': 0.0,
            'sine_percentage': 0.0,
            'dna_percentage': 0.0,
            'unknown_percentage': 0.0,
            'total_consensus_count': 0,
            'total_repeat_length': 0
        }
        
        consensus_files = [
            "TAREAN_consensus_rank_1.fasta",  # Satellites
            "TAREAN_consensus_rank_2.fasta",  # Low confidence satellites
            "TAREAN_consensus_rank_3.fasta",  # LTR elements
            "TAREAN_consensus_rank_4.fasta",  # rDNA
        ]
        
        total_length = 0
        consensus_count = 0
        
        for i, filename in enumerate(consensus_files):
            file_path = tarean_dir / filename
            if file_path.exists():
                file_length = 0
                file_count = 0
                
                for record in SeqIO.parse(file_path, "fasta"):
                    file_length += len(record.seq)
                    file_count += 1
                
                total_length += file_length
                consensus_count += file_count
                
                # Estimate percentages based on rank
                if i == 0:  # Rank 1 - satellites
                    repeat_content['satellite_percentage'] = (file_length / total_length) * 100 if total_length > 0 else 0
                elif i == 2:  # Rank 3 - LTR
                    repeat_content['ltr_percentage'] = (file_length / total_length) * 100 if total_length > 0 else 0
                else:  # Other ranks
                    repeat_content['unknown_percentage'] += (file_length / total_length) * 100 if total_length > 0 else 0
        
        repeat_content['total_consensus_count'] = consensus_count
        repeat_content['total_repeat_length'] = total_length
        
        return repeat_content
    
    def correlate_genome_size_repeats(self) -> Dict:
        """
        Correlate genome size with repeat content
        
        Returns:
            Dictionary with correlation results and statistics
        """
        self.logger.info("Correlating genome size with repeat content")
        
        if not self.cytometry_data or not self.repeat_data:
            self.logger.error("No cytometry or repeat data available")
            return {}
        
        # Create combined dataset
        combined_data = []
        
        for cytometry in self.cytometry_data:
            # Find matching repeat data
            repeat_match = None
            for repeat in self.repeat_data:
                if repeat.sample_id == cytometry.sample_id:
                    repeat_match = repeat
                    break
            
            if repeat_match:
                combined_data.append({
                    'sample_id': cytometry.sample_id,
                    'genome_size_gb': cytometry.genome_size_gb,
                    'total_repeat_percentage': repeat_match.total_repeat_percentage,
                    'satellite_percentage': repeat_match.satellite_percentage,
                    'ltr_percentage': repeat_match.ltr_percentage,
                    'line_percentage': repeat_match.line_percentage,
                    'sine_percentage': repeat_match.sine_percentage,
                    'dna_percentage': repeat_match.dna_percentage,
                    'unknown_percentage': repeat_match.unknown_percentage,
                    'total_consensus_count': repeat_match.total_consensus_count,
                    'total_repeat_length': repeat_match.total_repeat_length
                })
        
        if not combined_data:
            self.logger.error("No matching samples found between cytometry and repeat data")
            self.logger.debug(f"Cytometry samples: {[c.sample_id for c in self.cytometry_data]}")
            self.logger.debug(f"Repeat samples: {[r.sample_id for r in self.repeat_data]}")
            return {}
        
        # Create DataFrame for analysis
        df = pd.DataFrame(combined_data)
        
        # Calculate correlations
        correlation_results = {}
        
        for method in self.correlation_methods:
            correlations = {}
            
            # Correlate genome size with different repeat types
            repeat_columns = [
                'total_repeat_percentage',
                'satellite_percentage',
                'ltr_percentage',
                'line_percentage',
                'sine_percentage',
                'dna_percentage',
                'unknown_percentage'
            ]
            
            for column in repeat_columns:
                if column in df.columns:
                    # Remove NaN values
                    clean_data = df[['genome_size_gb', column]].dropna()
                    
                    self.logger.info(f"Column: {column}, Data points: {len(clean_data)}")
                    if len(clean_data) >= 2:  # Reduced minimum sample size for testing
                        if method == 'pearson':
                            corr, p_value = pearsonr(clean_data['genome_size_gb'], clean_data[column])
                        else:  # spearman
                            corr, p_value = spearmanr(clean_data['genome_size_gb'], clean_data[column])
                        
                        correlations[column] = {
                            'correlation': corr,
                            'p_value': p_value,
                            'significant': p_value < self.significance_level,
                            'sample_count': len(clean_data)
                        }
                        self.logger.info(f"Correlation for {column}: {corr:.3f}, p={p_value:.3f}")
            
            correlation_results[method] = correlations
        
        self.correlation_results = correlation_results
        self.logger.info(f"Calculated correlations for {len(combined_data)} samples")
        self.logger.info(f"Combined data: {combined_data}")
        
        self.correlation_results = correlation_results
        self.combined_data = combined_data
        
        return {
            'correlation_results': correlation_results,
            'combined_data': combined_data,
            'sample_count': len(combined_data),
            'dataframe': df
        }
    
    def generate_correlation_plots(self, output_dir: str = None) -> str:
        """
        Generate correlation plots and reports
        
        Args:
            output_dir: Output directory for plots and reports
            
        Returns:
            Path to generated plots directory
        """
        if not self.correlation_results:
            self.logger.error("No correlation results available")
            return ""
        
        # Use path utilities for output directory
        import sys
        sys.path.append('.')
        from post_tarean.path_utils import get_project_path, ensure_output_dir
        
        if output_dir is None:
            output_dir = get_project_path(self.project_id, "cytometry") / "validation"
        else:
            output_dir = Path(output_dir)
        
        ensure_output_dir(output_dir)
        
        # Create plots
        plots_created = []
        
        # 1. Genome size vs Total repeat percentage
        if 'total_repeat_percentage' in self.correlation_results.get('pearson', {}):
            fig, ax = plt.subplots(figsize=self.figure_size)
            
            df = pd.DataFrame(self.combined_data)
            
            # Create scatter plot
            ax.scatter(df['genome_size_gb'], df['total_repeat_percentage'], 
                      alpha=0.7, s=100, edgecolors='black')
            
            # Add trend line
            z = np.polyfit(df['genome_size_gb'], df['total_repeat_percentage'], 1)
            p = np.poly1d(z)
            ax.plot(df['genome_size_gb'], p(df['genome_size_gb']), "r--", alpha=0.8)
            
            # Add correlation info
            corr_info = self.correlation_results['pearson']['total_repeat_percentage']
            ax.text(0.05, 0.95, f"r = {corr_info['correlation']:.3f}\np = {corr_info['p_value']:.3f}", 
                   transform=ax.transAxes, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
            
            ax.set_xlabel('Genome Size (GB)')
            ax.set_ylabel('Total Repeat Percentage (%)')
            ax.set_title('Genome Size vs Total Repeat Content')
            ax.grid(True, alpha=0.3)
            
            plot_path = output_dir / "genome_size_vs_total_repeats.png"
            plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
            plt.close()
            plots_created.append(plot_path)
        
        # 2. Repeat type breakdown correlation
        repeat_types = ['satellite_percentage', 'ltr_percentage', 'line_percentage', 
                       'sine_percentage', 'dna_percentage', 'unknown_percentage']
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, repeat_type in enumerate(repeat_types):
            if repeat_type in self.correlation_results.get('pearson', {}):
                ax = axes[i]
                
                df = pd.DataFrame(self.combined_data)
                
                # Create scatter plot
                ax.scatter(df['genome_size_gb'], df[repeat_type], 
                          alpha=0.7, s=80, edgecolors='black')
                
                # Add trend line
                z = np.polyfit(df['genome_size_gb'], df[repeat_type], 1)
                p = np.poly1d(z)
                ax.plot(df['genome_size_gb'], p(df['genome_size_gb']), "r--", alpha=0.8)
                
                # Add correlation info
                corr_info = self.correlation_results['pearson'][repeat_type]
                ax.text(0.05, 0.95, f"r = {corr_info['correlation']:.3f}\np = {corr_info['p_value']:.3f}", 
                       transform=ax.transAxes, fontsize=8,
                       bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
                
                ax.set_xlabel('Genome Size (GB)')
                ax.set_ylabel(f'{repeat_type.replace("_", " ").title()} (%)')
                ax.set_title(f'Genome Size vs {repeat_type.replace("_", " ").title()}')
                ax.grid(True, alpha=0.3)
            else:
                axes[i].text(0.5, 0.5, 'No data', ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(f'{repeat_type.replace("_", " ").title()}')
        
        plt.tight_layout()
        plot_path = output_dir / "repeat_type_correlations.png"
        plt.savefig(plot_path, dpi=self.dpi, bbox_inches='tight')
        plt.close()
        plots_created.append(plot_path)
        
        # 3. Generate correlation report
        report_path = output_dir / "cytometry_correlation_report.txt"
        with open(report_path, 'w') as f:
            f.write("Flow Cytometry Correlation Report\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Project: {self.project_id}\n")
            f.write(f"Sample count: {len(self.combined_data)}\n\n")
            
            f.write("Correlation Results:\n")
            f.write("-" * 30 + "\n")
            
            for method, correlations in self.correlation_results.items():
                f.write(f"\n{method.upper()} Correlations:\n")
                for repeat_type, stats in correlations.items():
                    f.write(f"  {repeat_type.replace('_', ' ').title()}:\n")
                    f.write(f"    Correlation: {stats['correlation']:.3f}\n")
                    f.write(f"    P-value: {stats['p_value']:.3f}\n")
                    f.write(f"    Significant: {stats['significant']}\n")
                    f.write(f"    Sample count: {stats['sample_count']}\n\n")
        
        plots_created.append(report_path)
        
        self.logger.info(f"Generated {len(plots_created)} plots and reports in {output_dir}")
        return str(output_dir)
    
    def validate_genome_size_predictions(self, predicted_sizes: Dict[str, float]) -> Dict:
        """
        Validate genome size predictions against flow cytometry measurements
        
        Args:
            predicted_sizes: Dictionary of sample_id -> predicted genome size
            
        Returns:
            Dictionary with validation results
        """
        self.logger.info("Validating genome size predictions")
        
        validation_results = {
            'predictions': {},
            'cytometry_measurements': {},
            'validation_metrics': {}
        }
        
        # Match predictions with cytometry data
        for sample_id, predicted_size in predicted_sizes.items():
            # Find matching cytometry measurement
            cytometry_match = None
            for cytometry in self.cytometry_data:
                if cytometry.sample_id == sample_id:
                    cytometry_match = cytometry
                    break
            
            if cytometry_match:
                error = abs(predicted_size - cytometry_match.genome_size_gb)
                percent_error = (error / cytometry_match.genome_size_gb) * 100
                
                validation_results['predictions'][sample_id] = {
                    'predicted_size': predicted_size,
                    'measured_size': cytometry_match.genome_size_gb,
                    'absolute_error': error,
                    'percent_error': percent_error,
                    'within_10_percent': percent_error <= 10,
                    'within_20_percent': percent_error <= 20
                }
        
        # Calculate validation metrics
        if validation_results['predictions']:
            errors = [pred['percent_error'] for pred in validation_results['predictions'].values()]
            
            validation_results['validation_metrics'] = {
                'mean_error': np.mean(errors),
                'median_error': np.median(errors),
                'std_error': np.std(errors),
                'within_10_percent': sum(1 for e in errors if e <= 10) / len(errors),
                'within_20_percent': sum(1 for e in errors if e <= 20) / len(errors),
                'total_predictions': len(errors)
            }
        
        return validation_results

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description="Flow cytometry correlation analysis")
    parser.add_argument("cytometry_file", help="Path to cytometry data CSV file")
    parser.add_argument("--tarean-dirs", nargs="+", help="TAREAN analysis directories")
    parser.add_argument("--project-id", help="Project identifier")
    parser.add_argument("--output-dir", help="Output directory for plots and reports")
    
    args = parser.parse_args()
    
    # Initialize correlator
    correlator = FlowCytometryCorrelator(args.project_id)
    
    # Load data
    correlator.load_cytometry_data(args.cytometry_file)
    
    if args.tarean_dirs:
        correlator.load_repeat_data_from_tarean(args.tarean_dirs)
    
    # Perform correlation analysis
    results = correlator.correlate_genome_size_repeats()
    
    if results:
        print(f"Correlation analysis completed for {results['sample_count']} samples")
        
        # Generate plots
        output_dir = correlator.generate_correlation_plots(args.output_dir)
        print(f"Plots and reports generated in: {output_dir}")
        
        # Print summary
        for method, correlations in results['correlation_results'].items():
            print(f"\n{method.upper()} Correlations:")
            for repeat_type, stats in correlations.items():
                print(f"  {repeat_type}: r={stats['correlation']:.3f}, p={stats['p_value']:.3f}")
    else:
        print("No correlation results generated")

if __name__ == "__main__":
    main()
