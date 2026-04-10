#!/usr/bin/env python3
"""
NCBI Sequence Validation Script for RepOrtR

This script validates gathered NCBI sequences for quality and completeness.

Author: RepOrtR Team
Date: 2025
"""

import os
import sys
import argparse
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

# Define GC function for compatibility
def GC(sequence):
    """Calculate GC content as percentage"""
    return gc_fraction(sequence) * 100

class NCBISequenceValidator:
    """
    NCBI Sequence Validator
    
    This class validates gathered NCBI sequences for quality and completeness.
    """
    
    def __init__(self):
        """Initialize the validator"""
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Quality thresholds
        self.min_length = 100
        self.max_length = 50000
        self.min_gc = 20.0
        self.max_gc = 80.0
        self.max_ambiguous = 0.1  # 10% ambiguous bases
        self.max_homopolymer = 50  # Maximum consecutive same base
    
    def validate_fasta_file(self, fasta_file: str) -> Dict:
        """
        Validate sequences in a FASTA file
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Dictionary with validation results
        """
        self.logger.info(f"Validating FASTA file: {fasta_file}")
        
        if not os.path.exists(fasta_file):
            return {"error": f"File not found: {fasta_file}"}
        
        sequences = []
        issues = []
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_data = {
                    'accession': record.id,
                    'length': len(record.seq),
                    'gc_content': GC(str(record.seq)),
                    'ambiguous_bases': str(record.seq).count('N') + str(record.seq).count('X'),
                    'max_homopolymer': self._calculate_max_homopolymer(str(record.seq)),
                    'issues': []
                }
                
                # Check length
                if seq_data['length'] < self.min_length:
                    seq_data['issues'].append(f"Too short: {seq_data['length']} bp")
                
                if seq_data['length'] > self.max_length:
                    seq_data['issues'].append(f"Too long: {seq_data['length']} bp")
                
                # Check GC content
                if seq_data['gc_content'] < self.min_gc or seq_data['gc_content'] > self.max_gc:
                    seq_data['issues'].append(f"Extreme GC: {seq_data['gc_content']:.1f}%")
                
                # Check ambiguous bases
                ambiguous_ratio = seq_data['ambiguous_bases'] / seq_data['length']
                if ambiguous_ratio > self.max_ambiguous:
                    seq_data['issues'].append(f"Too many ambiguous: {ambiguous_ratio:.1%}")
                
                # Check homopolymers
                if seq_data['max_homopolymer'] > self.max_homopolymer:
                    seq_data['issues'].append(f"Long homopolymer: {seq_data['max_homopolymer']} bases")
                
                sequences.append(seq_data)
                
                if seq_data['issues']:
                    issues.append(seq_data)
        
        except Exception as e:
            return {"error": f"Failed to parse FASTA file: {e}"}
        
        # Calculate statistics
        if sequences:
            lengths = [seq['length'] for seq in sequences]
            gc_contents = [seq['gc_content'] for seq in sequences]
            valid_sequences = [seq for seq in sequences if not seq['issues']]
            
            stats = {
                'total_sequences': len(sequences),
                'valid_sequences': len(valid_sequences),
                'problematic_sequences': len(issues),
                'mean_length': sum(lengths) / len(lengths),
                'min_length': min(lengths),
                'max_length': max(lengths),
                'mean_gc': sum(gc_contents) / len(gc_contents),
                'min_gc': min(gc_contents),
                'max_gc': max(gc_contents),
                'issues': issues
            }
        else:
            stats = {"error": "No sequences found"}
        
        return stats
    
    def validate_metadata_file(self, metadata_file: str) -> Dict:
        """
        Validate metadata CSV file
        
        Args:
            metadata_file: Path to metadata CSV file
            
        Returns:
            Dictionary with validation results
        """
        self.logger.info(f"Validating metadata file: {metadata_file}")
        
        if not os.path.exists(metadata_file):
            return {"error": f"File not found: {metadata_file}"}
        
        try:
            df = pd.read_csv(metadata_file)
            
            # Check required columns
            required_columns = ['accession', 'title', 'length', 'gc_content', 'organism']
            missing_columns = [col for col in required_columns if col not in df.columns]
            
            if missing_columns:
                return {"error": f"Missing required columns: {missing_columns}"}
            
            # Basic statistics
            stats = {
                'total_records': len(df),
                'unique_organisms': df['organism'].nunique(),
                'mean_length': df['length'].mean(),
                'mean_gc': df['gc_content'].mean(),
                'length_range': (df['length'].min(), df['length'].max()),
                'gc_range': (df['gc_content'].min(), df['gc_content'].max()),
                'organisms': df['organism'].value_counts().to_dict()
            }
            
            return stats
            
        except Exception as e:
            return {"error": f"Failed to parse metadata file: {e}"}
    
    def _calculate_max_homopolymer(self, sequence: str) -> int:
        """Calculate maximum homopolymer length"""
        max_homopolymer = 0
        for base in 'ATCG':
            max_homopolymer = max(max_homopolymer, len(max(sequence.split(base), key=len)))
        return max_homopolymer
    
    def generate_validation_report(self, fasta_file: str, metadata_file: str, output_file: str) -> str:
        """
        Generate comprehensive validation report
        
        Args:
            fasta_file: Path to FASTA file
            metadata_file: Path to metadata file
            output_file: Output report file
            
        Returns:
            Path to generated report
        """
        self.logger.info("Generating validation report")
        
        # Validate FASTA file
        fasta_stats = self.validate_fasta_file(fasta_file)
        
        # Validate metadata file
        metadata_stats = self.validate_metadata_file(metadata_file)
        
        # Generate report
        with open(output_file, 'w') as f:
            f.write("NCBI Sequence Validation Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("FASTA File Validation:\n")
            f.write("-" * 25 + "\n")
            if "error" in fasta_stats:
                f.write(f"ERROR: {fasta_stats['error']}\n")
            else:
                f.write(f"Total sequences: {fasta_stats['total_sequences']}\n")
                f.write(f"Valid sequences: {fasta_stats['valid_sequences']}\n")
                f.write(f"Problematic sequences: {fasta_stats['problematic_sequences']}\n")
                f.write(f"Mean length: {fasta_stats['mean_length']:.0f} bp\n")
                f.write(f"Length range: {fasta_stats['min_length']}-{fasta_stats['max_length']} bp\n")
                f.write(f"Mean GC content: {fasta_stats['mean_gc']:.1f}%\n")
                f.write(f"GC range: {fasta_stats['min_gc']:.1f}-{fasta_stats['max_gc']:.1f}%\n\n")
                
                if fasta_stats['issues']:
                    f.write("Issues Found:\n")
                    for seq in fasta_stats['issues']:
                        f.write(f"  {seq['accession']}: {', '.join(seq['issues'])}\n")
                    f.write("\n")
            
            f.write("Metadata File Validation:\n")
            f.write("-" * 28 + "\n")
            if "error" in metadata_stats:
                f.write(f"ERROR: {metadata_stats['error']}\n")
            else:
                f.write(f"Total records: {metadata_stats['total_records']}\n")
                f.write(f"Unique organisms: {metadata_stats['unique_organisms']}\n")
                f.write(f"Mean length: {metadata_stats['mean_length']:.0f} bp\n")
                f.write(f"Mean GC content: {metadata_stats['mean_gc']:.1f}%\n")
                f.write(f"Length range: {metadata_stats['length_range'][0]}-{metadata_stats['length_range'][1]} bp\n")
                f.write(f"GC range: {metadata_stats['gc_range'][0]:.1f}-{metadata_stats['gc_range'][1]:.1f}%\n\n")
                
                f.write("Organism Distribution:\n")
                for organism, count in metadata_stats['organisms'].items():
                    f.write(f"  {organism}: {count} sequences\n")
            
            # Overall assessment
            f.write("\nOverall Assessment:\n")
            f.write("-" * 20 + "\n")
            
            if "error" in fasta_stats or "error" in metadata_stats:
                f.write("❌ VALIDATION FAILED - Check errors above\n")
            elif fasta_stats.get('problematic_sequences', 0) == 0:
                f.write("✅ ALL SEQUENCES PASSED QUALITY CHECKS\n")
            else:
                f.write(f"⚠️  {fasta_stats['problematic_sequences']} SEQUENCES HAVE ISSUES\n")
                f.write("Consider reviewing problematic sequences before use.\n")
        
        self.logger.info(f"Validation report generated: {output_file}")
        return output_file

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description="Validate NCBI gathered sequences")
    parser.add_argument("--fasta", required=True, help="Path to FASTA file")
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV file")
    parser.add_argument("--output", required=True, help="Output validation report file")
    
    args = parser.parse_args()
    
    # Initialize validator
    validator = NCBISequenceValidator()
    
    # Generate validation report
    report_file = validator.generate_validation_report(
        args.fasta, 
        args.metadata, 
        args.output
    )
    
    print(f"Validation report generated: {report_file}")

if __name__ == "__main__":
    main()

