#!/usr/bin/env python3
"""
PCR Validation Module for RepOrtR

This module designs PCR primers for repeat family validation and generates
comprehensive protocols for experimental validation.

Author: RepOrtR Team
Date: 2025
"""

import os
import sys
import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from datetime import datetime
import pandas as pd

# BioPython imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

# Define GC function for compatibility
def GC(sequence):
    """Calculate GC content as percentage"""
    return gc_fraction(sequence) * 100

@dataclass
class PCRPrimer:
    """Represents a designed PCR primer"""
    sequence: str
    length: int
    gc_content: float
    melting_temp: float
    primer_type: str  # 'forward', 'reverse'
    target_sequence: str
    target_family: str
    specificity_score: float
    
    def __post_init__(self):
        if not hasattr(self, 'gc_content'):
            self.gc_content = GC(self.sequence)
        if not hasattr(self, 'melting_temp'):
            self.melting_temp = self._estimate_melting_temp()
    
    def _estimate_melting_temp(self) -> float:
        """Estimate melting temperature using Wallace rule"""
        at_count = self.sequence.count('A') + self.sequence.count('T')
        gc_count = self.sequence.count('G') + self.sequence.count('C')
        return 2 * at_count + 4 * gc_count

@dataclass
class PCRProtocol:
    """Represents a complete PCR protocol"""
    protocol_id: str
    target_family: str
    forward_primer: PCRPrimer
    reverse_primer: PCRPrimer
    expected_product_size: int
    annealing_temp: float
    created_at: str
    
    def __post_init__(self):
        if not hasattr(self, 'created_at'):
            self.created_at = datetime.now().isoformat()

class PCRValidator:
    """
    PCR Validation Module
    
    This class designs PCR primers for repeat family validation and generates
    comprehensive protocols for experimental validation.
    """
    
    def __init__(self, repeat_sequences: List[str] = None, project_id: str = None):
        """
        Initialize the PCR validator
        
        Args:
            repeat_sequences: List of repeat sequences to design primers for
            project_id: Project identifier for output organization
        """
        self.repeat_sequences = repeat_sequences or []
        self.project_id = project_id or "default"
        self.primers = []
        self.protocols = []
        
        # Primer design parameters
        self.primer_params = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_MIN_GC_PERCENT': 40.0,
            'PRIMER_MAX_GC_PERCENT': 60.0,
            'PRIMER_MAX_POLY_X': 4
        }
        
        # PCR protocol parameters
        self.pcr_params = {
            'denaturation_temp': 94.0,
            'denaturation_time': 30,
            'annealing_time': 30,
            'extension_temp': 72.0,
            'extension_time': 60,
            'final_extension_time': 300,
            'cycle_count': 35
        }
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def load_repeat_sequences(self, fasta_file: str) -> List[str]:
        """
        Load repeat sequences from FASTA file
        
        Args:
            fasta_file: Path to FASTA file with repeat sequences
            
        Returns:
            List of repeat sequences
        """
        self.logger.info(f"Loading repeat sequences from {fasta_file}")
        
        sequences = []
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append(str(record.seq))
            
            self.repeat_sequences = sequences
            self.logger.info(f"Loaded {len(sequences)} repeat sequences")
            return sequences
            
        except Exception as e:
            self.logger.error(f"Failed to load repeat sequences: {e}")
            return []
    
    def design_pcr_primers(self, repeat_family: str, sequence: str) -> List[PCRPrimer]:
        """
        Design PCR primers for a repeat family
        
        Args:
            repeat_family: Name of the repeat family
            sequence: Repeat sequence to design primers for
            
        Returns:
            List of designed PCR primers
        """
        self.logger.info(f"Designing PCR primers for {repeat_family}")
        
        primers = []
        window_size = self.primer_params['PRIMER_OPT_SIZE']
        
        # Design forward primers (5' end)
        for i in range(0, len(sequence) - window_size + 1, window_size // 2):
            primer_seq = sequence[i:i + window_size]
            
            if self._is_valid_primer_sequence(primer_seq):
                forward_primer = PCRPrimer(
                    sequence=primer_seq,
                    length=len(primer_seq),
                    gc_content=GC(primer_seq),
                    melting_temp=2 * (primer_seq.count('A') + primer_seq.count('T')) + 4 * (primer_seq.count('G') + primer_seq.count('C')),
                    primer_type='forward',
                    target_sequence=sequence,
                    target_family=repeat_family,
                    specificity_score=0.8
                )
                primers.append(forward_primer)
        
        # Design reverse primers (3' end)
        for i in range(window_size, len(sequence) + 1, window_size // 2):
            primer_seq = sequence[i-window_size:i]
            
            if self._is_valid_primer_sequence(primer_seq):
                reverse_primer = PCRPrimer(
                    sequence=primer_seq,
                    length=len(primer_seq),
                    gc_content=GC(primer_seq),
                    melting_temp=2 * (primer_seq.count('A') + primer_seq.count('T')) + 4 * (primer_seq.count('G') + primer_seq.count('C')),
                    primer_type='reverse',
                    target_sequence=sequence,
                    target_family=repeat_family,
                    specificity_score=0.8
                )
                primers.append(reverse_primer)
        
        self.logger.info(f"Designed {len(primers)} primers for {repeat_family}")
        return primers
    
    def _is_valid_primer_sequence(self, sequence: str) -> bool:
        """Check if a sequence is valid for primer design"""
        if len(sequence) < self.primer_params['PRIMER_MIN_SIZE']:
            return False
        
        gc_content = GC(sequence)
        if gc_content < self.primer_params['PRIMER_MIN_GC_PERCENT'] or \
           gc_content > self.primer_params['PRIMER_MAX_GC_PERCENT']:
            return False
        
        # Check for homopolymers
        for base in 'ATCG':
            if base * self.primer_params['PRIMER_MAX_POLY_X'] in sequence:
                return False
        
        # Check for ambiguous bases
        if 'N' in sequence:
            return False
        
        return True
    
    def create_pcr_protocol(self, forward_primer: PCRPrimer, reverse_primer: PCRPrimer) -> PCRProtocol:
        """
        Create a PCR protocol for a primer pair
        
        Args:
            forward_primer: Forward primer
            reverse_primer: Reverse primer
            
        Returns:
            PCRProtocol object
        """
        # Calculate expected product size
        target_seq = forward_primer.target_sequence
        expected_size = len(target_seq)
        
        # Calculate annealing temperature
        annealing_temp = (forward_primer.melting_temp + reverse_primer.melting_temp) / 2 - 5
        
        # Ensure annealing temperature is reasonable
        if annealing_temp < 50:
            annealing_temp = 55
        elif annealing_temp > 70:
            annealing_temp = 65
        
        protocol = PCRProtocol(
            protocol_id=f"PCR_{forward_primer.target_family}_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
            target_family=forward_primer.target_family,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer,
            expected_product_size=expected_size,
            annealing_temp=annealing_temp,
            created_at=datetime.now().isoformat()
        )
        
        return protocol
    
    def generate_pcr_protocols(self, output_dir: str = None) -> str:
        """
        Generate PCR protocols for all repeat families
        
        Args:
            output_dir: Output directory for protocols
            
        Returns:
            Path to generated protocols directory
        """
        self.logger.info("Generating PCR protocols for repeat families")
        
        # Use path utilities for output directory
        import sys
        sys.path.append('.')
        from post_tarean.path_utils import get_project_path, ensure_output_dir
        
        if output_dir is None:
            output_dir = get_project_path(self.project_id, "pcr_validation") / "validation"
        else:
            output_dir = Path(output_dir)
        
        ensure_output_dir(output_dir)
        
        protocols_created = []
        
        # Generate protocols for each repeat sequence
        for i, sequence in enumerate(self.repeat_sequences):
            family_name = f"repeat_family_{i+1}"
            
            # Design primers
            primers = self.design_pcr_primers(family_name, sequence)
            
            # Group primers by type
            forward_primers = [p for p in primers if p.primer_type == 'forward']
            reverse_primers = [p for p in primers if p.primer_type == 'reverse']
            
            # Create protocols for primer pairs
            for fw_primer in forward_primers[:2]:  # Limit to 2 forward primers
                for rv_primer in reverse_primers[:2]:  # Limit to 2 reverse primers
                    protocol = self.create_pcr_protocol(fw_primer, rv_primer)
                    self.protocols.append(protocol)
                    
                    # Generate protocol file
                    protocol_file = self._generate_protocol_file(protocol, output_dir)
                    protocols_created.append(protocol_file)
        
        # Generate summary report
        summary_file = self._generate_summary_report(output_dir)
        protocols_created.append(summary_file)
        
        self.logger.info(f"Generated {len(protocols_created)} PCR protocol files in {output_dir}")
        return str(output_dir)
    
    def _generate_protocol_file(self, protocol: PCRProtocol, output_dir: Path) -> Path:
        """Generate a detailed PCR protocol file"""
        protocol_file = output_dir / f"{protocol.protocol_id}_protocol.txt"
        
        with open(protocol_file, 'w') as f:
            f.write("PCR Protocol for Repeat Family Validation\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Protocol ID: {protocol.protocol_id}\n")
            f.write(f"Target Family: {protocol.target_family}\n")
            f.write(f"Created: {protocol.created_at}\n\n")
            
            f.write("Primer Information:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Forward Primer: {protocol.forward_primer.sequence}\n")
            f.write(f"  Length: {protocol.forward_primer.length} bp\n")
            f.write(f"  GC Content: {protocol.forward_primer.gc_content:.1f}%\n")
            f.write(f"  Melting Temp: {protocol.forward_primer.melting_temp:.1f}°C\n")
            f.write(f"  Specificity Score: {protocol.forward_primer.specificity_score:.2f}\n\n")
            
            f.write(f"Reverse Primer: {protocol.reverse_primer.sequence}\n")
            f.write(f"  Length: {protocol.reverse_primer.length} bp\n")
            f.write(f"  GC Content: {protocol.reverse_primer.gc_content:.1f}%\n")
            f.write(f"  Melting Temp: {protocol.reverse_primer.melting_temp:.1f}°C\n")
            f.write(f"  Specificity Score: {protocol.reverse_primer.specificity_score:.2f}\n\n")
            
            f.write("PCR Conditions:\n")
            f.write("-" * 15 + "\n")
            f.write(f"Expected Product Size: {protocol.expected_product_size} bp\n")
            f.write(f"Denaturation: 94°C for 30s\n")
            f.write(f"Annealing: {protocol.annealing_temp:.1f}°C for 30s\n")
            f.write(f"Extension: 72°C for 60s\n")
            f.write(f"Final Extension: 72°C for 5 min\n")
            f.write(f"Cycles: 35\n\n")
            
            f.write("Reaction Setup:\n")
            f.write("-" * 15 + "\n")
            f.write("Component\t\tVolume\t\tFinal Concentration\n")
            f.write("Template DNA\t\t1-5 μL\t\t10-50 ng\n")
            f.write("Forward Primer\t\t1 μL\t\t0.5 μM\n")
            f.write("Reverse Primer\t\t1 μL\t\t0.5 μM\n")
            f.write("dNTPs\t\t\t2 μL\t\t200 μM each\n")
            f.write("10x Buffer\t\t2 μL\t\t1x\n")
            f.write("MgCl2\t\t\t1.5 μL\t\t1.5 mM\n")
            f.write("Taq Polymerase\t\t0.2 μL\t\t1 U\n")
            f.write("H2O\t\t\tto 20 μL\n")
            f.write("Total Volume\t\t20 μL\n")
        
        return protocol_file
    
    def _generate_summary_report(self, output_dir: Path) -> Path:
        """Generate a summary report of all PCR protocols"""
        summary_file = output_dir / "pcr_protocols_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("PCR Protocols Summary Report\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Project: {self.project_id}\n")
            f.write(f"Total Protocols: {len(self.protocols)}\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n\n")
            
            f.write("Protocol Summary:\n")
            f.write("-" * 20 + "\n")
            
            for i, protocol in enumerate(self.protocols, 1):
                f.write(f"{i}. {protocol.protocol_id}\n")
                f.write(f"   Target: {protocol.target_family}\n")
                f.write(f"   Forward: {protocol.forward_primer.sequence}\n")
                f.write(f"   Reverse: {protocol.reverse_primer.sequence}\n")
                f.write(f"   Product Size: {protocol.expected_product_size} bp\n")
                f.write(f"   Annealing Temp: {protocol.annealing_temp:.1f}°C\n")
                f.write(f"   Specificity: {protocol.forward_primer.specificity_score:.2f}/{protocol.reverse_primer.specificity_score:.2f}\n\n")
        
        return summary_file

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description="PCR validation module")
    parser.add_argument("repeat_fasta", help="Path to repeat sequences FASTA file")
    parser.add_argument("--project-id", help="Project identifier")
    parser.add_argument("--output-dir", help="Output directory for protocols")
    
    args = parser.parse_args()
    
    # Initialize PCR validator
    validator = PCRValidator(project_id=args.project_id)
    
    # Load repeat sequences
    sequences = validator.load_repeat_sequences(args.repeat_fasta)
    
    if sequences:
        print(f"Loaded {len(sequences)} repeat sequences")
        
        # Generate PCR protocols
        output_dir = validator.generate_pcr_protocols(args.output_dir)
        print(f"PCR protocols generated in: {output_dir}")
    else:
        print("No repeat sequences loaded")

if __name__ == "__main__":
    main()
