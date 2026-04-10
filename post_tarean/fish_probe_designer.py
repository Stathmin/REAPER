#!/usr/bin/env python3
"""
Optimized FISH Probe Designer for RepOrtR

This module designs FISH probes from TAREAN consensus sequences using BLAST
for fast cross-sample comparison and alignment.

Author: RepOrtR Team
Date: 2025
"""

import os
import sys
import argparse
import json
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
import pandas as pd

# BioPython imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Blast.Applications import NcbiblastnCommandline

# Import existing BLAST analyzer
from blast_consolidated import BLASTAnalyzer

# Define GC function for compatibility
def GC(sequence):
    """Calculate GC content as percentage"""
    return gc_fraction(sequence) * 100

@dataclass
class ConsensusSequence:
    """Represents a TAREAN consensus sequence with metadata"""
    cluster_id: str
    sequence: str
    length: int
    rank: int  # TAREAN rank (1=high confidence, 2=low confidence, etc.)
    gc_content: float
    project: str
    sample: str
    
    def __post_init__(self):
        if not hasattr(self, 'gc_content'):
            self.gc_content = GC(self.sequence)

@dataclass
class FISHProbe:
    """Represents a designed FISH probe"""
    sequence: str
    length: int
    gc_content: float
    melting_temp: float
    specificity_score: float
    probe_type: str  # 'specific', 'general', 'control'
    target_consensus: str
    target_cluster: str
    cross_reactivity: List[str] = None
    blast_hits: List[Dict] = None
    
    def __post_init__(self):
        if self.cross_reactivity is None:
            self.cross_reactivity = []
        if self.blast_hits is None:
            self.blast_hits = []

class OptimizedFISHProbeDesigner:
    """
    Optimized FISH probe designer using BLAST for fast cross-sample comparison
    
    This class analyzes consensus sequences from TAREAN analysis and designs
    FISH probes for experimental validation using BLAST for fast alignment.
    """
    
    def __init__(self, tarean_output_dir: str, project_id: str = None):
        """
        Initialize the optimized FISH probe designer
        
        Args:
            tarean_output_dir: Directory containing TAREAN analysis results
            project_id: Project identifier for output organization
        """
        self.tarean_dir = Path(tarean_output_dir)
        self.project_id = project_id or "default"
        self.consensus_sequences = []
        self.probes = []
        self.blast_analyzer = BLASTAnalyzer()
        
        # Probe design parameters
        self.probe_params = {
            'min_length': 20,
            'max_length': 30,
            'opt_length': 25,
            'min_gc': 40.0,
            'max_gc': 60.0,
            'opt_gc': 50.0,
            'min_tm': 55.0,
            'max_tm': 65.0,
            'opt_tm': 60.0,
            'max_poly_x': 4,  # Maximum homopolymer length
            'max_repeat': 3,   # Maximum repeat length
        }
        
        # BLAST parameters for cross-sample comparison
        self.blast_params = {
            'evalue_threshold': 1e-5,
            'identity_threshold': 80.0,
            'coverage_threshold': 80.0,
            'word_size': 6,
            'num_threads': 8,
        }
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
    
    def load_consensus_sequences(self) -> List[ConsensusSequence]:
        """
        Load consensus sequences from TAREAN output files
        
        Returns:
            List of ConsensusSequence objects
        """
        self.logger.info(f"Loading consensus sequences from {self.tarean_dir}")
        
        consensus_files = {
            1: "TAREAN_consensus_rank_1.fasta",  # High confidence satellites
            2: "TAREAN_consensus_rank_2.fasta",  # Low confidence satellites
            3: "TAREAN_consensus_rank_3.fasta",  # Putative LTR elements
            4: "TAREAN_consensus_rank_4.fasta",  # rDNA
        }
        
        loaded_sequences = []
        
        for rank, filename in consensus_files.items():
            filepath = self.tarean_dir / filename
            if filepath.exists():
                self.logger.info(f"Loading {filename} (rank {rank})")
                
                for record in SeqIO.parse(filepath, "fasta"):
                    # Parse cluster ID from header (e.g., "CL1_TR_1_x_568nt")
                    cluster_id = record.id.split('_')[0]
                    
                    consensus = ConsensusSequence(
                        cluster_id=cluster_id,
                        sequence=str(record.seq),
                        length=len(record.seq),
                        rank=rank,
                        gc_content=GC(str(record.seq)),
                        project=self.project_id,
                        sample=self.tarean_dir.name
                    )
                    
                    loaded_sequences.append(consensus)
                    self.logger.debug(f"Loaded {cluster_id}: {len(record.seq)}bp, GC: {consensus.gc_content:.1f}%")
        
        self.consensus_sequences = loaded_sequences
        self.logger.info(f"Loaded {len(loaded_sequences)} consensus sequences")
        return loaded_sequences
    
    def create_consensus_database(self, output_dir: str = None) -> str:
        """
        Create a BLAST database from consensus sequences for cross-sample comparison
        
        Args:
            output_dir: Output directory for database
            
        Returns:
            Path to created database
        """
        if output_dir is None:
            output_dir = Path(f"projects/{self.project_id}/samples/{self.tarean_dir.name}/blast_db")
        else:
            output_dir = Path(output_dir)
        
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create FASTA file with consensus sequences
        fasta_file = output_dir / "consensus_sequences.fasta"
        with open(fasta_file, 'w') as f:
            for consensus in self.consensus_sequences:
                f.write(f">{consensus.cluster_id}_{consensus.sample}\n")
                f.write(f"{consensus.sequence}\n")
        
        # Create BLAST database
        db_file = output_dir / "consensus_db"
        try:
            cmd = f"makeblastdb -in {fasta_file} -dbtype nucl -out {db_file}"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            self.logger.info(f"Created BLAST database: {db_file}")
            return str(db_file)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create BLAST database: {e}")
            return None
    
    def compare_samples_with_blast(self, other_sample_dirs: List[str]) -> Dict:
        """
        Compare consensus sequences across samples using BLAST
        
        Args:
            other_sample_dirs: List of paths to other sample TAREAN directories
            
        Returns:
            Dictionary with cross-sample comparison results
        """
        self.logger.info("Comparing consensus sequences across samples using BLAST")
        
        if not self.consensus_sequences:
            self.load_consensus_sequences()
        
        # Create database from current sample
        current_db = self.create_consensus_database()
        if not current_db:
            self.logger.error("Failed to create BLAST database for current sample")
            return {}
        
        comparison_results = {
            'samples_compared': 1 + len(other_sample_dirs),
            'total_consensus': len(self.consensus_sequences),
            'shared_motifs': {},
            'sample_specific': {},
            'conserved_families': {},
            'blast_results': {}
        }
        
        # Compare with each other sample
        for sample_dir in other_sample_dirs:
            sample_path = Path(sample_dir)
            if sample_path.exists():
                self.logger.info(f"Comparing with sample: {sample_path.name}")
                
                # Load consensus from other sample
                other_designer = OptimizedFISHProbeDesigner(sample_dir, sample_path.name)
                other_consensus = other_designer.load_consensus_sequences()
                
                # Create FASTA file for other sample
                other_fasta = sample_path / "consensus_for_blast.fasta"
                with open(other_fasta, 'w') as f:
                    for consensus in other_consensus:
                        f.write(f">{consensus.cluster_id}_{consensus.sample}\n")
                        f.write(f"{consensus.sequence}\n")
                
                # Run BLAST comparison
                blast_results = self._run_blast_comparison(str(other_fasta), current_db, sample_path.name)
                comparison_results['blast_results'][sample_path.name] = blast_results
                
                # Analyze BLAST results
                shared_motifs = self._analyze_blast_results(blast_results, sample_path.name)
                comparison_results['shared_motifs'].update(shared_motifs)
                
                # Update total consensus count
                comparison_results['total_consensus'] += len(other_consensus)
        
        self.logger.info(f"Cross-sample comparison complete: {comparison_results['samples_compared']} samples")
        return comparison_results
    
    def _run_blast_comparison(self, query_fasta: str, db_path: str, sample_name: str) -> pd.DataFrame:
        """
        Run BLAST comparison between query sequences and database
        
        Args:
            query_fasta: Path to query FASTA file
            db_path: Path to BLAST database
            sample_name: Name of the sample being compared
            
        Returns:
            DataFrame with BLAST results
        """
        self.logger.info(f"Running BLAST comparison for {sample_name}")
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.blast', delete=False) as tmp_file:
            blast_output = tmp_file.name
        
        try:
            # Run BLAST
            blast_cmd = NcbiblastnCommandline(
                query=query_fasta,
                db=db_path,
                evalue=self.blast_params['evalue_threshold'],
                word_size=self.blast_params['word_size'],
                num_threads=self.blast_params['num_threads'],
                outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
                out=blast_output
            )
            
            blast_cmd()
            
            # Parse BLAST results
            blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                           'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            
            df = pd.read_csv(blast_output, sep='\t', names=blast_columns)
            
            # Filter results based on thresholds
            filtered_df = df[
                (df['pident'] >= self.blast_params['identity_threshold']) &
                (df['evalue'] <= self.blast_params['evalue_threshold'])
            ]
            
            self.logger.info(f"BLAST comparison for {sample_name}: {len(filtered_df)} significant hits")
            return filtered_df
            
        except Exception as e:
            self.logger.error(f"BLAST comparison failed for {sample_name}: {e}")
            return pd.DataFrame()
        
        finally:
            # Clean up temporary file
            if os.path.exists(blast_output):
                os.remove(blast_output)
    
    def _analyze_blast_results(self, blast_df: pd.DataFrame, sample_name: str) -> Dict[str, List[str]]:
        """
        Analyze BLAST results to identify shared motifs and conserved families
        
        Args:
            blast_df: DataFrame with BLAST results
            sample_name: Name of the sample
            
        Returns:
            Dictionary with shared motifs and conserved families
        """
        if blast_df.empty:
            return {}
        
        shared_motifs = defaultdict(list)
        
        # Group by query-subject pairs
        for (qseqid, sseqid), group in blast_df.groupby(['qseqid', 'sseqid']):
            # Calculate alignment statistics
            avg_identity = group['pident'].mean()
            avg_length = group['length'].mean()
            total_hits = len(group)
            
            # Consider as shared motif if significant alignment
            if avg_identity >= 80 and avg_length >= 20:
                motif_key = f"{qseqid}_vs_{sseqid}"
                shared_motifs[motif_key] = {
                    'query': qseqid,
                    'subject': sseqid,
                    'avg_identity': avg_identity,
                    'avg_length': avg_length,
                    'total_hits': total_hits,
                    'sample': sample_name
                }
        
        return dict(shared_motifs)
    
    def design_probes_from_blast_results(self, blast_results: Dict) -> List[FISHProbe]:
        """
        Design probes based on BLAST comparison results
        
        Args:
            blast_results: Dictionary with BLAST comparison results
            
        Returns:
            List of FISHProbe objects
        """
        self.logger.info("Designing probes based on BLAST comparison results")
        
        probes = []
        
        # Design specific probes for unique consensus sequences
        specific_probes = self._design_specific_probes_from_blast(blast_results)
        probes.extend(specific_probes)
        
        # Design general probes for shared motifs
        general_probes = self._design_general_probes_from_blast(blast_results)
        probes.extend(general_probes)
        
        self.logger.info(f"Designed {len(probes)} probes from BLAST results")
        return probes
    
    def _design_specific_probes_from_blast(self, blast_results: Dict) -> List[FISHProbe]:
        """Design specific probes for consensus sequences with low cross-reactivity"""
        specific_probes = []
        
        # Find consensus sequences with few or no BLAST hits (specific to sample)
        for consensus in self.consensus_sequences:
            cluster_id = consensus.cluster_id
            blast_hits = []
            
            # Count BLAST hits for this consensus
            for sample_name, sample_results in blast_results.get('blast_results', {}).items():
                if not sample_results.empty:
                    hits = sample_results[
                        (sample_results['qseqid'].str.contains(cluster_id)) |
                        (sample_results['sseqid'].str.contains(cluster_id))
                    ]
                    if not hits.empty:
                        blast_hits.extend(hits.to_dict('records'))
            
            # If few hits, design specific probes
            if len(blast_hits) <= 2:  # Low cross-reactivity
                probes = self._design_probes_for_consensus(consensus, "specific", blast_hits)
                specific_probes.extend(probes)
        
        return specific_probes
    
    def _design_general_probes_from_blast(self, blast_results: Dict) -> List[FISHProbe]:
        """Design general probes for shared motifs identified by BLAST"""
        general_probes = []
        
        # Extract shared motifs from BLAST results
        shared_motifs = blast_results.get('shared_motifs', {})
        
        for motif_key, motif_info in shared_motifs.items():
            if motif_info['avg_identity'] >= 85 and motif_info['avg_length'] >= 25:
                # Design probe for conserved motif
                probe = self._design_probe_from_motif_info(motif_info, "general")
                if probe:
                    general_probes.append(probe)
        
        return general_probes
    
    def _design_probes_for_consensus(self, consensus: ConsensusSequence, 
                                   probe_type: str, blast_hits: List[Dict]) -> List[FISHProbe]:
        """Design probes for a specific consensus sequence"""
        probes = []
        sequence = consensus.sequence
        window_size = self.probe_params['opt_length']
        
        # Slide window through sequence
        for i in range(0, len(sequence) - window_size + 1, window_size // 2):
            probe_seq = sequence[i:i + window_size]
            
            # Check basic criteria
            if self._is_valid_probe(probe_seq):
                gc = GC(probe_seq)
                tm = self._estimate_melting_temp(probe_seq)
                
                # Calculate specificity based on BLAST hits
                specificity_score = self._calculate_specificity_from_blast(probe_seq, blast_hits)
                
                probe = FISHProbe(
                    sequence=probe_seq,
                    length=len(probe_seq),
                    gc_content=gc,
                    melting_temp=tm,
                    specificity_score=specificity_score,
                    probe_type=probe_type,
                    target_consensus=consensus.cluster_id,
                    target_cluster=consensus.cluster_id,
                    blast_hits=blast_hits
                )
                probes.append(probe)
        
        return probes
    
    def _design_probe_from_motif_info(self, motif_info: Dict, probe_type: str) -> Optional[FISHProbe]:
        """Design a probe from shared motif information"""
        # This is a simplified version - in practice would need sequence context
        # For now, create a representative probe
        probe_seq = "N" * self.probe_params['opt_length']  # Placeholder
        
        if self._is_valid_probe(probe_seq):
            return FISHProbe(
                sequence=probe_seq,
                length=len(probe_seq),
                gc_content=GC(probe_seq),
                melting_temp=self._estimate_melting_temp(probe_seq),
                specificity_score=0.5,  # General probes have lower specificity
                probe_type=probe_type,
                target_consensus="shared_motif",
                target_cluster=motif_info['query']
            )
        
        return None
    
    def _calculate_specificity_from_blast(self, probe_seq: str, blast_hits: List[Dict]) -> float:
        """Calculate specificity score based on BLAST hits"""
        if not blast_hits:
            return 1.0  # No hits = high specificity
        
        # Count hits and calculate specificity
        hit_count = len(blast_hits)
        specificity = max(0.0, 1.0 - (hit_count * 0.1))  # Penalize for each hit
        
        return specificity
    
    def _is_valid_probe(self, sequence: str) -> bool:
        """Check if a sequence meets basic probe criteria"""
        gc = GC(sequence)
        
        # Check GC content
        if gc < self.probe_params['min_gc'] or gc > self.probe_params['max_gc']:
            return False
        
        # Check for homopolymers
        for base in 'ATCG':
            if base * self.probe_params['max_poly_x'] in sequence:
                return False
        
        # Check for simple repeats
        for i in range(2, self.probe_params['max_repeat'] + 1):
            for j in range(len(sequence) - i + 1):
                repeat = sequence[j:j+i]
                if sequence.count(repeat) > 2:
                    return False
        
        return True
    
    def _estimate_melting_temp(self, sequence: str) -> float:
        """Estimate melting temperature using Wallace rule"""
        # Wallace rule: Tm = 2°C(A+T) + 4°C(G+C)
        at_count = sequence.count('A') + sequence.count('T')
        gc_count = sequence.count('G') + sequence.count('C')
        
        return 2 * at_count + 4 * gc_count
    
    def generate_optimized_report(self, blast_results: Dict, output_dir: str = None) -> str:
        """
        Generate comprehensive probe design report based on BLAST results
        
        Args:
            blast_results: Dictionary with BLAST comparison results
            output_dir: Output directory for reports
            
        Returns:
            Path to generated report
        """
        self.logger.info("Generating optimized probe design report")
        
        # Use path utilities for output directory
        import sys
        sys.path.append('.')
        from post_tarean.path_utils import get_project_path, ensure_output_dir
        
        if output_dir is None:
            output_dir = get_project_path(self.project_id, "fish_probe") / "validation"
        else:
            output_dir = Path(output_dir)
        
        ensure_output_dir(output_dir)
        
        # Design probes based on BLAST results
        probes = self.design_probes_from_blast_results(blast_results)
        
        # Separate specific and general probes
        specific_probes = [p for p in probes if p.probe_type == "specific"]
        general_probes = [p for p in probes if p.probe_type == "general"]
        
        # Generate report files
        report_files = []
        
        # 1. Specific probes FASTA
        specific_fasta = output_dir / "specific_probes_blast.fasta"
        with open(specific_fasta, 'w') as f:
            for i, probe in enumerate(specific_probes):
                f.write(f">specific_probe_{i+1}_cluster_{probe.target_cluster}\n")
                f.write(f"{probe.sequence}\n")
        report_files.append(specific_fasta)
        
        # 2. General probes FASTA
        general_fasta = output_dir / "general_probes_blast.fasta"
        with open(general_fasta, 'w') as f:
            for i, probe in enumerate(general_probes):
                f.write(f">general_probe_{i+1}_motif_{probe.target_cluster}\n")
                f.write(f"{probe.sequence}\n")
        report_files.append(general_fasta)
        
        # 3. Probe specifications CSV with BLAST information
        specs_csv = output_dir / "probe_specifications_blast.csv"
        with open(specs_csv, 'w') as f:
            f.write("Probe_ID,Sequence,Length,GC_Content,Melting_Temp,Specificity_Score,Type,Target,BLAST_Hits\n")
            for i, probe in enumerate(specific_probes):
                blast_hit_count = len(probe.blast_hits) if probe.blast_hits else 0
                f.write(f"specific_{i+1},{probe.sequence},{probe.length},{probe.gc_content:.1f},"
                       f"{probe.melting_temp:.1f},{probe.specificity_score:.3f},{probe.probe_type},"
                       f"{probe.target_cluster},{blast_hit_count}\n")
            for i, probe in enumerate(general_probes):
                blast_hit_count = len(probe.blast_hits) if probe.blast_hits else 0
                f.write(f"general_{i+1},{probe.sequence},{probe.length},{probe.gc_content:.1f},"
                       f"{probe.melting_temp:.1f},{probe.specificity_score:.3f},{probe.probe_type},"
                       f"{probe.target_cluster},{blast_hit_count}\n")
        report_files.append(specs_csv)
        
        # 4. BLAST analysis summary
        blast_summary = output_dir / "blast_analysis_summary.txt"
        with open(blast_summary, 'w') as f:
            f.write("BLAST-Based Probe Design Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Project: {self.project_id}\n")
            f.write(f"Sample: {self.tarean_dir.name}\n")
            f.write(f"Consensus sequences analyzed: {len(self.consensus_sequences)}\n\n")
            
            f.write("BLAST Comparison Results:\n")
            f.write(f"- Samples compared: {blast_results.get('samples_compared', 0)}\n")
            f.write(f"- Total consensus: {blast_results.get('total_consensus', 0)}\n")
            f.write(f"- Shared motifs: {len(blast_results.get('shared_motifs', {}))}\n\n")
            
            f.write("Probe Design Results:\n")
            f.write(f"- Specific probes: {len(specific_probes)}\n")
            f.write(f"- General probes: {len(general_probes)}\n")
            f.write(f"- Total probes: {len(probes)}\n\n")
            
            # Add BLAST statistics for each sample
            for sample_name, sample_results in blast_results.get('blast_results', {}).items():
                if not sample_results.empty:
                    f.write(f"BLAST Results for {sample_name}:\n")
                    f.write(f"- Significant hits: {len(sample_results)}\n")
                    f.write(f"- Average identity: {sample_results['pident'].mean():.1f}%\n")
                    f.write(f"- Average length: {sample_results['length'].mean():.1f}bp\n\n")
        
        report_files.append(blast_summary)
        
        self.logger.info(f"Generated optimized probe design report in {output_dir}")
        return str(output_dir)

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description="Optimized FISH probe designer using BLAST")
    parser.add_argument("tarean_dir", help="Directory containing TAREAN analysis results")
    parser.add_argument("--project-id", help="Project identifier")
    parser.add_argument("--output-dir", help="Output directory for reports")
    parser.add_argument("--compare-samples", nargs="+", help="Other sample directories for comparison")
    
    args = parser.parse_args()
    
    # Initialize designer
    designer = OptimizedFISHProbeDesigner(args.tarean_dir, args.project_id)
    
    # Load consensus sequences
    designer.load_consensus_sequences()
    
    # Compare across samples if specified
    if args.compare_samples:
        blast_results = designer.compare_samples_with_blast(args.compare_samples)
        if blast_results:
            print(f"BLAST comparison: {blast_results.get('samples_compared', 0)} samples")
            print(f"Shared motifs: {len(blast_results.get('shared_motifs', {}))}")
        else:
            print("BLAST comparison failed")
            blast_results = {'blast_results': {}}
    else:
        blast_results = {'blast_results': {}}
    
    # Generate optimized report
    output_dir = designer.generate_optimized_report(blast_results, args.output_dir)
    print(f"Optimized probe design report generated in: {output_dir}")

if __name__ == "__main__":
    main()
