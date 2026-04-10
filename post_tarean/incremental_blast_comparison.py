#!/usr/bin/env python3
"""
Incremental BLAST Comparison System for RepOrtR

This module maintains a persistent database of BLAST comparison results
to avoid re-computing when samples are added or removed from projects.

Author: RepOrtR Team
Date: 2025
"""

import os
import sys
import json
import sqlite3
import hashlib
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, asdict
from datetime import datetime
import pandas as pd

# BioPython imports
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

# Import existing BLAST analyzer
from blast_consolidated import BLASTAnalyzer

@dataclass
class BLASTComparisonRecord:
    """Represents a BLAST comparison result with metadata"""
    comparison_id: str
    project_id: str
    sample_a: str
    sample_b: str
    blast_hits: int
    avg_identity: float
    avg_length: float
    shared_motifs: int
    created_at: str
    updated_at: str
    blast_results_path: str
    status: str  # 'completed', 'failed', 'pending'
    
    def __post_init__(self):
        if not hasattr(self, 'created_at'):
            self.created_at = datetime.now().isoformat()
        if not hasattr(self, 'updated_at'):
            self.updated_at = datetime.now().isoformat()

class IncrementalBLASTComparison:
    """
    Incremental BLAST comparison system with persistent results database
    
    This class maintains a SQLite database of BLAST comparison results
    to avoid re-computing when samples are added or removed.
    """
    
    def __init__(self, project_id: str, database_path: str = None):
        """
        Initialize the incremental BLAST comparison system
        
        Args:
            project_id: Project identifier
            database_path: Path to SQLite database (optional)
        """
        self.project_id = project_id
        
        # Setup database path using path utilities
        import sys
        sys.path.append('.')
        from post_tarean.path_utils import get_project_path, ensure_output_dir
        
        if database_path is None:
            project_dir = get_project_path(project_id, "blast_comparison")
            database_path = project_dir / "blast_comparisons.db"
        self.database_path = Path(database_path)
        ensure_output_dir(self.database_path.parent)
        
        # BLAST parameters
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
        
        # Initialize database
        self._init_database()
    
    def _init_database(self):
        """Initialize SQLite database with required tables"""
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            
            # Create comparisons table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS blast_comparisons (
                    comparison_id TEXT PRIMARY KEY,
                    project_id TEXT NOT NULL,
                    sample_a TEXT NOT NULL,
                    sample_b TEXT NOT NULL,
                    blast_hits INTEGER,
                    avg_identity REAL,
                    avg_length REAL,
                    shared_motifs INTEGER,
                    created_at TEXT,
                    updated_at TEXT,
                    blast_results_path TEXT,
                    status TEXT DEFAULT 'pending',
                    UNIQUE(project_id, sample_a, sample_b)
                )
            ''')
            
            # Create sample metadata table
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS sample_metadata (
                    sample_id TEXT PRIMARY KEY,
                    project_id TEXT NOT NULL,
                    consensus_count INTEGER,
                    total_length INTEGER,
                    created_at TEXT,
                    updated_at TEXT,
                    fasta_hash TEXT
                )
            ''')
            
            # Create indexes for performance
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_project_samples ON blast_comparisons(project_id, sample_a, sample_b)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_sample_metadata ON sample_metadata(sample_id, project_id)')
            
            conn.commit()
    
    def _get_comparison_id(self, sample_a: str, sample_b: str) -> str:
        """Generate unique comparison ID"""
        # Sort samples to ensure consistent ID regardless of order
        sorted_samples = sorted([sample_a, sample_b])
        return f"{self.project_id}_{sorted_samples[0]}_{sorted_samples[1]}"
    
    def _get_sample_hash(self, sample_dir: str) -> str:
        """Calculate hash of sample consensus sequences to detect changes"""
        sample_path = Path(sample_dir)
        consensus_files = [
            "TAREAN_consensus_rank_1.fasta",
            "TAREAN_consensus_rank_2.fasta", 
            "TAREAN_consensus_rank_3.fasta",
            "TAREAN_consensus_rank_4.fasta"
        ]
        
        hasher = hashlib.md5()
        for filename in consensus_files:
            filepath = sample_path / filename
            if filepath.exists():
                hasher.update(filepath.read_bytes())
        
        return hasher.hexdigest()
    
    def _update_sample_metadata(self, sample_id: str, sample_dir: str):
        """Update sample metadata in database"""
        sample_path = Path(sample_dir)
        consensus_count = 0
        total_length = 0
        
        # Count consensus sequences and total length
        consensus_files = [
            "TAREAN_consensus_rank_1.fasta",
            "TAREAN_consensus_rank_2.fasta",
            "TAREAN_consensus_rank_3.fasta", 
            "TAREAN_consensus_rank_4.fasta"
        ]
        
        for filename in consensus_files:
            filepath = sample_path / filename
            if filepath.exists():
                for record in SeqIO.parse(filepath, "fasta"):
                    consensus_count += 1
                    total_length += len(record.seq)
        
        fasta_hash = self._get_sample_hash(sample_dir)
        
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                INSERT OR REPLACE INTO sample_metadata 
                (sample_id, project_id, consensus_count, total_length, created_at, updated_at, fasta_hash)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            ''', (
                sample_id,
                self.project_id,
                consensus_count,
                total_length,
                datetime.now().isoformat(),
                datetime.now().isoformat(),
                fasta_hash
            ))
            conn.commit()
    
    def _check_sample_changed(self, sample_id: str, sample_dir: str) -> bool:
        """Check if sample has changed since last comparison"""
        current_hash = self._get_sample_hash(sample_dir)
        
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT fasta_hash FROM sample_metadata WHERE sample_id = ?', (sample_id,))
            result = cursor.fetchone()
            
            if result is None:
                return True  # New sample
            
            stored_hash = result[0]
            return current_hash != stored_hash
    
    def _get_existing_comparison(self, sample_a: str, sample_b: str) -> Optional[BLASTComparisonRecord]:
        """Get existing comparison from database"""
        comparison_id = self._get_comparison_id(sample_a, sample_b)
        
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT * FROM blast_comparisons 
                WHERE comparison_id = ?
            ''', (comparison_id,))
            
            result = cursor.fetchone()
            if result:
                return BLASTComparisonRecord(
                    comparison_id=result[0],
                    project_id=result[1],
                    sample_a=result[2],
                    sample_b=result[3],
                    blast_hits=result[4],
                    avg_identity=result[5],
                    avg_length=result[6],
                    shared_motifs=result[7],
                    created_at=result[8],
                    updated_at=result[9],
                    blast_results_path=result[10],
                    status=result[11]
                )
        
        return None
    
    def _save_comparison_result(self, record: BLASTComparisonRecord):
        """Save comparison result to database"""
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            cursor.execute('''
                INSERT OR REPLACE INTO blast_comparisons 
                (comparison_id, project_id, sample_a, sample_b, blast_hits, avg_identity, 
                 avg_length, shared_motifs, created_at, updated_at, blast_results_path, status)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                record.comparison_id,
                record.project_id,
                record.sample_a,
                record.sample_b,
                record.blast_hits,
                record.avg_identity,
                record.avg_length,
                record.shared_motifs,
                record.created_at,
                record.updated_at,
                record.blast_results_path,
                record.status
            ))
            conn.commit()
    
    def _run_blast_comparison(self, sample_a_dir: str, sample_b_dir: str) -> Dict:
        """Run BLAST comparison between two samples"""
        sample_a = Path(sample_a_dir).name
        sample_b = Path(sample_b_dir).name
        
        self.logger.info(f"Running BLAST comparison: {sample_a} vs {sample_b}")
        
        # Create BLAST database from sample A
        db_path = self._create_sample_database(sample_a_dir, sample_a)
        if not db_path:
            return {}
        
        # Create FASTA file from sample B
        query_fasta = self._create_sample_fasta(sample_b_dir, sample_b)
        if not query_fasta:
            return {}
        
        # Run BLAST
        blast_results = self._execute_blast(query_fasta, db_path, sample_b)
        
        # Clean up temporary files
        if os.path.exists(query_fasta):
            os.remove(query_fasta)
        
        return blast_results
    
    def _create_sample_database(self, sample_dir: str, sample_name: str) -> Optional[str]:
        """Create BLAST database from sample consensus sequences"""
        sample_path = Path(sample_dir)
        db_dir = Path(f"projects/{self.project_id}/blast_db/{sample_name}")
        db_dir.mkdir(parents=True, exist_ok=True)
        
        # Create FASTA file with consensus sequences
        fasta_file = db_dir / "consensus_sequences.fasta"
        consensus_files = [
            "TAREAN_consensus_rank_1.fasta",
            "TAREAN_consensus_rank_2.fasta",
            "TAREAN_consensus_rank_3.fasta",
            "TAREAN_consensus_rank_4.fasta"
        ]
        
        with open(fasta_file, 'w') as f:
            for filename in consensus_files:
                filepath = sample_path / filename
                if filepath.exists():
                    for record in SeqIO.parse(filepath, "fasta"):
                        f.write(f">{record.id}_{sample_name}\n")
                        f.write(f"{record.seq}\n")
        
        # Create BLAST database
        db_file = db_dir / "consensus_db"
        try:
            cmd = f"makeblastdb -in {fasta_file} -dbtype nucl -out {db_file}"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            return str(db_file)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create BLAST database: {e}")
            return None
    
    def _create_sample_fasta(self, sample_dir: str, sample_name: str) -> Optional[str]:
        """Create FASTA file from sample consensus sequences"""
        sample_path = Path(sample_dir)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_file:
            consensus_files = [
                "TAREAN_consensus_rank_1.fasta",
                "TAREAN_consensus_rank_2.fasta",
                "TAREAN_consensus_rank_3.fasta",
                "TAREAN_consensus_rank_4.fasta"
            ]
            
            for filename in consensus_files:
                filepath = sample_path / filename
                if filepath.exists():
                    for record in SeqIO.parse(filepath, "fasta"):
                        tmp_file.write(f">{record.id}_{sample_name}\n")
                        tmp_file.write(f"{record.seq}\n")
            
            return tmp_file.name
    
    def _execute_blast(self, query_fasta: str, db_path: str, sample_name: str) -> Dict:
        """Execute BLAST search and return results"""
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
            
            # Calculate statistics
            blast_hits = len(filtered_df)
            avg_identity = filtered_df['pident'].mean() if not filtered_df.empty else 0.0
            avg_length = filtered_df['length'].mean() if not filtered_df.empty else 0.0
            
            # Count shared motifs (unique query-subject pairs)
            shared_motifs = len(filtered_df.groupby(['qseqid', 'sseqid'])) if not filtered_df.empty else 0
            
            # Save BLAST results to file
            results_path = f"projects/{self.project_id}/blast_results/{sample_name}_results.csv"
            Path(results_path).parent.mkdir(parents=True, exist_ok=True)
            filtered_df.to_csv(results_path, index=False)
            
            return {
                'blast_hits': blast_hits,
                'avg_identity': avg_identity,
                'avg_length': avg_length,
                'shared_motifs': shared_motifs,
                'blast_results_path': results_path,
                'dataframe': filtered_df
            }
            
        except Exception as e:
            self.logger.error(f"BLAST execution failed: {e}")
            return {}
        
        finally:
            # Clean up temporary file
            if os.path.exists(blast_output):
                os.remove(blast_output)
    
    def compare_samples_incremental(self, sample_dirs: List[str], force_recompute: bool = False) -> Dict:
        """
        Compare samples incrementally, reusing existing results when possible
        
        Args:
            sample_dirs: List of sample directory paths
            force_recompute: Force recomputation of all comparisons
            
        Returns:
            Dictionary with comparison results and statistics
        """
        self.logger.info(f"Starting incremental BLAST comparison for {len(sample_dirs)} samples")
        
        # Update sample metadata
        for sample_dir in sample_dirs:
            sample_id = Path(sample_dir).name
            self._update_sample_metadata(sample_id, sample_dir)
        
        # Generate all pairwise comparisons
        comparisons = []
        for i, sample_a_dir in enumerate(sample_dirs):
            for sample_b_dir in sample_dirs[i+1:]:
                sample_a = Path(sample_a_dir).name
                sample_b = Path(sample_b_dir).name
                comparisons.append((sample_a, sample_b, sample_a_dir, sample_b_dir))
        
        # Track statistics
        total_comparisons = len(comparisons)
        reused_comparisons = 0
        new_comparisons = 0
        failed_comparisons = 0
        
        all_results = {}
        
        for sample_a, sample_b, sample_a_dir, sample_b_dir in comparisons:
            comparison_id = self._get_comparison_id(sample_a, sample_b)
            
            # Check if comparison exists and is still valid
            existing_record = self._get_existing_comparison(sample_a, sample_b)
            
            if existing_record and not force_recompute:
                # Check if samples have changed
                sample_a_changed = self._check_sample_changed(sample_a, sample_a_dir)
                sample_b_changed = self._check_sample_changed(sample_b, sample_b_dir)
                
                if not sample_a_changed and not sample_b_changed:
                    # Reuse existing result
                    self.logger.info(f"Reusing existing comparison: {sample_a} vs {sample_b}")
                    reused_comparisons += 1
                    
                    # Load existing BLAST results
                    if os.path.exists(existing_record.blast_results_path):
                        df = pd.read_csv(existing_record.blast_results_path)
                        all_results[comparison_id] = {
                            'blast_hits': existing_record.blast_hits,
                            'avg_identity': existing_record.avg_identity,
                            'avg_length': existing_record.avg_length,
                            'shared_motifs': existing_record.shared_motifs,
                            'blast_results_path': existing_record.blast_results_path,
                            'dataframe': df,
                            'reused': True
                        }
                        continue
            
            # Run new comparison
            self.logger.info(f"Running new comparison: {sample_a} vs {sample_b}")
            new_comparisons += 1
            
            try:
                blast_results = self._run_blast_comparison(sample_a_dir, sample_b_dir)
                
                if blast_results:
                    # Create record
                    record = BLASTComparisonRecord(
                        comparison_id=comparison_id,
                        project_id=self.project_id,
                        sample_a=sample_a,
                        sample_b=sample_b,
                        blast_hits=blast_results['blast_hits'],
                        avg_identity=blast_results['avg_identity'],
                        avg_length=blast_results['avg_length'],
                        shared_motifs=blast_results['shared_motifs'],
                        blast_results_path=blast_results['blast_results_path'],
                        status='completed',
                        created_at=datetime.now().isoformat(),
                        updated_at=datetime.now().isoformat()
                    )
                    
                    # Save to database
                    self._save_comparison_result(record)
                    
                    # Add to results
                    all_results[comparison_id] = {
                        **blast_results,
                        'reused': False
                    }
                else:
                    failed_comparisons += 1
                    
            except Exception as e:
                self.logger.error(f"Comparison failed: {sample_a} vs {sample_b}: {e}")
                failed_comparisons += 1
        
        # Generate summary
        summary = {
            'project_id': self.project_id,
            'total_comparisons': total_comparisons,
            'reused_comparisons': reused_comparisons,
            'new_comparisons': new_comparisons,
            'failed_comparisons': failed_comparisons,
            'comparison_results': all_results,
            'database_path': str(self.database_path)
        }
        
        self.logger.info(f"Incremental comparison complete:")
        self.logger.info(f"- Total comparisons: {total_comparisons}")
        self.logger.info(f"- Reused: {reused_comparisons}")
        self.logger.info(f"- New: {new_comparisons}")
        self.logger.info(f"- Failed: {failed_comparisons}")
        
        return summary
    
    def add_sample(self, sample_dir: str) -> Dict:
        """
        Add a new sample and compare with all existing samples
        
        Args:
            sample_dir: Path to new sample directory
            
        Returns:
            Dictionary with comparison results
        """
        sample_id = Path(sample_dir).name
        self.logger.info(f"Adding new sample: {sample_id}")
        
        # Get all existing samples from database
        existing_samples = self._get_existing_samples()
        
        # Add new sample to list
        all_samples = existing_samples + [sample_dir]
        
        # Run incremental comparison
        return self.compare_samples_incremental(all_samples)
    
    def remove_sample(self, sample_id: str) -> Dict:
        """
        Remove a sample and clean up related comparisons
        
        Args:
            sample_id: ID of sample to remove
            
        Returns:
            Dictionary with cleanup results
        """
        self.logger.info(f"Removing sample: {sample_id}")
        
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            
            # Remove comparisons involving this sample
            cursor.execute('''
                DELETE FROM blast_comparisons 
                WHERE sample_a = ? OR sample_b = ?
            ''', (sample_id, sample_id))
            
            # Remove sample metadata
            cursor.execute('DELETE FROM sample_metadata WHERE sample_id = ?', (sample_id,))
            
            conn.commit()
        
        return {
            'removed_sample': sample_id,
            'database_path': str(self.database_path)
        }
    
    def _get_existing_samples(self) -> List[str]:
        """Get list of existing sample directories"""
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            cursor.execute('SELECT sample_id FROM sample_metadata WHERE project_id = ?', (self.project_id,))
            results = cursor.fetchall()
            
            # Convert to sample directory paths
            sample_dirs = []
            for (sample_id,) in results:
                sample_dir = f"projects/{self.project_id}/samples/{sample_id}"
                if Path(sample_dir).exists():
                    sample_dirs.append(sample_dir)
            
            return sample_dirs
    
    def get_comparison_statistics(self) -> Dict:
        """Get statistics about stored comparisons"""
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.cursor()
            
            # Get total comparisons
            cursor.execute('SELECT COUNT(*) FROM blast_comparisons WHERE project_id = ?', (self.project_id,))
            total_comparisons = cursor.fetchone()[0]
            
            # Get completed comparisons
            cursor.execute('SELECT COUNT(*) FROM blast_comparisons WHERE project_id = ? AND status = "completed"', (self.project_id,))
            completed_comparisons = cursor.fetchone()[0]
            
            # Get failed comparisons
            cursor.execute('SELECT COUNT(*) FROM blast_comparisons WHERE project_id = ? AND status = "failed"', (self.project_id,))
            failed_comparisons = cursor.fetchone()[0]
            
            # Get total samples
            cursor.execute('SELECT COUNT(*) FROM sample_metadata WHERE project_id = ?', (self.project_id,))
            total_samples = cursor.fetchone()[0]
            
            return {
                'project_id': self.project_id,
                'total_samples': total_samples,
                'total_comparisons': total_comparisons,
                'completed_comparisons': completed_comparisons,
                'failed_comparisons': failed_comparisons,
                'database_path': str(self.database_path)
            }

def main():
    """Main function for command-line usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Incremental BLAST comparison system")
    parser.add_argument("project_id", help="Project identifier")
    parser.add_argument("--sample-dirs", nargs="+", help="Sample directories to compare")
    parser.add_argument("--add-sample", help="Add new sample directory")
    parser.add_argument("--remove-sample", help="Remove sample by ID")
    parser.add_argument("--force-recompute", action="store_true", help="Force recomputation of all comparisons")
    parser.add_argument("--stats", action="store_true", help="Show comparison statistics")
    
    args = parser.parse_args()
    
    # Initialize incremental comparison system
    incremental_blast = IncrementalBLASTComparison(args.project_id)
    
    if args.stats:
        # Show statistics
        stats = incremental_blast.get_comparison_statistics()
        print("Comparison Statistics:")
        print(f"- Project: {stats['project_id']}")
        print(f"- Total samples: {stats['total_samples']}")
        print(f"- Total comparisons: {stats['total_comparisons']}")
        print(f"- Completed: {stats['completed_comparisons']}")
        print(f"- Failed: {stats['failed_comparisons']}")
        print(f"- Database: {stats['database_path']}")
    
    elif args.add_sample:
        # Add new sample
        results = incremental_blast.add_sample(args.add_sample)
        print(f"Added sample: {Path(args.add_sample).name}")
        print(f"- New comparisons: {results['new_comparisons']}")
        print(f"- Reused comparisons: {results['reused_comparisons']}")
    
    elif args.remove_sample:
        # Remove sample
        results = incremental_blast.remove_sample(args.remove_sample)
        print(f"Removed sample: {results['removed_sample']}")
    
    elif args.sample_dirs:
        # Compare samples
        results = incremental_blast.compare_samples_incremental(args.sample_dirs, args.force_recompute)
        print(f"Comparison Results:")
        print(f"- Total comparisons: {results['total_comparisons']}")
        print(f"- Reused: {results['reused_comparisons']}")
        print(f"- New: {results['new_comparisons']}")
        print(f"- Failed: {results['failed_comparisons']}")
    
    else:
        print("Please specify an action: --sample-dirs, --add-sample, --remove-sample, or --stats")

if __name__ == "__main__":
    main()
