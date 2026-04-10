#!/usr/bin/env python3
"""
Real Post-TAREAN Analysis
Processes actual TAREAN outputs to generate meaningful reports
"""

import os
import sys
import glob
import pandas as pd
import sqlite3
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RealTAREANAnalyzer:
    """Analyzes real TAREAN outputs"""
    
    def __init__(self, project_id, sample_id):
        self.project_id = project_id
        self.sample_id = sample_id
        self.tarean_dir = f"projects/{project_id}/samples/{sample_id}/tarean"
        self.output_dir = f"projects/{project_id}/samples/{sample_id}/post_tarean"
        
    def analyze_clusters(self):
        """Analyze cluster information from TAREAN outputs"""
        clusters_dir = f"{self.tarean_dir}/seqclust/clustering/clusters"
        
        # Check if we're running from the post_tarean directory
        if not os.path.exists(clusters_dir):
            # Try relative path from post_tarean directory
            clusters_dir = f"../{self.tarean_dir}/seqclust/clustering/clusters"
        
        if not os.path.exists(clusters_dir):
            logger.warning(f"Clusters directory not found: {clusters_dir}")
            return pd.DataFrame()
        
        cluster_data = []
        cluster_dirs = glob.glob(f"{clusters_dir}/dir_CL*")
        
        for cluster_dir in cluster_dirs:
            cluster_num = os.path.basename(cluster_dir).replace('dir_CL', '').replace('000', '')
            
            # Get cluster size from hitsort.cls
            cls_file = f"{cluster_dir}/hitsort.cls"
            cluster_size = 0
            if os.path.exists(cls_file):
                with open(cls_file, 'r') as f:
                    cluster_size = len([line for line in f if line.startswith('>')])
            
            # Check for TAREAN report
            tarean_report = f"{cluster_dir}/tarean/report.html"
            has_tarean = os.path.exists(tarean_report)
            
            cluster_data.append({
                'cluster_id': f"CL{cluster_num}",
                'cluster_size': cluster_size,
                'has_tarean_report': has_tarean,
                'cluster_path': cluster_dir
            })
        
        return pd.DataFrame(cluster_data)
    
    def analyze_sequences(self):
        """Analyze sequence information from TAREAN databases"""
        sequences_db = f"{self.tarean_dir}/sequences.db"
        hitsort_db = f"{self.tarean_dir}/hitsort.db"
        
        seq_stats = {}
        
        # Analyze sequences database
        if not os.path.exists(sequences_db):
            # Try relative path from post_tarean directory
            sequences_db = f"../{sequences_db}"
        if os.path.exists(sequences_db):
            try:
                conn = sqlite3.connect(sequences_db)
                cursor = conn.cursor()
                
                # Get table names
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
                tables = cursor.fetchall()
                seq_stats['database_tables'] = [table[0] for table in tables]
                
                # Get sequence count if table exists
                if ('sequences',) in tables:
                    cursor.execute("SELECT COUNT(*) FROM sequences")
                    seq_count = cursor.fetchone()[0]
                    seq_stats['total_sequences'] = seq_count
                
                conn.close()
            except Exception as e:
                logger.error(f"Error analyzing sequences database: {e}")
                seq_stats['database_error'] = str(e)
        
        # Analyze hitsort database
        if not os.path.exists(hitsort_db):
            # Try relative path from post_tarean directory
            hitsort_db = f"../{hitsort_db}"
        if os.path.exists(hitsort_db):
            try:
                conn = sqlite3.connect(hitsort_db)
                cursor = conn.cursor()
                
                # Get table names
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
                tables = cursor.fetchall()
                seq_stats['hitsort_tables'] = [table[0] for table in tables]
                
                # Get cluster count if table exists
                if ('clusters',) in tables:
                    cursor.execute("SELECT COUNT(*) FROM clusters")
                    cluster_count = cursor.fetchone()[0]
                    seq_stats['total_clusters'] = cluster_count
                
                conn.close()
            except Exception as e:
                logger.error(f"Error analyzing hitsort database: {e}")
                seq_stats['hitsort_error'] = str(e)
        
        return seq_stats
    
    def analyze_contigs(self):
        """Analyze contig information"""
        contigs_file = f"{self.tarean_dir}/contigs.fasta"
        
        contig_stats = {}
        
        if not os.path.exists(contigs_file):
            # Try relative path from post_tarean directory
            contigs_file = f"../{contigs_file}"
        if os.path.exists(contigs_file):
            try:
                with open(contigs_file, 'r') as f:
                    content = f.read()
                    contigs = content.split('>')[1:]  # Skip first empty element
                    contig_stats['total_contigs'] = len(contigs)
                    
                    # Calculate total length
                    total_length = 0
                    for contig in contigs:
                        lines = contig.strip().split('\n')
                        if len(lines) > 1:
                            sequence = ''.join(lines[1:])
                            total_length += len(sequence)
                    
                    contig_stats['total_length'] = total_length
                    contig_stats['avg_length'] = total_length / len(contigs) if contigs else 0
                    
            except Exception as e:
                logger.error(f"Error analyzing contigs: {e}")
                contig_stats['error'] = str(e)
        
        return contig_stats
    
    def generate_reports(self):
        """Generate comprehensive reports from real TAREAN data"""
        logger.info(f"Analyzing TAREAN outputs for {self.sample_id}")
        
        # Analyze different components
        clusters_df = self.analyze_clusters()
        seq_stats = self.analyze_sequences()
        contig_stats = self.analyze_contigs()
        
        # Create output directory
        if not os.path.exists(self.output_dir):
            # Try relative path from post_tarean directory
            self.output_dir = f"../{self.output_dir}"
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Ensure we're using absolute paths for output
        if not os.path.isabs(self.output_dir):
            self.output_dir = os.path.abspath(self.output_dir)
        
        # Generate CSV report
        csv_report = f"{self.output_dir}/real_analysis.csv"
        print(f"Creating CSV report at: {csv_report}")
        if not clusters_df.empty:
            clusters_df.to_csv(csv_report, index=False)
            logger.info(f"Generated cluster analysis CSV: {csv_report}")
            print(f"CSV file created successfully")
        else:
            print("No cluster data to write to CSV")
        
        # Generate summary report
        summary_report = f"{self.output_dir}/real_summary.txt"
        with open(summary_report, 'w') as f:
            f.write(f"Real TAREAN Analysis Summary for {self.sample_id}\n")
            f.write("=" * 50 + "\n")
            f.write(f"Project: {self.project_id}\n")
            f.write(f"Sample: {self.sample_id}\n")
            f.write(f"Analysis Date: {pd.Timestamp.now()}\n\n")
            
            f.write("CLUSTER ANALYSIS:\n")
            f.write("-" * 20 + "\n")
            if not clusters_df.empty:
                f.write(f"Total clusters found: {len(clusters_df)}\n")
                f.write(f"Clusters with TAREAN reports: {clusters_df['has_tarean_report'].sum()}\n")
                f.write(f"Total sequences in clusters: {clusters_df['cluster_size'].sum()}\n")
                f.write(f"Average cluster size: {clusters_df['cluster_size'].mean():.1f}\n")
            else:
                f.write("No cluster data available\n")
            
            f.write("\nSEQUENCE ANALYSIS:\n")
            f.write("-" * 20 + "\n")
            for key, value in seq_stats.items():
                f.write(f"{key}: {value}\n")
            
            f.write("\nCONTIG ANALYSIS:\n")
            f.write("-" * 20 + "\n")
            for key, value in contig_stats.items():
                f.write(f"{key}: {value}\n")
        
        logger.info(f"Generated summary report: {summary_report}")
        
        # Generate Excel-like report (CSV with multiple sheets)
        excel_report = f"{self.output_dir}/real_analysis.xlsx"
        with pd.ExcelWriter(excel_report, engine='openpyxl') as writer:
            if not clusters_df.empty:
                clusters_df.to_excel(writer, sheet_name='Clusters', index=False)
            
            # Create summary sheet
            summary_data = []
            for key, value in seq_stats.items():
                summary_data.append(['Sequence', key, str(value)])
            for key, value in contig_stats.items():
                summary_data.append(['Contig', key, str(value)])
            
            summary_df = pd.DataFrame(summary_data, columns=['Category', 'Metric', 'Value'])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
        
        logger.info(f"Generated Excel report: {excel_report}")
        
        # Mark analysis as complete
        complete_file = f"{self.output_dir}/analysis_complete.txt"
        with open(complete_file, 'w') as f:
            f.write(f"Real TAREAN analysis completed for {self.sample_id}\n")
            f.write(f"Timestamp: {pd.Timestamp.now()}\n")
        
        logger.info(f"Analysis complete: {complete_file}")
        
        return {
            'csv_report': csv_report,
            'summary_report': summary_report,
            'excel_report': excel_report,
            'complete_file': complete_file
        }

def main():
    if len(sys.argv) != 3:
        print("Usage: python real_analysis.py <project_id> <sample_id>")
        sys.exit(1)
    
    project_id = sys.argv[1]
    sample_id = sys.argv[2]
    
    analyzer = RealTAREANAnalyzer(project_id, sample_id)
    results = analyzer.generate_reports()
    
    print(f"Real TAREAN analysis completed for {sample_id}")
    print(f"Reports generated: {list(results.keys())}")

if __name__ == "__main__":
    main() 