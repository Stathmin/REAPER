#!/usr/bin/env python3
"""
NCBI Data Gatherer for RepOrtR

This module fetches relevant repeat sequences from NCBI based on taxonomic groups
and updates project-specific BLAST databases for annotation purposes.

Author: RepOrtR Team
Date: 2025
"""

import os
import sys
import argparse
import json
import logging
import time
import requests  # type: ignore[import-not-found]
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from datetime import datetime
import pandas as pd  # type: ignore[import-not-found]
from Bio import Entrez, SeqIO  # type: ignore[import-not-found]
from Bio.Seq import Seq  # type: ignore[import-not-found]
from Bio.SeqUtils import gc_fraction  # type: ignore[import-not-found]
import tempfile
import subprocess

# Define GC function for compatibility
def GC(sequence):
    """Calculate GC content as percentage"""
    return gc_fraction(sequence) * 100

@dataclass
class NCBISequence:
    """Represents a sequence from NCBI"""
    accession: str
    title: str
    sequence: str
    length: int
    gc_content: float
    organism: str
    taxonomy: str
    features_count: int
    description: str
    date_created: str
    
    def __post_init__(self):
        if not hasattr(self, 'length'):
            self.length = len(self.sequence)
        if not hasattr(self, 'gc_content'):
            self.gc_content = GC(self.sequence)
        if not hasattr(self, 'features_count'):
            self.features_count = 0
        if not hasattr(self, 'date_created'):
            self.date_created = datetime.now().isoformat()

class NCBIDataGatherer:
    """
    NCBI Data Gatherer Module
    
    This class fetches relevant repeat sequences from NCBI based on taxonomic groups
    and updates project-specific BLAST databases for annotation purposes.
    """
    
    def __init__(self, project_id: str = None, email: str = None):
        """
        Initialize the NCBI data gatherer
        
        Args:
            project_id: Project identifier for data organization
            email: Email for NCBI Entrez (required for high-throughput access)
        """
        self.project_id = project_id or "default"
        
        global_cfg_path = os.environ.get("REPORTR_GLOBAL_CONFIG", "projects/global_config.yaml")

        # Load email from global config if not provided
        if email is None:
            try:
                import yaml  # type: ignore[import-not-found]
                with open(global_cfg_path, 'r') as f:
                    config = yaml.safe_load(f)
                self.email = config.get('global', {}).get('ncbi_gathering', {}).get('email')
                if not self.email:
                    raise ValueError(
                        f"Email not found in global config. Please set ncbi_gathering.email in {global_cfg_path}"
                    )
            except Exception as e:
                raise ValueError(
                    f"Failed to load email from config: {e}. Please provide email parameter or set ncbi_gathering.email in {global_cfg_path}"
                )
        else:
            self.email = email
        
        # Set up Entrez
        Entrez.email = self.email
        
        # Search parameters (config-driven; no hardcoded thresholds).
        # These values must come from the global config.
        try:
            import yaml  # type: ignore[import-not-found]
            with open(global_cfg_path, "r") as f:
                config = yaml.safe_load(f) or {}
            ncbi_cfg = (config.get("global") or {}).get("ncbi_gathering") or {}
        except Exception as e:
            raise ValueError(f"Failed to load ncbi_gathering config: {e}")

        # Required config keys.
        required_keys = [
            "max_sequence_length",
            "min_sequence_length",
            "max_features",
            "batch_size",
            "repeat_terms",
        ]
        for k in required_keys:
            if k not in ncbi_cfg:
                raise KeyError(f"Missing required global.ncbi_gathering.{k} in {global_cfg_path}")

        self.max_sequence_length = int(ncbi_cfg["max_sequence_length"])
        self.min_sequence_length = int(ncbi_cfg["min_sequence_length"])
        self.max_features = int(ncbi_cfg["max_features"])
        self.batch_size = int(ncbi_cfg["batch_size"])

        self.repeat_terms = list(ncbi_cfg.get("repeat_terms") or [])
        if not self.repeat_terms:
            raise ValueError("global.ncbi_gathering.repeat_terms is empty")
        
        # Setup logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)
        
        # Create project directories using path utilities
        import sys
        sys.path.append('.')
        from post_tarean.path_utils import get_project_path, ensure_output_dir
        
        self.project_dir = get_project_path(self.project_id, "ncbi_gathering")
        self.blast_dir = self.project_dir / "blast_db"
        self.sequences_dir = self.project_dir / "sequences"
        self.metadata_dir = self.project_dir / "metadata"
        
        for directory in [self.blast_dir, self.sequences_dir, self.metadata_dir]:
            ensure_output_dir(directory)
    
    def build_search_query(self, taxid: str, additional_terms: List[str] = None) -> str:
        """
        Build NCBI search query for a taxonomic group
        
        Args:
            taxid: NCBI taxonomy ID
            additional_terms: Additional search terms to include
            
        Returns:
            Formatted search query string
        """
        # Base query with taxonomy constraint.
        #
        # Use NCBI's standard taxid syntax. The previous `porgn:__txid...` form
        # can trigger opaque "Search Backend failed" errors on E-utilities.
        query_parts = [f"txid{taxid}[Organism:exp]"]
        
        # Add repeat-related terms (config-driven base + optional extra terms).
        repeat_terms = list(self.repeat_terms)
        if additional_terms:
            repeat_terms.extend(additional_terms)

        # Deduplicate while preserving order.
        deduped: List[str] = []
        seen = set()
        for t in repeat_terms:
            if t in seen:
                continue
            seen.add(t)
            deduped.append(t)
        repeat_terms = deduped
        
        # Build title search terms
        title_terms = []
        for term in repeat_terms:
            title_terms.append(f'"{term}"[Title]')
        
        # Combine with OR operator
        title_query = " OR ".join(title_terms)
        query_parts.append(f"({title_query})")
        
        # Combine all parts with AND (no nucleotide filter for now)
        final_query = " AND ".join(query_parts)
        
        self.logger.info(f"Built query: {final_query}")
        return final_query
    
    def search_ncbi(self, query: str) -> List[str]:
        """
        Search NCBI for sequences matching the query
        
        Args:
            query: NCBI search query
            
        Returns:
            List of sequence accession numbers
        """
        self.logger.info(f"Searching NCBI with query: {query}")
        
        # NCBI E-utilities intermittently returns invalid JSON (raw tabs/newlines inside
        # string fields) on backend errors. Use XML output for robust parsing.
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

        param_variants = [
            {"tool": "RepOrtR", "email": self.email},
            {},
        ]

        import xml.etree.ElementTree as ET

        # Use history server so we can page through all IDs (no hard cap).
        for attempt in range(1, 8):
            for extra in param_variants:
                params = {
                    "db": "nucleotide",
                    "term": query,
                    "usehistory": "y",
                    "retmax": "0",
                    "retmode": "xml",
                    **extra,
                }
                try:
                    resp = requests.get(base_url, params=params, timeout=30)
                    resp.raise_for_status()
                    txt = resp.text

                    # Backend failures often come back as XML with <ERROR>...</ERROR>,
                    # but sometimes it's not XML at all. Handle both.
                    if "Search Backend failed" in txt or "<ERROR>" in txt:
                        # Best-effort extraction of the error text.
                        err = None
                        try:
                            root = ET.fromstring(txt)
                            err_el = root.find(".//ERROR")
                            if err_el is not None and (err_el.text or "").strip():
                                err = (err_el.text or "").strip()
                        except Exception:
                            err = None

                        self.logger.error(
                            f"NCBI search backend error (attempt {attempt}/7): {err or 'Search Backend failed'}"
                        )
                        time.sleep(1.5 * attempt)
                        continue

                    root = ET.fromstring(txt)
                    count_el = root.find(".//Count")
                    count = int((count_el.text or "0").strip()) if count_el is not None else 0
                    qk_el = root.find(".//QueryKey")
                    we_el = root.find(".//WebEnv")
                    if qk_el is None or we_el is None or not (qk_el.text and we_el.text):
                        raise RuntimeError("NCBI esearch did not return QueryKey/WebEnv (usehistory=y failed)")

                    query_key = (qk_el.text or "").strip()
                    webenv = (we_el.text or "").strip()
                    self.logger.info(f"Found {count} sequences (usehistory), paging IDs in batches of {self.batch_size}")

                    ids: List[str] = []
                    retstart = 0
                    while retstart < count:
                        page_params = {
                            "db": "nucleotide",
                            "query_key": query_key,
                            "WebEnv": webenv,
                            "retstart": str(retstart),
                            "retmax": str(int(self.batch_size)),
                            "retmode": "xml",
                            **extra,
                        }
                        page_txt = None
                        for page_attempt in range(1, 6):
                            try:
                                resp2 = requests.get(base_url, params=page_params, timeout=30)
                                resp2.raise_for_status()
                                page_txt = resp2.text
                                if "Search Backend failed" in page_txt or "<ERROR>" in page_txt:
                                    time.sleep(1.5 * (attempt + page_attempt))
                                    continue
                                break
                            except Exception:
                                time.sleep(1.5 * (attempt + page_attempt))
                                continue

                        if not page_txt:
                            break

                        root2 = ET.fromstring(page_txt)
                        page_ids = [el.text for el in root2.findall(".//IdList/Id") if el.text]
                        if not page_ids:
                            break
                        ids.extend(page_ids)
                        retstart += len(page_ids)
                        # Be polite to NCBI
                        time.sleep(0.34)

                    self.logger.info(f"Retrieved {len(ids)} IDs (of Count={count})")
                    return ids
                except Exception as e:
                    self.logger.error(f"NCBI search failed (attempt {attempt}/7): {e!r}")
                    time.sleep(1.5 * attempt)

        self.logger.error(
            "NCBI search failed after retries. This is often due to transient NCBI backend issues, "
            "throttling, or network/DNS problems."
        )
        return []
    
    def fetch_sequence_details(self, id_list: List[str]) -> List[NCBISequence]:
        """
        Fetch detailed sequence information from NCBI
        
        Args:
            id_list: List of NCBI sequence IDs
            
        Returns:
            List of NCBISequence objects
        """
        self.logger.info(f"Fetching details for {len(id_list)} sequences")
        
        sequences = []
        
        # Process in batches with retries (NCBI backend can be transiently flaky)
        for i in range(0, len(id_list), self.batch_size):
            batch_ids = id_list[i:i + self.batch_size]
            
            try:
                # Fetch sequence records
                records = None
                for attempt in range(1, 8):
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="text")
                        records = list(SeqIO.parse(handle, "genbank"))
                        handle.close()
                        break
                    except Exception as e:
                        self.logger.error(
                            f"Entrez.efetch failed (batch {i//self.batch_size + 1} attempt {attempt}/7): {e!r}"
                        )
                        time.sleep(1.5 * attempt)

                if records is None:
                    continue
                
                for record in records:
                    # Skip sequences that are too long or over-annotated
                    if len(record.seq) > self.max_sequence_length:
                        self.logger.debug(f"Skipping {record.id}: too long ({len(record.seq)} bp)")
                        continue
                    
                    if len(record.seq) < self.min_sequence_length:
                        self.logger.debug(f"Skipping {record.id}: too short ({len(record.seq)} bp)")
                        continue
                    
                    # Count features (annotations)
                    features_count = len(record.features)
                    if features_count > self.max_features:
                        self.logger.debug(f"Skipping {record.id}: too many features ({features_count})")
                        continue
                    
                    # Extract taxonomy
                    taxonomy = "Unknown"
                    if hasattr(record, 'annotations') and 'taxonomy' in record.annotations:
                        taxonomy = "; ".join(record.annotations['taxonomy'])
                    
                    # Extract organism
                    organism = "Unknown"
                    if hasattr(record, 'annotations') and 'organism' in record.annotations:
                        organism = record.annotations['organism']
                    
                    # Create sequence object
                    seq_obj = NCBISequence(
                        accession=record.id,
                        title=record.name,
                        sequence=str(record.seq),
                        length=len(record.seq),
                        gc_content=GC(str(record.seq)),
                        organism=organism,
                        taxonomy=taxonomy,
                        features_count=features_count,
                        description=record.description,
                        date_created=datetime.now().isoformat()
                    )
                    
                    sequences.append(seq_obj)
                
                self.logger.info(f"Processed batch {i//self.batch_size + 1}, got {len(sequences)} valid sequences")
                
                # Rate limiting (keep under NCBI guidelines)
                time.sleep(0.34)
                
            except Exception as e:
                self.logger.error(f"Failed to fetch batch {i//self.batch_size + 1}: {e}")
                time.sleep(1.0)
                continue
        
        self.logger.info(f"Total valid sequences: {len(sequences)}")
        return sequences
    
    def filter_sequences(self, sequences: List[NCBISequence]) -> List[NCBISequence]:
        """
        Filter sequences based on quality criteria
        
        Args:
            sequences: List of NCBISequence objects
            
        Returns:
            Filtered list of sequences
        """
        self.logger.info(f"Filtering {len(sequences)} sequences")
        
        filtered = []
        
        for seq in sequences:
            # Check GC content (reasonable range)
            if seq.gc_content < 20 or seq.gc_content > 80:
                self.logger.info(f"Skipping {seq.accession}: extreme GC content ({seq.gc_content:.1f}%)")
                continue
            
            self.logger.info(f"Keeping {seq.accession}: length={seq.length}, gc={seq.gc_content:.1f}%, features={seq.features_count}")
            filtered.append(seq)
        
        self.logger.info(f"Filtered to {len(filtered)} high-quality sequences")
        return filtered
    
    def save_sequences(self, sequences: List[NCBISequence], output_file: str = None) -> str:
        """
        Save sequences to FASTA file
        
        Args:
            sequences: List of NCBISequence objects
            output_file: Output file path (optional)
            
        Returns:
            Path to saved FASTA file
        """
        if output_file is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = self.sequences_dir / f"ncbi_repeats_{timestamp}.fasta"
        
        with open(output_file, 'w') as f:
            for seq in sequences:
                # Create a BLAST-friendly header.
                # Use a standard seqid prefix so `makeblastdb -parse_seqids` enables `blastdbcmd -entry <accession>`.
                # Keep the descriptive fields after whitespace so they remain in the title but don't break ID parsing.
                header = (
                    f">ref|{seq.accession}| "
                    f"organism={seq.organism} title={seq.title} len={seq.length} gc={seq.gc_content:.1f}%"
                )
                f.write(f"{header}\n{seq.sequence}\n")
        
        self.logger.info(f"Saved {len(sequences)} sequences to {output_file}")
        return str(output_file)
    
    def create_blast_database(self, fasta_file: str, db_name: str = None) -> str:
        """
        Create BLAST database from FASTA file
        
        Args:
            fasta_file: Path to FASTA file
            db_name: Database name (optional)
            
        Returns:
            Path to BLAST database
        """
        if db_name is None:
            db_name = Path(fasta_file).stem
        
        db_path = self.blast_dir / db_name
        
        try:
            # Create BLAST database
            cmd = [
                "makeblastdb",
                "-in", fasta_file,
                "-dbtype", "nucl",
                "-parse_seqids",
                "-out", str(db_path),
                "-title", f"{db_name}_repeat_database"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info(f"Created BLAST database: {db_path}")
            
            return str(db_path)
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to create BLAST database: {e}")
            self.logger.error(f"Error output: {e.stderr}")
            return ""
    
    def save_metadata(self, sequences: List[NCBISequence], taxid: str) -> str:
        """
        Save sequence metadata to CSV file
        
        Args:
            sequences: List of NCBISequence objects
            taxid: Taxonomy ID used for search
            
        Returns:
            Path to metadata file
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        metadata_file = self.metadata_dir / f"ncbi_metadata_{taxid}_{timestamp}.csv"
        
        # Convert to DataFrame
        data = []
        for seq in sequences:
            data.append({
                'accession': seq.accession,
                'title': seq.title,
                'length': seq.length,
                'gc_content': seq.gc_content,
                'organism': seq.organism,
                'taxonomy': seq.taxonomy,
                'features_count': seq.features_count,
                'description': seq.description,
                'date_created': seq.date_created
            })
        
        df = pd.DataFrame(data)
        df.to_csv(metadata_file, index=False)
        
        self.logger.info(f"Saved metadata to {metadata_file}")
        return str(metadata_file)
    
    def gather_sequences_for_taxonomy(self, taxid: str, additional_terms: List[str] = None) -> Dict:
        """
        Main method to gather sequences for a taxonomic group
        
        Args:
            taxid: NCBI taxonomy ID
            additional_terms: Additional search terms
            
        Returns:
            Dictionary with results summary
        """
        self.logger.info(f"Gathering sequences for taxonomy ID: {taxid}")
        
        # Build search query
        query = self.build_search_query(taxid, additional_terms)
        
        # Search NCBI
        id_list = self.search_ncbi(query)
        if not id_list:
            return {
                "error": (
                    "No sequences found (or NCBI search failed). "
                    "See the gather_ncbi_sequences log for the exact E-utilities error."
                )
            }
        
        # Fetch sequence details
        sequences = self.fetch_sequence_details(id_list)
        if not sequences:
            return {"error": "No valid sequences retrieved"}
        
        # Filter sequences
        filtered_sequences = self.filter_sequences(sequences)
        if not filtered_sequences:
            return {"error": "No sequences passed quality filters"}
        
        # Save sequences
        fasta_file = self.save_sequences(filtered_sequences)
        
        # Create BLAST database
        blast_db = self.create_blast_database(fasta_file)
        
        # Save metadata
        metadata_file = self.save_metadata(filtered_sequences, taxid)
        
        # Prepare results summary
        results = {
            "taxid": taxid,
            "query": query,
            "total_found": len(id_list),
            "valid_sequences": len(filtered_sequences),
            "fasta_file": fasta_file,
            "blast_database": blast_db,
            "metadata_file": metadata_file,
            "sequences": filtered_sequences
        }
        
        # Calculate statistics
        if filtered_sequences:
            lengths = [seq.length for seq in filtered_sequences]
            gc_contents = [seq.gc_content for seq in filtered_sequences]
            
            results["statistics"] = {
                "mean_length": sum(lengths) / len(lengths),
                "min_length": min(lengths),
                "max_length": max(lengths),
                "mean_gc": sum(gc_contents) / len(gc_contents),
                "min_gc": min(gc_contents),
                "max_gc": max(gc_contents)
            }
        
        self.logger.info(f"Gathering complete: {len(filtered_sequences)} sequences saved")
        return results

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(description="NCBI data gatherer for repeat sequences")
    parser.add_argument("taxid", help="NCBI taxonomy ID (e.g., 147389 for Triticeae)")
    parser.add_argument("--project-id", help="Project identifier")
    parser.add_argument("--email", help="Email for NCBI Entrez access")
    parser.add_argument("--additional-terms", nargs="+", help="Additional search terms")
    parser.add_argument("--output-dir", help="Output directory")
    
    args = parser.parse_args()
    
    # Initialize gatherer
    gatherer = NCBIDataGatherer(args.project_id, args.email)
    
    # Gather sequences
    results = gatherer.gather_sequences_for_taxonomy(args.taxid, args.additional_terms)
    
    if "error" in results:
        print(f"Error: {results['error']}")
        sys.exit(1)
    
    # Print summary
    print(f"\nNCBI Data Gathering Complete!")
    print(f"Taxonomy ID: {results['taxid']}")
    print(f"Query: {results['query']}")
    print(f"Total found: {results['total_found']}")
    print(f"Valid sequences: {results['valid_sequences']}")
    print(f"FASTA file: {results['fasta_file']}")
    print(f"BLAST database: {results['blast_database']}")
    print(f"Metadata file: {results['metadata_file']}")
    
    if "statistics" in results:
        stats = results["statistics"]
        print(f"\nStatistics:")
        print(f"  Length: {stats['mean_length']:.0f} bp (range: {stats['min_length']}-{stats['max_length']})")
        print(f"  GC content: {stats['mean_gc']:.1f}% (range: {stats['min_gc']:.1f}-{stats['max_gc']:.1f}%)")

if __name__ == "__main__":
    main()
