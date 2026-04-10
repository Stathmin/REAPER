#!/usr/bin/env python3
"""
Consolidated BLAST Analysis Module
Combines functionality from blast.py, blast_safe.py, blast_Kroupin.py, and blastold.py
"""

import os
import re
import glob
import tempfile
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import logging
from typing import List, Dict, Optional, Tuple

# Optional imports
try:
    from rapidfuzz import fuzz
    FUZZ_AVAILABLE = True
except ImportError:
    FUZZ_AVAILABLE = False

logger = logging.getLogger(__name__)


class BLASTAnalyzer:
    """Consolidated BLAST analysis class"""
    
    def __init__(self, config=None, config_manager=None, project_id: Optional[str] = None):
        """
        Initialise BLAST analyzer.

        Parameters
        ----------
        config
            Either a PipelineConfig instance or a project_id string (legacy).
        config_manager
            Optional ConfigManager instance providing `.config` and `.project_id`.
        project_id
            Optional explicit project identifier; mainly used when resolving
            TAREAN paths for post-TAREAN integration.
        """
        self.config_manager = None
        self.config = None
        self.project_id = project_id

        # Prefer an explicit ConfigManager instance when provided – this is
        # the path used by the post-TAREAN pipeline integration.
        if config_manager is not None:
            self.config_manager = config_manager
            self.config = config_manager.config
            if self.project_id is None:
                self.project_id = getattr(config_manager, "project_id", None)
        # Legacy-style initialisation by project_id string.
        elif isinstance(config, str):
            from config import create_config_for_project
            self.config_manager = create_config_for_project(config)
            self.config = self.config_manager.config
            if self.project_id is None:
                # Fall back to the project_id argument or the string itself.
                self.project_id = getattr(self.config_manager, "project_id", config)
        # Direct config object (PipelineConfig) – typically used in tests or
        # lower-level callers that manage ConfigManager themselves.
        elif config is not None:
            self.config = config

        # db_dict maps logical DB names (e.g. "local") to BLAST database
        # prefixes suitable for -db; fasta_dict maps the same keys to FASTA
        # paths used for sequence slicing (legacy self-BLAST mode).
        self.db_dict, self.fasta_dict = self._discover_databases()
        
    def _discover_databases(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        """Discover available BLAST databases and their FASTA sources."""
        db_dict: Dict[str, str] = {}
        fasta_dict: Dict[str, str] = {}
        
        # Use config database paths if available
        if self.config:
            db_paths = {
                'ncbi': self.config.databases.ncbi_db_path,
                'local': self.config.databases.local_db_path,
                'comp': self.config.databases.comparative_db_path,
                'ref': self.config.databases.reference_db_path,
            }
            
            # NCBI database (if present)
            ncbi_pattern = os.path.join(db_paths['ncbi'], '*.nhr')
            if glob.glob(ncbi_pattern):
                db_prefix = os.path.join(db_paths['ncbi'], self.config.databases.database_names['ncbi'])
                db_dict['ncbi'] = db_prefix
                db_dict['ncbi_x3'] = db_prefix  # same BLAST DB, different query treatment
                fasta_path = db_prefix if db_prefix.endswith('.fasta') else f"{db_prefix}.fasta"
                fasta_dict['ncbi'] = fasta_path
                fasta_dict['ncbi_x3'] = fasta_path
            
            # Local database
            local_pattern = os.path.join(db_paths['local'], '*.nhr')
            if glob.glob(local_pattern):
                db_prefix = os.path.join(db_paths['local'], self.config.databases.database_names['local'])
                db_dict['local'] = db_prefix
                db_dict['local_x3'] = db_prefix
                fasta_path = db_prefix if db_prefix.endswith('.fasta') else f"{db_prefix}.fasta"
                fasta_dict['local'] = fasta_path
                fasta_dict['local_x3'] = fasta_path
            
            # Comparative database
            comp_pattern = os.path.join(db_paths['comp'], '*.nhr')
            if glob.glob(comp_pattern):
                db_prefix = os.path.join(db_paths['comp'], self.config.databases.database_names['comp'])
                db_dict['comp'] = db_prefix
                db_dict['comp_x3'] = db_prefix
                fasta_path = db_prefix if db_prefix.endswith('.fasta') else f"{db_prefix}.fasta"
                fasta_dict['comp'] = fasta_path
                fasta_dict['comp_x3'] = fasta_path
            
            # Reference database
            ref_pattern = os.path.join(db_paths['ref'], '*.nhr')
            if glob.glob(ref_pattern):
                db_prefix = os.path.join(db_paths['ref'], self.config.databases.database_names['ref'])
                db_dict['ref'] = db_prefix
                fasta_path = db_prefix if db_prefix.endswith('.fasta') else f"{db_prefix}.fasta"
                fasta_dict['ref'] = fasta_path
        else:
            # Fallback to default paths (used mainly in standalone mode)
            if glob.glob('./ncbi_repeats_db/*.nhr'):
                db_prefix = './ncbi_repeats_db/ncbi_repeats'
                db_dict['ncbi'] = db_prefix
                db_dict['ncbi_x3'] = db_prefix
                fasta_path = f"{db_prefix}.fasta"
                fasta_dict['ncbi'] = fasta_path
                fasta_dict['ncbi_x3'] = fasta_path
            
            if glob.glob('./local_db_solo/*.nhr'):
                db_prefix = './local_db_solo/multifasta'
                db_dict['local'] = db_prefix
                db_dict['local_x3'] = db_prefix
                fasta_path = f"{db_prefix}.fasta"
                fasta_dict['local'] = fasta_path
                fasta_dict['local_x3'] = fasta_path
            
            if glob.glob('./comparatives_db/*.nhr'):
                db_prefix = './comparatives_db/COMPBASE'
                db_dict['comp'] = db_prefix
                db_dict['comp_x3'] = db_prefix
                fasta_path = f"{db_prefix}.fasta"
                fasta_dict['comp'] = fasta_path
                fasta_dict['comp_x3'] = fasta_path
            
            if glob.glob('./important_db/*.nhr'):
                db_prefix = './important_db/reference'
                db_dict['ref'] = db_prefix
                fasta_dict['ref'] = f"{db_prefix}.fasta"
        
        logger.info(f"Discovered databases: {list(db_dict.keys())}")
        return db_dict, fasta_dict
    
    def _get_last_name(self, filename: str) -> str:
        """Extract filename without extension from path"""
        return re.sub(r'\.\w+', '', filename.split('/')[-1])
    
    def _create_x3_variant(self, filename: str, output_dir: str = None) -> str:
        """Create a version of given fasta with tripled sequences"""
        if output_dir is None:
            output_dir = os.getcwd()
        
        with open(filename, 'r') as file:
            content = file.read()
            # Split by '>' and process each sequence
            sequences = content.split('>')
            processed_sequences = []
            
            for seq in filter(None, sequences):
                lines = seq.split('\n')
                header = lines[0]
                sequence = ''.join(lines[1:])
                # Triple the sequence
                tripled_seq = sequence * 3
                processed_sequences.append(f'>{header}\n{tripled_seq}')
            
            x3_content = '\n'.join(processed_sequences)
            x3_path = os.path.join(output_dir, f'x3_{self._get_last_name(filename)}.fasta')
            
            with open(x3_path, 'w') as outfile:
                outfile.write(x3_content)
            
            return x3_path
    
    def _run_blast_command(self, input_fasta: str, db: str, task: str = 'blastn',
                          blastn_path: str = None, num_threads: int = None,
                          verbose: bool = None, tmpdirname: str = None) -> str:
        """Run BLAST command and return output file path"""
        # Use config values if available and not provided
        if self.config:
            if blastn_path is None:
                blastn_path = self.config.blast.blastn_path
            if num_threads is None:
                num_threads = self.config.blast.default_threads
            if verbose is None:
                verbose = self.config.blast.verbose
        else:
            # Fallback to defaults
            if blastn_path is None:
                blastn_path = 'blastn'
            if num_threads is None:
                num_threads = 20
            if verbose is None:
                verbose = False
        
        with tempfile.NamedTemporaryFile(dir=tmpdirname, delete=False) as out_file, \
             tempfile.NamedTemporaryFile(dir=tmpdirname, delete=False) as err_file:
            
            # Build command
            word_size_param = "-word_size 6" if task == "blastn" else ""
            command = (f'{blastn_path} -query "{input_fasta}" -db {db} -task {task} '
                      f'{word_size_param} -num_threads {num_threads} -outfmt 6')
            
            if verbose:
                logger.info(f"Running command: {command}")
            
            # Execute command
            result = subprocess.run(
                command,
                shell=True,
                stdout=out_file,
                stderr=err_file,
                text=True
            )
        
        # Check for errors
        if result.returncode != 0:
            with open(err_file.name) as f:
                error_msg = f.read()
            raise RuntimeError(f"BLAST error: {error_msg}")
        
        return out_file.name
    
    def _parse_fasta(self, fasta_path: str) -> Dict[str, str]:
        """Parse FASTA file into dictionary"""
        parsed_db = {}
        current_key = None
        current_seq = []
        
        with open(fasta_path, 'r') as file:
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_key:
                        parsed_db[current_key] = ''.join(current_seq)
                    current_key = line[1:]
                    current_seq = []
                else:
                    current_seq.append(line)
        
        if current_key:
            parsed_db[current_key] = ''.join(current_seq)
        
        return parsed_db
    
    def _get_from_db(self, db: str, item: str) -> str:
        """Get specific item from FASTA database"""
        parsed_db = self._parse_fasta(db)
        return parsed_db.get(item, '')
    
    def _get_from_db_multiple(self, db: str, pattern: str) -> str:
        """Get multiple items from FASTA database matching pattern"""
        parsed_db = self._parse_fasta(db)
        patterned_keys = [x for x in parsed_db.keys() if re.match(pattern, x)]
        
        if not patterned_keys:
            logger.warning(f"No matches found for pattern {pattern} in {db}")
            return ""
        
        # Return full FASTA format with headers
        fasta_entries = []
        for key in patterned_keys:
            fasta_entries.append(f">{key}\n{parsed_db[key]}")
        
        return "\n".join(fasta_entries)
    
    def _protect_brackets(self, string: str) -> str:
        """Escape brackets for regex matching"""
        return string.replace('(', r'\(').replace(')', r'\)')
    
    def _slice_databases(self, patterns_list: List[str], dbs_list: List[str], 
                        output_path: str) -> str:
        """Create FASTA file with sequences from given databases"""
        logger.info(f"Slicing databases: patterns={patterns_list}, dbs={dbs_list}, output={output_path}")
        content_parts = []
        
        for pattern, db_name in zip(patterns_list, dbs_list):
            if db_name not in self.db_dict:
                logger.warning(f"Database {db_name} not found")
                continue
            
            fasta_path = self.fasta_dict.get(db_name)
            if not fasta_path:
                logger.warning(f"No FASTA source registered for database {db_name}")
                continue
            
            logger.info(f"Extracting from {fasta_path} with pattern {pattern}")
            content = self._get_from_db_multiple(fasta_path, pattern)
            logger.info(f"Extracted content length: {len(content)}")
            if content:
                content_parts.append(content)
        
        full_content = "\n".join(content_parts)
        logger.info(f"Total content length: {len(full_content)}")
        
        with open(output_path, 'w') as file:
            file.write(full_content)
        
        logger.info(f"Created sequence file: {output_path}")
        return output_path
    
    def _blast_coordinator(self, tasks: List[str],
                          input_fasta_abs: str,
                          num_threads: int,
                          verbose: bool,
                          fmt_header: List[str],
                          dbs: List[str]) -> pd.DataFrame:
        """
        Coordinate BLAST runs across multiple databases and tasks.

        Parameters
        ----------
        tasks
            BLAST tasks to run (e.g. ["megablast", "dc-megablast"]).
        input_fasta_abs
            Absolute path to the query FASTA file.
        num_threads
            Number of threads for BLAST.
        verbose
            Whether to log full BLAST command lines.
        fmt_header
            Column names for the `-outfmt 6` table.
        dbs
            Logical database names to query (subset of `self.db_dict.keys()`).
        """
        dfs = []
        
        for db_name in dbs:
            if db_name not in self.db_dict:
                logger.warning(f"Requested database {db_name} is not available; skipping.")
                continue

            for task in tasks:
                with tempfile.TemporaryDirectory() as tmpdirname:
                    # Handle x3 variants
                    if '_x3' in db_name:
                        current_input_fasta = self._create_x3_variant(input_fasta_abs, tmpdirname)
                    else:
                        current_input_fasta = input_fasta_abs
                    
                    # Run BLAST
                    result_file = self._run_blast_command(
                        input_fasta=current_input_fasta,
                        db=self.db_dict[db_name],
                        task=task,
                        num_threads=num_threads,
                        verbose=verbose,
                        tmpdirname=tmpdirname
                    )
                    
                    # Parse results
                    df = pd.read_csv(result_file, names=fmt_header, delimiter='\t')
                    df['task'] = task
                    df['db'] = db_name
                    dfs.append(df)
        
        if not dfs:
            logger.warning(
                "BLAST coordinator received no usable databases or produced no "
                "results; returning empty DataFrame."
            )
            return pd.DataFrame(columns=fmt_header + ['task', 'db'])

        return pd.concat(dfs, ignore_index=True)

    def _build_query_fasta_for_subjects(
        self,
        subjects: List[str],
        tmpdirname: Optional[str] = None,
    ) -> Optional[str]:
        """
        Build a query FASTA from TAREAN/post-TAREAN consensus sequences.

        This is the new code path used by the post-TAREAN pipeline when
        querying external databases (e.g. NCBI-derived repeat DBs). It
        treats `subjects` as sample identifiers (e.g. KA12) rather than
        trying to slice database FASTA files by KA* patterns.
        """
        if tmpdirname is None:
            tmpdirname = os.getcwd()

        if self.config is None or self.config_manager is None:
            logger.warning(
                "ConfigManager/config missing on BLASTAnalyzer; cannot resolve "
                "TAREAN outputs for query construction."
            )
            return None

        try:
            from post_tarean.io_helpers import resolve_tarean_path
            from post_tarean.utils import get_tareans
        except ImportError as exc:
            logger.error(f"Failed to import post_tarean helpers for query construction: {exc}")
            return None

        records: List[str] = []

        for subject in subjects:
            tarean_path, candidates = resolve_tarean_path(self.config_manager, self.config, subject)
            if not str(tarean_path):
                logger.warning(
                    f"Could not find TAREAN data for subject {subject}. "
                    f"Tried: {', '.join(str(c) for c in candidates)}"
                )
                continue

            try:
                df = get_tareans(tarean_path, subject)
            except Exception as exc:  # pragma: no cover - defensive
                logger.error(f"Failed to load TAREAN consensus data for {subject}: {exc}")
                continue

            # The `cluster` column already carries headers like
            #   KA12_CL1_TR_1_x_568nt
            # which encode query length information used later by
            # `_clean_dataframe` via the `extract_length` helper.
            for _, row in df.iterrows():
                header = str(row.get('cluster', '')).strip()
                seq = str(row.get('seq', '')).replace('\n', '')
                if not header or not seq:
                    continue
                records.append(f">{header}\n{seq}")

        if not records:
            logger.warning("No query sequences could be constructed from TAREAN outputs.")
            return None

        safe_subjects = [re.sub(r'[^A-Za-z0-9_.-]+', '_', s) for s in subjects]
        output_path = os.path.join(tmpdirname, f'post_tarean_queries_{"_".join(safe_subjects)}.fasta')

        with open(output_path, 'w') as out_f:
            out_f.write("\n".join(records) + "\n")

        logger.info(
            f"Constructed query FASTA for subjects {subjects}: "
            f"{output_path} ({len(records)} sequences)"
        )
        return output_path
    
    def _clean_dataframe(self, df: pd.DataFrame, tasks: List[str]) -> pd.DataFrame:
        """Clean and filter BLAST results dataframe"""
        df = df.copy()
        
        # Get thresholds from config
        evalue_threshold = 0.1  # Default
        pident_threshold = 60.0  # Default
        if self.config:
            evalue_threshold = self.config.blast.evalue_threshold
            pident_threshold = self.config.blast.pident_threshold
        
        # Filter by E-value and percent identity
        df = df[df['evalue'] <= evalue_threshold].reset_index(drop=True)
        df = df[df['pident'] >= pident_threshold].reset_index(drop=True)
        
        # Set categorical types
        df['task'] = pd.Categorical(df['task'], tasks)
        df['db'] = pd.Categorical(df['db'], list(self.db_dict.keys()))
        
        # Extract query length from sequence names
        def extract_length(seq_name):
            numbers = re.findall(r'\d+', seq_name)
            if numbers:
                return int(numbers[-1])
            else:
                # Fallback: try to extract length from the sequence name format
                # Example: KP5_CL1_TR_1_x_568nt -> 568
                length_match = re.search(r'_x_(\d+)nt', seq_name)
                if length_match:
                    return int(length_match.group(1))
                else:
                    # Default length if no pattern matches
                    return 100
        
        df['qlength'] = df['qseqid'].apply(extract_length)
        
        # Adjust lengths for x3 variants
        x3_mask = df['db'].apply(lambda x: '_x3' in x)
        df.loc[x3_mask, 'qlength'] = df.loc[x3_mask, 'qlength'] * 3
        
        # Convert length to integer
        df['length'] = df['length'].astype(int)
        
        return df
    
    def _merge_intervals(self, intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
        """Merge overlapping intervals"""
        if not intervals:
            return []
        
        intervals = sorted(intervals, key=lambda x: x[0])
        merged = [intervals[0]]
        
        for current in intervals[1:]:
            prev = merged[-1]
            if current[0] <= prev[1] + 1:
                merged[-1] = (prev[0], max(prev[1], current[1]))
            else:
                merged.append(current)
        
        return merged
    
    def _sum_lengths(self, df: pd.DataFrame) -> int:
        """Calculate cumulative coverage of query by all HSPs with subject"""
        if df.empty:
            return 0
        
        intervals = list(zip(df['qstart'], df['qend']))
        merged = self._merge_intervals(intervals)
        return sum(end - start + 1 for start, end in merged)
    
    def _return_best_examples(self, df: pd.DataFrame, evalue_max_olig: float = 0.001) -> pd.DataFrame:
        """Return best examples from BLAST results"""
        if df.empty:
            return df
        
        # Calculate coverage
        total_coverage = self._sum_lengths(df)
        qlength = df['qlength'].iloc[0]
        coverage_ratio = total_coverage / qlength
        
        # Determine coverage type
        if coverage_ratio >= 0.95:
            coverage_type = 'near-full'
        elif coverage_ratio >= 0.8:
            coverage_type = 'composite'
        elif coverage_ratio >= 0.6:
            coverage_type = 'partial'
        else:
            coverage_type = 'weak'
        
        # Add coverage information
        df['total_coverage'] = total_coverage
        df['coverage_ratio'] = coverage_ratio
        df['coverage_type'] = coverage_type
        
        # Filter by E-value for oligo hits
        if coverage_type in ['weak', 'partial']:
            df = df[df['evalue'] <= evalue_max_olig]
        
        # Return longest interval if multiple hits
        if len(df) > 1:
            df = df.loc[df['length'].idxmax():df['length'].idxmax()]
        
        return df
    
    def _sort_for_rules(self, df: pd.DataFrame, subjects: List[str]) -> pd.DataFrame:
        """Sort results according to priority rules"""
        if df.empty:
            return df
        
        # Rule 1: near-full > composite > partial > weak
        coverage_order = {'near-full': 0, 'composite': 1, 'partial': 2, 'weak': 3}
        df['coverage_priority'] = df['coverage_type'].map(coverage_order)
        
        # Rule 2: megablast > dc-megablast > blastn
        task_order = {'megablast': 0, 'dc-megablast': 1, 'blastn': 2}
        df['task_priority'] = df['task'].map(task_order)
        
        # Rule 3: x1 > x3
        df['x3_priority'] = df['db'].apply(lambda x: 1 if '_x3' in x else 0)
        
        # Rule 4: subjects sort order
        df['subject_priority'] = df['qseqid'].apply(lambda x: 
            next((i for i, s in enumerate(subjects) if x.startswith(s)), len(subjects)))
        
        # Sort by all priorities
        df = df.sort_values([
            'coverage_priority', 'task_priority', 'x3_priority', 
            'subject_priority', 'evalue'
        ]).reset_index(drop=True)
        
        # Drop priority columns
        df = df.drop(['coverage_priority', 'task_priority', 'x3_priority', 'subject_priority'], axis=1)
        
        return df
    
    def _apply_levenshtein_filter(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply Levenshtein distance filtering to remove similar NCBI accessions"""
        if df.empty:
            return df
        
        # Check if filtering is enabled
        if self.config and not self.config.filtering.remove_self_blasts:
            return df
        
        # Only apply to NCBI database results
        ncbi_df = df[df['db'].str.contains('ncbi', case=False)]
        other_df = df[~df['db'].str.contains('ncbi', case=False)]
        
        if ncbi_df.empty:
            return df
        
        logger.info("Applying Levenshtein distance filtering to NCBI results...")
        
        # Get threshold from config
        threshold = 1  # Default
        if self.config:
            threshold = int(100 - self.config.filtering.levenshtein_threshold)  # Convert percentage to distance
        
        filtered_ncbi = []
        for qseqid, group in ncbi_df.groupby('qseqid'):
            # Get unique subject names
            subject_names = sorted(group['sseqid'].unique())
            
            if len(subject_names) <= 1:
                filtered_ncbi.append(group)
                continue
            
            # Keep the first name and filter similar ones
            names_to_keep = [subject_names[0]]
            for name in subject_names[1:]:
                # Check if this name is too similar to any kept name
                too_similar = False
                for kept_name in names_to_keep:
                    distance = self._levenshtein_distance(kept_name, name)
                    if distance <= threshold:  # Use config threshold
                        too_similar = True
                        break
                
                if not too_similar:
                    names_to_keep.append(name)
            
            # Filter the group to keep only the selected names
            filtered_group = group[group['sseqid'].isin(names_to_keep)]
            filtered_ncbi.append(filtered_group)
        
        if filtered_ncbi:
            filtered_ncbi_df = pd.concat(filtered_ncbi, ignore_index=True)
            result_df = pd.concat([other_df, filtered_ncbi_df], ignore_index=True)
            logger.info(f"Levenshtein filtering: {len(ncbi_df)} -> {len(filtered_ncbi_df)} NCBI hits")
            return result_df
        else:
            return df
    
    def _levenshtein_distance(self, seq1: str, seq2: str) -> int:
        """Calculate Levenshtein distance between two strings"""
        size_x = len(seq1) + 1
        size_y = len(seq2) + 1
        matrix = np.zeros((size_x, size_y))
        
        for x in range(size_x):
            matrix[x, 0] = x
        for y in range(size_y):
            matrix[0, y] = y
        
        for x in range(1, size_x):
            for y in range(1, size_y):
                if seq1[x-1] == seq2[y-1]:
                    matrix[x, y] = min(
                        matrix[x-1, y] + 1,
                        matrix[x-1, y-1],
                        matrix[x, y-1] + 1
                    )
                else:
                    matrix[x, y] = min(
                        matrix[x-1, y] + 1,
                        matrix[x-1, y-1] + 1,
                        matrix[x, y-1] + 1
                    )
        
        return int(matrix[size_x - 1, size_y - 1])
    
    def _filter_by_e_value(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter results by E-value based on coverage type"""
        if df.empty:
            return df
        
        # Check if filtering is enabled
        if self.config and not self.config.filtering.filter_weak_by_evalue:
            return df
        
        # Get threshold from config
        threshold = 1.0  # Default
        if self.config:
            threshold = self.config.filtering.weak_evalue_threshold
        
        # For weak or partial alignments, require E-value < threshold
        mask = ~((df['coverage_type'].isin(['weak', 'partial'])) & (df['evalue'] >= threshold))
        return df[mask].reset_index(drop=True)
    
    def _remove_self_blasts(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove self-blast results"""
        if df.empty:
            return df
        
        # Check if filtering is enabled
        if self.config and not self.config.filtering.remove_self_blasts:
            return df
        
        # Remove where query and subject are the same
        return df[df['qseqid'] != df['sseqid']].reset_index(drop=True)
    
    def _stairway_view(self, df: pd.DataFrame) -> pd.DataFrame:
        """Create stairway view for multiple subjects"""
        if df.empty:
            return df
        
        # Group by subject and keep best hit per query
        best_hits = []
        
        for sseqid, group in df.groupby('sseqid'):
            best_hit = group.loc[group['evalue'].idxmin()]
            best_hits.append(best_hit)
        
        return pd.DataFrame(best_hits)
    
    def run_blast(self, subjects: List[str], dbs: List[str] = None, 
                  tasks: List[str] = None, num_threads: int = None, 
                  verbose: bool = None) -> pd.DataFrame:
        """Main method to run BLAST analysis"""
        # Use config values if available and not provided
        if self.config:
            if dbs is None:
                dbs = self.config.blast.default_dbs
            if tasks is None:
                tasks = self.config.blast.default_tasks
            if num_threads is None:
                num_threads = self.config.blast.default_threads
            if verbose is None:
                verbose = self.config.blast.verbose
        else:
            # Fallback to defaults
            if dbs is None:
                dbs = ['local']
            if tasks is None:
                tasks = ['megablast', 'dc-megablast', 'blastn']
            if num_threads is None:
                num_threads = 20
            if verbose is None:
                verbose = True
        
        logger.info(f"Starting BLAST analysis for subjects: {subjects}")
        logger.info(f"Using databases: {dbs}")
        logger.info(f"Using tasks: {tasks}")
        
        fmt_header = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                     'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        
        try:
            # ------------------------------------------------------------------
            # Query construction strategy
            # ------------------------------------------------------------------
            #
            # For post-TAREAN integration (where a ConfigManager is present),
            # we treat `subjects` as sample IDs and build a dedicated query
            # FASTA from TAREAN consensus sequences. The BLAST databases
            # (e.g. the Triticeae NCBI-derived repeat DB) are queried as-is,
            # without attempting to slice them by KA* patterns.
            #
            # For legacy / standalone usage (no ConfigManager attached), we
            # preserve the historical behaviour of slicing local FASTA DBs
            # based on `subjects`-derived regex patterns.
            # ------------------------------------------------------------------

            if self.config_manager is not None:
                # New, config-driven query construction path.
                with tempfile.TemporaryDirectory() as tmpdirname:
                    input_fasta_abs = self._build_query_fasta_for_subjects(
                        subjects=subjects,
                        tmpdirname=tmpdirname,
                    )
                    if not input_fasta_abs:
                        logger.warning("BLAST skipped: no query FASTA could be built.")
                        return pd.DataFrame()

                    # Run BLAST analysis against the configured databases.
                    df = self._blast_coordinator(
                        tasks=tasks,
                        input_fasta_abs=input_fasta_abs,
                        num_threads=num_threads,
                        verbose=verbose,
                        fmt_header=fmt_header,
                        dbs=dbs,
                    )
            else:
                # Legacy self/local mode – slice local FASTA DBs by patterns.
                # Check if required databases are available; degrade gracefully if not.
                missing_dbs = [db for db in dbs if db not in self.db_dict]
                if missing_dbs:
                    logger.warning(
                        f"Skipping BLAST analysis because databases are missing: {missing_dbs}"
                    )
                    return pd.DataFrame()

                # Create patterns for sequence extraction
                patterns = [self._protect_brackets(x) + '_.*' for x in subjects]
                logger.info(f"Search patterns: {patterns}")

                # Extract sequences from databases into a temporary FASTA.
                with tempfile.TemporaryDirectory() as tmpdirname:
                    output_path = os.path.join(
                        tmpdirname,
                        f'{"_".join(patterns)}-sequences.fasta',
                    )
                    self._slice_databases(patterns, dbs, output_path)
                    input_fasta_abs = os.path.abspath(output_path)

                    df = self._blast_coordinator(
                        tasks=tasks,
                        input_fasta_abs=input_fasta_abs,
                        num_threads=num_threads,
                        verbose=verbose,
                        fmt_header=fmt_header,
                        dbs=dbs,
                    )
            
            # Clean and process results
            df = self._clean_dataframe(df, tasks)
            
            if df.empty:
                logger.warning("BLAST returned no hits, skipping downstream processing.")
                return pd.DataFrame()
            
            # Process each query-subject-task-db combination
            processed_dfs = []
            total_combinations = len(df.groupby(['qseqid', 'sseqid', 'task', 'db']))
            
            logger.info(f"Processing {total_combinations} query-subject combinations...")
            
            # Get evalue threshold from config
            evalue_max_olig = 0.001  # Default
            if self.config:
                evalue_max_olig = self.config.blast.evalue_max_olig
            
            for (qseqid, sseqid, task, db), group in df.groupby(['qseqid', 'sseqid', 'task', 'db']):
                processed_group = self._return_best_examples(group, evalue_max_olig=evalue_max_olig)
                processed_dfs.append(processed_group)
            
            df = pd.concat(processed_dfs, ignore_index=True)
            
            # Apply final filters and sorting
            df = self._sort_for_rules(df, subjects)
            df = self._remove_self_blasts(df)
            df = self._filter_by_e_value(df)
            
            # Apply Levenshtein filtering if available
            if FUZZ_AVAILABLE:
                df = self._apply_levenshtein_filter(df)
            
            # Create stairway view for multiple subjects
            if len(subjects) > 1:
                df = self._stairway_view(df)
            
            logger.info(f"BLAST analysis completed. Found {len(df)} significant hits.")
            return df
            
        except Exception as e:
            logger.error(f"Error during BLAST analysis: {e}")
            # Degrade gracefully on BLAST-level failures.
            return pd.DataFrame()


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Run consolidated BLAST analysis (dev helper)")
    parser.add_argument("subject", nargs="?", default="sample", help="Subject/sample identifier")
    parser.add_argument("--db", dest="dbs", action="append", default=["local"], help="DB name (repeatable)")
    parser.add_argument("--task", dest="tasks", action="append", default=["megablast"], help="BLAST task (repeatable)")
    args = parser.parse_args()

    analyzer = BLASTAnalyzer()
    results = analyzer.run_blast([args.subject], dbs=args.dbs, tasks=args.tasks)
    print(f"Found {len(results)} BLAST hits")
    try:
        print(results.head())
    except Exception:
        print(results)


if __name__ == "__main__":
    main() 