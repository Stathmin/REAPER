#!/usr/bin/env python3
"""
Prepare comparative reads for analysis with genome size-based proportions
Calculates read ratios automatically from genome sizes to ensure equal coverage per genome size

HOLY WORKFLOW COMPLIANCE: This module supports the modular workflow's comparative analysis framework.
All configuration comes from projects/global_config.yaml.
"""

import os
import sys
import argparse
import yaml
import subprocess
import tempfile
from pathlib import Path
import logging

# Ensure repo root is on sys.path when invoked as a script.
# Snakemake executes this via `python3 workflows/smk_scripts/prepare_comparative_reads.py ...`,
# so `workflows.*` imports require the repository root to be importable.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflows.smk_scripts.iterative_kmer_filter_paired_fastq import filter_paired_reads_to_unmapped

from workflows.smk_scripts.reportr_config import resolve_value as resolve_cfg_value
from workflows.smk_scripts.reportr_config import sample_metadata as sample_cfg

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_config():
    """Load full configuration YAML."""
    config_path = "projects/global_config.yaml"
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def calculate_read_ratios(config, project_id, comparative_analysis_id):
    """
    Calculate read ratios based on genome sizes to ensure equal coverage per genome size
    
    Args:
        project_config: Project configuration dictionary
        comparative_analysis_id: ID of the comparative analysis
    
    Returns:
        dict: Sample name -> read count mapping
    """
    # Get the comparative analysis configuration
    project_config = config["projects"][project_id]
    comparative_analysis = project_config["comparative_analyses"][comparative_analysis_id]
    sample_names = comparative_analysis["samples"]
    
    # Get genome sizes for each sample
    genome_sizes = {}
    for sample_name in sample_names:
        genome_sizes[sample_name] = sample_cfg(config, project=project_id, sample=sample_name)["genome_size"]
    
    # In RepOrtR, total_reads_per_assembly is treated as total FASTA records
    # (i.e., mates are separate records). For paired-end interleaved FASTA,
    # the number of *pairs* is total_reads_per_assembly // 2.
    total_reads_per_assembly = int(
        resolve_cfg_value(
            config,
            project=project_id,
            comparative_analysis=comparative_analysis_id,
            allow_sample_overrides=False,
            path=["total_reads_per_assembly"],
        )
    )
    total_pairs_per_assembly = max(1, total_reads_per_assembly // 2)
    
    # Calculate the target coverage per genome size unit
    # We want equal coverage per genome size across all samples
    total_genome_size = sum(genome_sizes.values())
    
    # Allocate PAIRS proportional to genome size, then convert to FASTA records
    # (2 records per pair) to preserve f/r alternation.
    pair_counts = {}
    for sample_name in sample_names:
        genome_size = genome_sizes[sample_name]
        sample_pairs = int((genome_size / total_genome_size) * total_pairs_per_assembly)
        pair_counts[sample_name] = sample_pairs
    
    # Ensure total pairs equals the target
    total_pairs_calculated = sum(pair_counts.values())
    if total_pairs_calculated != total_pairs_per_assembly:
        largest_sample = max(pair_counts.keys(), key=lambda x: pair_counts[x])
        adjustment = total_pairs_per_assembly - total_pairs_calculated
        pair_counts[largest_sample] += adjustment

    read_counts = {k: v * 2 for k, v in pair_counts.items()}
    
    logger.info(f"Calculated read ratios for {comparative_analysis_id}:")
    for sample_name in sample_names:
        genome_size = genome_sizes[sample_name]
        pairs = pair_counts[sample_name]
        reads = read_counts[sample_name]
        coverage_per_genome = pairs / genome_size
        logger.info(
            f"  {sample_name}: {pairs} pairs ({reads} FASTA records), genome_size={genome_size}, "
            f"pairs/genome={coverage_per_genome:.2f}"
        )
    
    return read_counts

def prepare_comparative_reads(
    project_id,
    comparative_analysis_id,
    seed=None,
    outdir=None,
    filter_unmapped_against=None,
    filter_threads=8,
    kmer_k=27,
    minkmerhits=1,
    removeifeitherbad=True,
    keep_sampled_fastq=True,
):
    """
    Prepare comparative reads with genome size-based proportions

    When ``filter_unmapped_against`` is set, k-mer filtering uses the same
    pipeline as ``deconseq_filter_reads_iter`` (see ``iterative_kmer_filter_paired_fastq``).

    Args:
        project_id: Project identifier
        comparative_analysis_id: ID of the comparative analysis to run
    """
    logger.info(f"Preparing comparative reads for project {project_id}, analysis {comparative_analysis_id}")
    
    # Load configuration
    config = load_config()
    project_config = config["projects"][project_id]
    
    # Calculate read ratios based on genome sizes
    read_counts = calculate_read_ratios(config, project_id, comparative_analysis_id)
    
    # Get sample names for this comparative analysis
    comparative_cfg = project_config["comparative_analyses"][comparative_analysis_id]
    sample_names = comparative_cfg["samples"]

    # Log prefix policy (comparative headers).
    base_prefixes = {s: str(project_config["samples"][s].get("prefix") or str(s).upper()) for s in sample_names}
    longest = max((len(p) for p in base_prefixes.values()), default=0)
    cfg_floor = int(
        ((comparative_cfg.get("tarean_params", {}) or {}).get("prefix_length"))
        or ((project_config.get("tarean_params", {}) or {}).get("prefix_length"))
        or 0
    )
    prefix_len = max(longest, cfg_floor)
    padded_prefixes = {s: p.ljust(prefix_len, "_") for s, p in base_prefixes.items()}
    logger.info(f"Comparative prefix_len={prefix_len} (longest={longest}, cfg_floor={cfg_floor})")
    logger.info("Comparative padded prefixes: " + ", ".join(f"{s}={p}" for s, p in padded_prefixes.items()))

    # Strict comparative operational parameters: resolved from
    # global → project → comparative_analyses.<id> only.
    if filter_unmapped_against:
        filter_threads = int(
            resolve_cfg_value(
                config,
                project=project_id,
                comparative_analysis=comparative_analysis_id,
                allow_sample_overrides=False,
                path=["deconseq", "threads"],
            )
        )
        kmer_k = int(
            resolve_cfg_value(
                config,
                project=project_id,
                comparative_analysis=comparative_analysis_id,
                allow_sample_overrides=False,
                path=["deconseq", "kmer_k"],
            )
        )
        minkmerhits = int(
            resolve_cfg_value(
                config,
                project=project_id,
                comparative_analysis=comparative_analysis_id,
                allow_sample_overrides=False,
                path=["deconseq", "minkmerhits"],
            )
        )
        removeifeitherbad = bool(
            resolve_cfg_value(
                config,
                project=project_id,
                comparative_analysis=comparative_analysis_id,
                allow_sample_overrides=False,
                path=["deconseq", "removeifeitherbad"],
            )
        )
    
    # Create output directory (allow iterative runs under a custom directory)
    comparative_dir = outdir or f"projects/{project_id}/comparative/{comparative_analysis_id}"
    os.makedirs(comparative_dir, exist_ok=True)

    # Write combined reads incrementally to avoid large in-memory lists.
    total_reads = 0
    per_sample_stats = []
    filter_meta = None
    if filter_unmapped_against:
        filter_meta = {
            "ref": str(filter_unmapped_against),
            "kmer_k": int(kmer_k),
            "threads": int(filter_threads),
            "minkmerhits": int(minkmerhits),
            "removeifeitherbad": bool(removeifeitherbad),
        }

    output_file = f"{comparative_dir}/comparative_reads.fasta"
    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as out_f:
        for sample_name in sample_names:
            target_reads = int(read_counts[sample_name])
            target_pairs = max(1, target_reads // 2)
            logger.info(f"Processing {sample_name}: target {target_pairs} pairs ({target_reads} FASTA records)")

            sample_prefix = padded_prefixes[sample_name]

            # Comparative preparation MUST use cleaned reads.
            # Do not fall back to raw_reads: cleaned FASTQs provide consistent quality gating
            # and avoid raw-format edge cases.
            r1 = Path(f"projects/{project_id}/samples/{sample_name}/filtered_reads/R1.fq")
            r2 = Path(f"projects/{project_id}/samples/{sample_name}/filtered_reads/R2.fq")
            if not (r1.exists() and r2.exists()):
                raise FileNotFoundError(
                    f"Missing cleaned reads for sample {sample_name}: expected {r1} and {r2}. "
                    "Run the clean_reads step first. raw_reads are not accepted for comparatives."
                )

            sample_tmp = Path(comparative_dir) / "tmp" / sample_name
            sample_tmp.mkdir(parents=True, exist_ok=True)
            tmp_fa = sample_tmp / f"{sample_name}_comparative_tmp.fa"

            # Optional iterative filtering: same pipeline as deconseq_filter_reads_iter.
            if filter_unmapped_against:
                ref = Path(filter_unmapped_against)
                if not ref.exists():
                    raise FileNotFoundError(f"filter_unmapped_against not found: {ref}")
                filt_r1 = sample_tmp / "unmapped_R1.fq"
                filt_r2 = sample_tmp / "unmapped_R2.fq"
                with tempfile.TemporaryDirectory(dir=str(sample_tmp), prefix="deconseq_") as td:
                    filter_paired_reads_to_unmapped(
                        r1,
                        r2,
                        ref,
                        filt_r1,
                        filt_r2,
                        Path(td),
                        kmer_k=int(kmer_k),
                        threads=int(filter_threads),
                        minkmerhits=int(minkmerhits),
                        removeifeitherbad=bool(removeifeitherbad),
                    )
                r1, r2 = filt_r1, filt_r2

            # Match solo prepare logic: deterministic reformat.sh subsampling.
            # Keep the provided seed unchanged across all seed-bearing paths.
            sample_seed = int(seed) if seed is not None else 0
            if keep_sampled_fastq:
                sampled_fastq_dir = Path(comparative_dir) / "sampled_fastq" / sample_name
                sampled_fastq_dir.mkdir(parents=True, exist_ok=True)
                sampled_r1 = sampled_fastq_dir / "R1.fq"
                sampled_r2 = sampled_fastq_dir / "R2.fq"

                # Step 1: sample paired FASTQ (used for FastQC after sampling)
                cmd = [
                    "reformat.sh",
                    f"in1={r1}",
                    f"in2={r2}",
                    f"out1={sampled_r1}",
                    f"out2={sampled_r2}",
                    f"samplereadstarget={target_pairs}",
                    f"sampleseed={sample_seed}",
                    "fastawrap=0",
                    "ow=t",
                ]
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode != 0:
                    raise RuntimeError(
                        f"reformat.sh failed for {sample_name} (exit {proc.returncode}):\n"
                        f"{proc.stderr}"
                    )

                # Step 2: convert sampled FASTQ pair to interleaved FASTA (no resampling)
                cmd = [
                    "reformat.sh",
                    f"in1={sampled_r1}",
                    f"in2={sampled_r2}",
                    f"out={tmp_fa}",
                    "interleaved=t",
                    "ow=t",
                ]
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode != 0:
                    raise RuntimeError(
                        f"reformat.sh FASTA conversion failed for {sample_name} (exit {proc.returncode}):\n"
                        f"{proc.stderr}"
                    )
            else:
                # Directly subsample into FASTA (skip writing sampled FASTQ files).
                cmd = [
                    "reformat.sh",
                    f"in1={r1}",
                    f"in2={r2}",
                    f"out={tmp_fa}",
                    "interleaved=t",
                    f"samplereadstarget={target_pairs}",
                    f"sampleseed={sample_seed}",
                    "fastawrap=0",
                    "ow=t",
                ]
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode != 0:
                    raise RuntimeError(
                        f"reformat.sh FASTA subsample failed for {sample_name} (exit {proc.returncode}):\n"
                        f"{proc.stderr}"
                    )

            # Rename headers exactly like scripts/prepare_reads.sh (odd=f, even=r).
            emitted = 0
            current_lines: list[str] = []
            rec_idx = 0
            with open(tmp_fa, "r") as f:
                for line in f:
                    if line.startswith(">") and current_lines:
                        rec_idx += 1
                        if emitted < target_reads:
                            suffix = "f" if rec_idx % 2 == 1 else "r"
                            pair_idx = (rec_idx + 1) // 2 if rec_idx % 2 == 1 else rec_idx // 2
                            header = f">{sample_prefix}read{pair_idx}_{suffix}\n"
                            out_f.write(header)
                            out_f.writelines(current_lines[1:])
                            emitted += 1
                        current_lines = []
                    current_lines.append(line)
                if current_lines:
                    rec_idx += 1
                    if emitted < target_reads:
                        suffix = "f" if rec_idx % 2 == 1 else "r"
                        pair_idx = (rec_idx + 1) // 2 if rec_idx % 2 == 1 else rec_idx // 2
                        header = f">{sample_prefix}read{pair_idx}_{suffix}\n"
                        out_f.write(header)
                        out_f.writelines(current_lines[1:])
                        emitted += 1

            if emitted < target_reads:
                raise RuntimeError(f"{sample_name}: prepared fewer records than requested ({emitted} < {target_reads})")
            if emitted % 2 != 0:
                raise RuntimeError(f"{sample_name}: odd emitted record count {emitted} (broken pairing)")

            total_reads += emitted
            logger.info(f"Added {emitted//2} pairs ({emitted} records) from {sample_name}")
            per_sample_stats.append(
                {
                    "sample": sample_name,
                    "target_records": target_reads,
                    "emitted_records": emitted,
                }
            )
    
    logger.info(f"Comparative assembly prepared: {total_reads} total reads written to {output_file}")
    return {
        "output_file": output_file,
        "per_sample": per_sample_stats,
        "filter": filter_meta,
    }

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Prepare comparative reads with genome size-based proportions")
    parser.add_argument("project_id", help="Project identifier")
    parser.add_argument("comparative_analysis_id", help="Comparative analysis identifier")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--seed", type=int, default=None, help="Deterministic RNG seed for sampling")
    parser.add_argument("--outdir", default=None, help="Override output directory (for iterative runs)")
    parser.add_argument(
        "--write-summary",
        default=None,
        help="If set, write a small text summary to this path.",
    )
    parser.add_argument(
        "--write-token",
        default=None,
        help="If set, write a preparation-complete token to this path.",
    )
    parser.add_argument(
        "--filter-unmapped-against",
        default=None,
        help="FASTA reference to filter reads against (keep unmapped) before sampling",
    )
    parser.add_argument("--filter-threads", type=int, default=8, help="Threads for filtering step")
    parser.add_argument(
        "--kmer-k",
        type=int,
        default=27,
        help="K-mer size for bbduk when --filter-unmapped-against is set (match deconseq.kmer_k)",
    )
    parser.add_argument(
        "--keep-sampled-fastq",
        action="store_true",
        help="Write comparative sampled_fastq/<sample>/R1.fq,R2.fq (needed for comparative_fastqc_sampled).",
    )
    parser.add_argument(
        "--minkmerhits",
        type=int,
        default=1,
        help="bbduk minkmerhits when filtering (match deconseq.minkmerhits)",
    )
    parser.add_argument(
        "--removeifeitherbad",
        choices=("t", "f"),
        default="t",
        help="bbduk removeifeitherbad: t=stricter (either mate matches removes pair), f=lenient",
    )

    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    try:
        result = prepare_comparative_reads(
            args.project_id,
            args.comparative_analysis_id,
            seed=args.seed,
            outdir=args.outdir,
            filter_unmapped_against=args.filter_unmapped_against,
            filter_threads=args.filter_threads,
            kmer_k=args.kmer_k,
            minkmerhits=args.minkmerhits,
            removeifeitherbad=(args.removeifeitherbad == "t"),
            keep_sampled_fastq=bool(args.keep_sampled_fastq),
        )
        output_file = result["output_file"]
        if args.write_summary:
            Path(os.path.dirname(args.write_summary)).mkdir(parents=True, exist_ok=True)
            with open(args.write_summary, "w") as f:
                f.write("Comparative Analysis Summary (iterative)\n")
                f.write(f"Project: {args.project_id}\n")
                f.write(f"Analysis: {args.comparative_analysis_id}\n")
                f.write(f"Output: {output_file}\n")
                if args.seed is not None:
                    f.write(f"Seed: {args.seed}\n")
                fm = result.get("filter")
                if fm:
                    f.write(f"FilterUnmappedAgainst: {fm['ref']}\n")
                    f.write(f"KmerK: {fm['kmer_k']}\n")
                    f.write(f"FilterThreads: {fm['threads']}\n")
                    f.write(f"Minkmerhits: {fm.get('minkmerhits', 1)}\n")
                    f.write(f"RemoveIfEitherBad: {fm.get('removeifeitherbad', True)}\n")
                f.write("Per-sample targets vs emitted FASTA records:\n")
                for row in result.get("per_sample", []):
                    f.write(
                        f"  {row['sample']}: target={row['target_records']} emitted={row['emitted_records']}\n"
                    )
        if args.write_token:
            Path(os.path.dirname(args.write_token)).mkdir(parents=True, exist_ok=True)
            Path(args.write_token).write_text("PREPARATION_READY\n")
        print(f"Comparative reads prepared successfully: {output_file}")
    except Exception as e:
        logger.error(f"Error preparing comparative reads: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 