#!/usr/bin/env bash
set -euo pipefail

project="$1"
sample="$2"
sample_prefix="$3"
pythonhashseed="$4"
temp_dir="$5"
reads_per_assembly="$6"
max_retries="$7"
in_r1="$8"
in_r2="$9"
prepared_out="${10}"
tarean_token="${11}"
sampled_fastqc_html1="${12}"
sampled_fastqc_html2="${13}"
fastqc_threads="${14}"
log_path="${15}"

LOG="$log_path"

mkdir -p "$(dirname "$prepared_out")"
mkdir -p "$(dirname "$LOG")"
mkdir -p "$temp_dir"
mkdir -p "$(dirname "$sampled_fastqc_html1")"

# Use a per-job temp directory under the configured base temp root.
# This avoids collisions when Snakemake runs jobs in parallel.
if [[ "$temp_dir" = /* ]]; then
  _base_tmp_root="$temp_dir"
else
  _base_tmp_root="$(pwd)/${temp_dir}"
fi

_job_tmp_parent="${_base_tmp_root}/jobs/prepare_reads/${project}/${sample}"
mkdir -p "$_job_tmp_parent"
job_tmp="$(mktemp -d "${_job_tmp_parent}/run_XXXXXXXX")"

_cleanup_job_tmp() {
  local exit_code="$?"
  if [[ "$exit_code" -eq 0 ]]; then
    rm -rf "$job_tmp" || true
    # Best-effort cleanup of empty job parent dirs so tmp/ is removable.
    rmdir "$_job_tmp_parent" 2>/dev/null || true
    rmdir "$(dirname "$_job_tmp_parent")" 2>/dev/null || true
    rmdir "$(dirname "$(dirname "$_job_tmp_parent")")" 2>/dev/null || true
    rmdir "$(dirname "$(dirname "$(dirname "$_job_tmp_parent")")")" 2>/dev/null || true
  else
    echo "$(date): PREPARE_READS_TMPDIR_RETAINED=${job_tmp}" >>"$LOG" || true
  fi
}
trap _cleanup_job_tmp EXIT

{
  echo "$(date): PREPARE_READS start for project=${project}, sample=${sample}, prefix=${sample_prefix}"
  echo "$(date): Input R1=${in_r1}"
  echo "$(date): Input R2=${in_r2}"
  echo "$(date): Output FASTA=${prepared_out}"
  echo "$(date): Sampled FastQC outdir=$(dirname "$sampled_fastqc_html1") (threads=${fastqc_threads})"
  echo "$(date): TEMP base dir=${temp_dir}"
  echo "$(date): JOB TMP dir=${job_tmp}"
  echo "$(date): reads_per_assembly=${reads_per_assembly}, max_retries=${max_retries}"
} >>"$LOG"

export TMPDIR="$job_tmp"
export TMP="$job_tmp"
export TEMP="$job_tmp"
# Keep seed deterministic for any Python tooling invoked from this script.
export PYTHONHASHSEED="$pythonhashseed"
mkdir -p "$job_tmp"

# FastQC output normalization: make rule-declared filenames deterministic
_fastqc_out_stem() {
  local fname
  fname="$(basename "$1")"
  if [[ "$fname" == *.gz ]]; then
    fname="${fname%.gz}"
  fi
  for ext in .fastq .fq .fasta .fa .txt; do
    if [[ "$fname" == *"$ext" ]]; then
      fname="${fname%$ext}"
      break
    fi
  done
  printf '%s' "$fname"
}

_move_fastqc_outputs() {
  local in_path="$1"
  local out_html="$2"
  local out_zip="${out_html%.html}.zip"
  local qc_dir
  qc_dir="$(dirname "$out_html")"
  local stem
  stem="$(_fastqc_out_stem "$in_path")"
  local produced_html="${qc_dir}/${stem}_fastqc.html"
  local produced_zip="${qc_dir}/${stem}_fastqc.zip"

  if [[ -f "$produced_html" && "$produced_html" != "$out_html" ]]; then
    mv -f "$produced_html" "$out_html"
  fi
  if [[ -f "$produced_zip" && "$produced_zip" != "$out_zip" ]]; then
    mv -f "$produced_zip" "$out_zip"
  fi
}

# Resume capability
if [[ -f "$tarean_token" ]]; then
  {
    echo "$(date): PREPARE_READS token found, skipping preparation"
    echo "PREPARE_READS_SUCCESS: token_exists=1"
  } >>"$LOG"
  exit 0
fi

retry_count=0

while [[ $retry_count -lt $max_retries ]]; do
  echo "$(date): Attempt $((retry_count + 1)) of $max_retries" >>"$LOG"

  # RepOrtR convention: reads_per_assembly means FASTA record count in prepared_forRE.fasta.
  # For paired-end interleaved FASTA, pairs = records / 2 and record count must be even.
  if (( reads_per_assembly % 2 != 0 )); then
    {
      echo "$(date): PREPARE_READS_ERROR: reads_per_assembly must be even for paired-end inputs."
      echo "PREPARE_READS_ERROR: odd_reads_per_assembly=${reads_per_assembly}"
    } >>"$LOG"
    exit 1
  fi
  target_records="$reads_per_assembly"
  target_pairs="$(( reads_per_assembly / 2 ))"

  echo "$(date): Running reformat.sh with target_pairs=${target_pairs} (target_records=${target_records})" >>"$LOG"
  tmp_fa="${job_tmp}/tmp_${project}_${sample}.fa"
  # IMPORTANT: temp filenames must be unique per job to avoid parallel collisions.
  sampled_r1="${job_tmp}/SAMPLED_${project}_${sample}_R1.fq.gz"
  sampled_r2="${job_tmp}/SAMPLED_${project}_${sample}_R2.fq.gz"

  # First: sample paired FASTQ for post-sampling FastQC
  if ! reformat.sh \
      in1="$in_r1" \
      in2="$in_r2" \
      out1="$sampled_r1" \
      out2="$sampled_r2" \
      samplereadstarget="$target_pairs" \
      sampleseed="$pythonhashseed" \
      ow=t >>"$LOG" 2>&1; then
    retry_count=$((retry_count + 1))
    echo "$(date): reformat.sh failed during sampled FASTQ generation (attempt ${retry_count}/${max_retries})" >>"$LOG"

    if [[ $retry_count -lt $max_retries ]]; then
      echo "$(date): Waiting 30 seconds before retry..." >>"$LOG"
      sleep 30
      continue
    else
      {
        echo "$(date): All retry attempts failed during sampled FASTQ generation."
        echo "PREPARE_READS_ERROR: step=reformat_sampled_fastq exit_code=1"
      } >>"$LOG"
      exit 1
    fi
  fi

  # Run FastQC on the sampled subset
  fastqc_outdir="$(dirname "$sampled_fastqc_html1")"
  {
    echo "$(date): PREPARE_READS running fastqc on sampled reads"
    echo "COMMAND: fastqc -t ${fastqc_threads} ${sampled_r1} ${sampled_r2} -o ${fastqc_outdir}"
  } >>"$LOG"
  if ! fastqc -t "$fastqc_threads" "$sampled_r1" "$sampled_r2" -o "$fastqc_outdir" --quiet >>"$LOG" 2>&1; then
    rc=$?
    {
      echo "$(date): PREPARE_READS_ERROR: fastqc failed on sampled reads"
      echo "PREPARE_READS_ERROR: step=fastqc_sampled exit_code=${rc}"
    } >>"$LOG"
    exit 1
  fi

  # Normalize FastQC outputs to the rule-declared paths (SAMPLED_R1_fastqc.html, etc.)
  _move_fastqc_outputs "$sampled_r1" "$sampled_fastqc_html1"
  _move_fastqc_outputs "$sampled_r2" "$sampled_fastqc_html2"

  # Then: convert sampled FASTQ pair to interleaved FASTA (no resampling)
  if ! reformat.sh \
      in1="$sampled_r1" \
      in2="$sampled_r2" \
      out="$tmp_fa" \
      interleaved=t \
      ow=t >>"$LOG" 2>&1; then
    retry_count=$((retry_count + 1))
    echo "$(date): reformat.sh failed during FASTA conversion (attempt ${retry_count}/${max_retries})" >>"$LOG"

    if [[ $retry_count -lt $max_retries ]]; then
      echo "$(date): Waiting 30 seconds before retry..." >>"$LOG"
      sleep 30
      continue
    else
      {
        echo "$(date): All retry attempts failed during FASTA conversion."
        echo "PREPARE_READS_ERROR: step=reformat_fasta_conversion exit_code=1"
      } >>"$LOG"
      exit 1
    fi
  fi

  # Optional cleanup of sampled FASTQs to conserve space
  rm -f "$sampled_r1" "$sampled_r2" >>"$LOG" 2>&1 || true

  # Header renaming
  if awk -v p="$sample_prefix" '
    /^>/ {
      c++
      if (c % 2 == 1) {
        print ">" p "read" int((c+1)/2) "_f"
      } else {
        print ">" p "read" int(c/2) "_r"
      }
      next
    }
    { print }
  ' "$tmp_fa" >"$prepared_out" 2>>"$LOG"; then
    rm -f "$tmp_fa"

    # Validate that we produced the expected number of FASTA records.
    # This is critical when total_reads_per_assembly changes: it forces failures
    # (and therefore reruns) instead of silently producing undersized assemblies.
    actual_records="$(grep -c '^>' "$prepared_out" || true)"
    expected_records="$target_records"
    {
      echo "$(date): PREPARE_READS validate: expected_records=${expected_records} actual_records=${actual_records}"
    } >>"$LOG"

    if [[ "$actual_records" -lt "$expected_records" ]]; then
      {
        echo "$(date): PREPARE_READS_ERROR: prepared_forRE.fasta has fewer records than requested."
        echo "PREPARE_READS_ERROR: expected_records=${expected_records} actual_records=${actual_records}"
      } >>"$LOG"
      exit 1
    fi

    # For paired-end interleaved FASTA, record count should be even.
    if (( actual_records % 2 != 0 )); then
      {
        echo "$(date): PREPARE_READS_ERROR: prepared_forRE.fasta record count is odd (broken pairing)."
        echo "PREPARE_READS_ERROR: odd_record_count=${actual_records}"
      } >>"$LOG"
      exit 1
    fi

    # If reformat produced more, keep exact target while preserving pair alternation.
    if [[ "$actual_records" -gt "$expected_records" ]]; then
      {
        echo "$(date): PREPARE_READS warning: prepared_forRE.fasta has more records than requested; truncating."
        echo "PREPARE_READS_TRUNCATE: expected_records=${expected_records} actual_records=${actual_records}"
      } >>"$LOG"
      awk -v n="$expected_records" '
        /^>/ { h++ }
        { if (h<=n) print }
      ' "$prepared_out" > "${prepared_out}.tmp"
      mv -f "${prepared_out}.tmp" "$prepared_out"
      actual_records="$(grep -c '^>' "$prepared_out" || true)"
      {
        echo "$(date): PREPARE_READS truncate validate: expected_records=${expected_records} actual_records=${actual_records}"
      } >>"$LOG"
      if [[ "$actual_records" -ne "$expected_records" ]]; then
        {
          echo "$(date): PREPARE_READS_ERROR: truncation did not yield exact record count."
          echo "PREPARE_READS_ERROR: expected_records=${expected_records} actual_records=${actual_records}"
        } >>"$LOG"
        exit 1
      fi
    fi

    {
      echo "$(date): Read preparation completed successfully with prefix ${sample_prefix}"
      echo "PREPARE_READS_SUCCESS: exit_code=0"
    } >>"$LOG"
    echo "PREPARATION_READY: $(date)" >"$tarean_token"
    exit 0
  else
    rm -f "$tmp_fa"
    retry_count=$((retry_count + 1))
    echo "$(date): Header renaming failed (attempt ${retry_count}/${max_retries})" >>"$LOG"

    if [[ $retry_count -lt $max_retries ]]; then
      echo "$(date): Waiting 30 seconds before retry..." >>"$LOG"
      sleep 30
    else
      {
        echo "$(date): All retry attempts failed during header renaming. Check input files and system resources."
        echo "PREPARE_READS_ERROR: step=header_renaming exit_code=1"
      } >>"$LOG"
      exit 1
    fi
  fi
done

