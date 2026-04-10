#!/usr/bin/env bash
set -euo pipefail

project="$1"
sample="$2"
pythonhashseed="$3"
in_r1="$4"
in_r2="$5"
out_r1="$6"
out_r2="$7"
cleaned_token="$8"
raw_fastqc_html1="$9"
raw_fastqc_html2="${10}"
cleaned_fastqc_html1="${11}"
cleaned_fastqc_html2="${12}"
adapters="${13}"
bbduk_params="${14}"
fastqc_threads="${15}"
temp_dir="${16}"
log_path="${17}"

LOG="$log_path"

# Ensure output and log directories exist
mkdir -p "$(dirname "$out_r1")"
mkdir -p "$(dirname "$raw_fastqc_html1")"
mkdir -p "$(dirname "$LOG")"

# Use a per-job temp directory under the configured base temp root.
# This avoids collisions when Snakemake runs jobs in parallel.
if [[ "$temp_dir" = /* ]]; then
  _base_tmp_root="$temp_dir"
else
  _base_tmp_root="$(pwd)/${temp_dir}"
fi

_job_tmp_parent="${_base_tmp_root}/jobs/clean_reads/${project}/${sample}"
mkdir -p "$_job_tmp_parent"
job_tmp="$(mktemp -d "${_job_tmp_parent}/run_XXXXXXXX")"

_cleanup_job_tmp() {
  local exit_code="$?"
  if [[ "$exit_code" -eq 0 ]]; then
    rm -rf "$job_tmp" || true
    # Prune empty tmp parents so no tmp/ is left behind on success.
    rmdir "$_job_tmp_parent" 2>/dev/null || true
    rmdir "$(dirname "$_job_tmp_parent")" 2>/dev/null || true
    rmdir "$(dirname "$(dirname "$_job_tmp_parent")")" 2>/dev/null || true
    rmdir "$(dirname "$(dirname "$(dirname "$_job_tmp_parent")")")" 2>/dev/null || true
  else
    echo "$(date): CLEAN_READS_TMPDIR_RETAINED=${job_tmp}" >>"$LOG" || true
  fi
}
trap _cleanup_job_tmp EXIT

{
  echo "$(date): CLEAN_READS start for project=${project}, sample=${sample}"
  echo "$(date): Input R1=${in_r1}"
  echo "$(date): Input R2=${in_r2}"
  echo "$(date): Output R1=${out_r1}"
  echo "$(date): Output R2=${out_r2}"
  echo "$(date): Raw FastQC outdir=$(dirname "$raw_fastqc_html1")"
  echo "$(date): Cleaned FastQC outdir=$(dirname "$cleaned_fastqc_html1")"
  echo "$(date): Adapters=${adapters}"
  echo "$(date): BBDuk params=${bbduk_params}"
  echo "$(date): FastQC threads=${fastqc_threads}"
  echo "$(date): TEMP base dir=${temp_dir}"
  echo "$(date): JOB TMP dir=${job_tmp}"
} >> "$LOG"

export TMPDIR="$job_tmp"
export TMP="$job_tmp"
export TEMP="$job_tmp"
export PYTHONHASHSEED="$pythonhashseed"

raw_qc_dir="$(dirname "$raw_fastqc_html1")"
clean_qc_dir="$(dirname "$cleaned_fastqc_html1")"

_fastqc_out_stem() {
  # FastQC typically strips .fastq/.fq and (if present) only strips the trailing .gz.
  # Compute a conservative expected output stem for the input filename.
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

{
  echo "$(date): CLEAN_READS running fastqc (raw inputs)"
  echo "COMMAND: fastqc -t ${fastqc_threads} ${in_r1} ${in_r2} -o ${raw_qc_dir}"
} >>"$LOG"

if ! fastqc -t "$fastqc_threads" "$in_r1" "$in_r2" -o "$raw_qc_dir" --quiet >>"$LOG" 2>&1; then
  {
    echo "$(date): CLEAN_READS_ERROR: fastqc failed (raw inputs)"
    echo "CLEAN_READS_ERROR: step=fastqc_raw exit_code=$?"
  } >>"$LOG"
  exit 1
fi

# Normalize FastQC output filenames to the rule-declared targets
_move_fastqc_outputs "$in_r1" "$raw_fastqc_html1"
_move_fastqc_outputs "$in_r2" "$raw_fastqc_html2"

# Clean reads with BBDuk (log command; rely on timestamps for duration)
{
  echo "$(date): CLEAN_READS running bbduk.sh"
  echo "COMMAND: bbduk.sh in=${in_r1} in2=${in_r2} out=${out_r1} out2=${out_r2} ref=${adapters} ${bbduk_params}"
} >>"$LOG"

if ! bbduk.sh in="$in_r1" in2="$in_r2" out="$out_r1" out2="$out_r2" ref="$adapters" $bbduk_params >>"$LOG" 2>&1; then
  {
    echo "$(date): CLEAN_READS_ERROR: bbduk.sh failed (see above for details)"
    echo "CLEAN_READS_ERROR: step=bbduk exit_code=$?"
  } >>"$LOG"
  exit 1
fi

# Normalize read headers
if ! sed -Ei 's#^(@[A-Za-z0-9]+)/[12]$#\1#' "$out_r1" "$out_r2" >>"$LOG" 2>&1; then
  {
    echo "$(date): CLEAN_READS_ERROR: sed header normalization failed"
    echo "CLEAN_READS_ERROR: step=sed exit_code=$?"
  } >>"$LOG"
  exit 1
fi

# Run FastQC on cleaned reads
{
  echo "$(date): CLEAN_READS running fastqc"
  echo "COMMAND: fastqc -t ${fastqc_threads} ${out_r1} ${out_r2} -o ${clean_qc_dir}"
} >>"$LOG"

if ! fastqc -t "$fastqc_threads" "$out_r1" "$out_r2" -o "$clean_qc_dir" --quiet >>"$LOG" 2>&1; then
  {
    echo "$(date): CLEAN_READS_ERROR: fastqc failed"
    echo "CLEAN_READS_ERROR: step=fastqc exit_code=$?"
  } >>"$LOG"
  exit 1
fi

_move_fastqc_outputs "$out_r1" "$cleaned_fastqc_html1"
_move_fastqc_outputs "$out_r2" "$cleaned_fastqc_html2"

touch "$cleaned_token"

{
  echo "$(date): CLEAN_READS completed successfully"
  echo "CLEAN_READS_SUCCESS: exit_code=0"
} >>"$LOG"

