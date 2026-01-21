#!/usr/bin/env bash


#1 Data Curation
set -euo pipefail

# merge_and_trim_fastp.sh
# Workflow A: merge duplicate FASTQs (raw), prepare work dir, run fastp trimming for each sample


THREADS=${THREADS:-4}     
RAW_DIR="$(pwd)"
MERGED_DIR="${RAW_DIR}/merged_fastq"
WORK_DIR="${RAW_DIR}/work_fastq"
TRIM_DIR="${RAW_DIR}/trimmed"
REPORT_DIR="${RAW_DIR}/fastp_reports"
LOG_DIR="${RAW_DIR}/logs"

mkdir -p "$MERGED_DIR" "$WORK_DIR" "$TRIM_DIR" "$REPORT_DIR" "$LOG_DIR"


DUPLICATES=(
  "WR_C05 WR_C10 WR_05_10"
  "WR_C7  WR_13  WR_07_13"
)


echo "RAW_DIR = $RAW_DIR"
echo "MERGED_DIR = $MERGED_DIR"
echo "WORK_DIR = $WORK_DIR"
echo "TRIM_DIR = $TRIM_DIR"
echo "REPORT_DIR = $REPORT_DIR"
echo "LOG_DIR = $LOG_DIR"
echo

# 1) Merge duplicate FASTQs (R1 and R2) using zcat -> gzip
for entry in "${DUPLICATES[@]}"; do
  read -r A B MERGED <<<"$entry"

  A_R1="${RAW_DIR}/${A}_R1.fastq.gz"
  A_R2="${RAW_DIR}/${A}_R2.fastq.gz"
  B_R1="${RAW_DIR}/${B}_R1.fastq.gz"
  B_R2="${RAW_DIR}/${B}_R2.fastq.gz"

  OUT_R1="${MERGED_DIR}/${MERGED}_R1.fastq.gz"
  OUT_R2="${MERGED_DIR}/${MERGED}_R2.fastq.gz"

  echo "Merging duplicates: $A + $B -> $MERGED"
  if [[ ! -f "$A_R1" || ! -f "$A_R2" || ! -f "$B_R1" || ! -f "$B_R2" ]]; then
    echo "ERROR: one of the files for $A or $B is missing. Check filenames." >&2
    echo "  Missing: " \
      && [[ ! -f "$A_R1" ]] && echo "  $A_R1" \
      && [[ ! -f "$A_R2" ]] && echo "  $A_R2" \
      && [[ ! -f "$B_R1" ]] && echo "  $B_R1" \
      && [[ ! -f "$B_R2" ]] && echo "  $B_R2"
    exit 1
  fi

  # merge R1
  echo "  -> merging R1 to $OUT_R1"
  zcat "$A_R1" "$B_R1" | gzip > "$OUT_R1"
  # merge R2
  echo "  -> merging R2 to $OUT_R2"
  zcat "$A_R2" "$B_R2" | gzip > "$OUT_R2"

  if [[ ! -s "$OUT_R1" || ! -s "$OUT_R2" ]]; then
    echo "ERROR: merged outputs appear empty for $MERGED" >&2
    exit 1
  fi
done

# Skip EV_24
SKIP_SAMPLE="EV_24"

# Symlinks to keep track of merged sample names
MERGED_NAMES=()
for entry in "${DUPLICATES[@]}"; do
  read -r A B MERGED <<<"$entry"
  MERGED_NAMES+=("$MERGED")
done

# Create symlinks: for merged samples, point to merged_fastq; for others, symlink original raw files
echo
echo "Preparing work_fastq/ (symlinks to merged or original FASTQs). Skipping $SKIP_SAMPLE."
for r1 in "${RAW_DIR}"/*_R1.fastq.gz; do
  sample=$(basename "$r1" _R1.fastq.gz)

  # skip EV_24 (high-missingness sample)
  if [[ "$sample" == "$SKIP_SAMPLE" ]]; then
    echo "  Skipping $sample"
    continue
  fi

 
  skip=false
  for entry in "${DUPLICATES[@]}"; do
    read -r A B MERGED <<<"$entry"
    if [[ "$sample" == "$A" || "$sample" == "$B" ]]; then
      skip=true
      break
    fi
  done
  if $skip; then
    echo "  Skipping original duplicate member $sample (merged version will be used)"
    continue
  fi

  # Skip if name already present
  is_merged=false
  for m in "${MERGED_NAMES[@]}"; do
    if [[ "$sample" == "$m" ]]; then
      is_merged=true
      break
    fi
  done
  if $is_merged; then
    # create symlink to merged file
    src_r1="${MERGED_DIR}/${sample}_R1.fastq.gz"
    src_r2="${MERGED_DIR}/${sample}_R2.fastq.gz"
    if [[ ! -f "$src_r1" || ! -f "$src_r2" ]]; then
      echo "ERROR: expected merged files for $sample not found in $MERGED_DIR" >&2
      exit 1
    fi
    ln -sf "$src_r1" "${WORK_DIR}/${sample}_R1.fastq.gz"
    ln -sf "$src_r2" "${WORK_DIR}/${sample}_R2.fastq.gz"
    echo "  Linked merged sample $sample -> work_fastq/"
    continue
  fi

  # otherwise, symlink original raw files
  src_r1="${RAW_DIR}/${sample}_R1.fastq.gz"
  src_r2="${RAW_DIR}/${sample}_R2.fastq.gz"
  if [[ ! -f "$src_r1" || ! -f "$src_r2" ]]; then
    echo "ERROR: missing pair for $sample: $src_r1 or $src_r2" >&2
    exit 1
  fi
  ln -sf "$src_r1" "${WORK_DIR}/${sample}_R1.fastq.gz"
  ln -sf "$src_r2" "${WORK_DIR}/${sample}_R2.fastq.gz"
  echo "  Linked $sample -> work_fastq/"
done


for entry in "${DUPLICATES[@]}"; do
  read -r A B MERGED <<<"$entry"
  src_r1="${MERGED_DIR}/${MERGED}_R1.fastq.gz"
  src_r2="${MERGED_DIR}/${MERGED}_R2.fastq.gz"
  if [[ ! -f "$src_r1" || ! -f "$src_r2" ]]; then
    echo "ERROR: merged files expected but not found for $MERGED" >&2
    exit 1
  fi
  ln -sf "$src_r1" "${WORK_DIR}/${MERGED}_R1.fastq.gz"
  ln -sf "$src_r2" "${WORK_DIR}/${MERGED}_R2.fastq.gz"
  echo "  Linked merged sample ${MERGED} -> work_fastq/"
done

echo
echo "Contents of work_fastq/:"
ls -1 "$WORK_DIR" | sed -n '1,200p'
echo

# 3) Run fastp for each sample in work_fastq/
echo "Running fastp trimming on samples in work_fastq/ using $THREADS threads per sample."
for r1 in "${WORK_DIR}"/*_R1.fastq.gz; do
  sample=$(basename "$r1" _R1.fastq.gz)
  r2="${WORK_DIR}/${sample}_R2.fastq.gz"
  if [[ ! -f "$r2" ]]; then
    echo "ERROR: read2 missing for $sample ($r2)" >&2
    exit 1
  fi

  out_r1="${TRIM_DIR}/${sample}.trim.R1.fastq.gz"
  out_r2="${TRIM_DIR}/${sample}.trim.R2.fastq.gz"
  html_report="${REPORT_DIR}/${sample}.fastp.html"
  json_report="${REPORT_DIR}/${sample}.fastp.json"
  log_file="${LOG_DIR}/${sample}_fastp.log"

  echo "Trimming $sample -> $out_r1 , $out_r2"
  fastp -i "$r1" -I "$r2" \
    -o "$out_r1" -O "$out_r2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --cut_mean_quality 20 \
    --cut_window_size 4 \
    --length_required 50 \
    -h "$html_report" -j "$json_report" -w "$THREADS" \
    2>&1 | tee "$log_file"

  # sanity check
  if [[ ! -s "$out_r1" || ! -s "$out_r2" ]]; then
    echo "ERROR: fastp produced empty outputs for $sample. See $log_file" >&2
    exit 1
  fi
done

echo
echo "All samples trimmed. Trimmed files are in: $TRIM_DIR"
echo "fastp reports in: $REPORT_DIR"
echo "Logs in: $LOG_DIR"
echo "Done at: $(date)"
