#!/usr/bin/env bash


#2 Mapping to a reference and mapping summary

set -euo pipefail

# map_and_stats_bwa.sh
# Maps paired-end trimmed FASTQs to a reference using bwa mem -> samtools, produces sorted BAMs and stats.

REF="/home/manas/Downloads/Vultures_popGen_analysis/results_modular/gatk_pipeline_20250924_2254/ref.fa"

THREADS="${THREADS:-16}"

# output directories
OUTDIR="../mapped_bam"          
LOGDIR="../mapping_logs"
STATDIR="../mapping_stats"

BWA_EXTRA="-M -K 100000000"

cd "$(dirname "$0")" || exit 1  # change to trimmed/ directory where script is located
mkdir -p "$OUTDIR" "$LOGDIR" "$STATDIR"

if [[ ! -f "$REF" ]]; then
  echo "ERROR: reference FASTA not found at: $REF" >&2
  exit 1
fi

echo "Using REF = $REF"
echo "Threads = $THREADS"
echo "Output BAM dir = $OUTDIR"
echo "Logs dir = $LOGDIR"
echo "Stats dir = $STATDIR"
echo

# check bwa index files; if missing, run bwa index
if [[ ! -f "${REF}.bwt" ]]; then
  echo "bwa index not found for $REF -> running bwa index"
  bwa index "$REF"
  echo "bwa index completed."
else
  echo "bwa index found for reference."
fi
echo

# helper function to extract numbers (handles formats like "12345 + 0")

num_only() {
  echo "$1" | awk '{print $1}'
}


# mapping loop

R1_LIST=( *.trim.R1.fastq.gz )

if [[ ${#R1_LIST[@]} -eq 0 ]]; then
  echo "No *.trim.R1.fastq.gz found in $(pwd). Aborting." >&2
  exit 1
fi

SUMMARY_FILE="${STATDIR}/mapping_summary.tsv"
echo -e "sample\ttotal_reads\tmapped_reads\tpercent_mapped\tproperly_paired\tmapped_pairs" > "$SUMMARY_FILE"

for R1 in "${R1_LIST[@]}"; do
  # derive sample name: remove .trim.R1.fastq.gz
  SAMPLE="${R1%%.trim.R1.fastq.gz}"
  R2="${SAMPLE}.trim.R2.fastq.gz"

  if [[ ! -f "$R2" ]]; then
    echo "WARNING: missing R2 for sample $SAMPLE -> skipping" >&2
    continue
  fi

  BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
  BAM_IDX="${BAM}.bai"
  LOG="${LOGDIR}/${SAMPLE}.bwa_samtools.log"
  FLAGSTAT="${STATDIR}/${SAMPLE}.flagstat.txt"
  IDXSTATS="${STATDIR}/${SAMPLE}.idxstats.txt"
  SAMTOOLS_STATS="${STATDIR}/${SAMPLE}.samtools.stats.txt"

  echo "Processing sample: $SAMPLE"
  echo "R1: $R1"
  echo "R2: $R2"
  echo "BAM: $BAM"
  echo "Log: $LOG"
  echo

    # Create a unique @RG ID and PU (using sample name and timestamp)
  RG_ID="${SAMPLE}.1"
  RG_PU="${SAMPLE}.1"

  RG_STR="@RG\\tID:${RG_ID}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tLB:lib1\\tPU:${RG_PU}"


  echo "[$(date +'%F %T')] bwa mem -> samtools sort (threads=$THREADS) for sample $SAMPLE" | tee "$LOG"
  set -x
  bwa mem $BWA_EXTRA -t "$THREADS" -R "$RG_STR" "$REF" "$R1" "$R2" \
    2>> "$LOG" | \
    samtools view -b -@ 4 -F 4 - 2>> "$LOG" | \
    samtools sort -@ "$THREADS" -o "$BAM" - 2>> "$LOG"
  set +x


  # Index BAM
  echo "Indexing BAM..." | tee -a "$LOG"
  samtools index "$BAM" 2>> "$LOG"

  # Generate stats
  echo "samtools flagstat -> $FLAGSTAT" | tee -a "$LOG"
  samtools flagstat "$BAM" > "$FLAGSTAT" 2>> "$LOG"

  echo "samtools idxstats -> $IDXSTATS" | tee -a "$LOG"
  samtools idxstats "$BAM" > "$IDXSTATS" 2>> "$LOG"

  echo "samtools stats -> $SAMTOOLS_STATS" | tee -a "$LOG"
  samtools stats "$BAM" > "$SAMTOOLS_STATS" 2>> "$LOG"

  # Parse key numbers from flagstat
  total_reads=$(awk 'NR==1{print $1}' "$FLAGSTAT")
  mapped_reads=$(awk '/mapped \(/ {print $1; exit}' "$FLAGSTAT")
  properly_paired=$(awk '/properly paired/ {print $1; exit}' "$FLAGSTAT")

  # fallback safety (if awk misses)
  total_reads=${total_reads:-0}
  mapped_reads=${mapped_reads:-0}
  properly_paired=${properly_paired:-0}

  if [[ "$total_reads" -gt 0 ]]; then
    percent_mapped=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{printf("%.2f", (m/t)*100)}')
  else
    percent_mapped="0.00"
  fi

  # mapped pairs approximated as mapped_reads / 2 (paired-end) - integer division
  if [[ "$mapped_reads" -gt 0 ]]; then
    mapped_pairs=$(( mapped_reads / 2 ))
  else
    mapped_pairs=0
  fi

  # Append to summary
  echo -e "${SAMPLE}\t${total_reads}\t${mapped_reads}\t${percent_mapped}\t${properly_paired}\t${mapped_pairs}" >> "$SUMMARY_FILE"

  echo "Sample $SAMPLE done. Summary appended to $SUMMARY_FILE"
  echo
done

echo "All samples processed. Final mapping summary: $SUMMARY_FILE"
echo "Per-sample stats are in: $STATDIR, logs in: $LOGDIR, BAMs in: $OUTDIR"
echo "Done at: $(date)"
