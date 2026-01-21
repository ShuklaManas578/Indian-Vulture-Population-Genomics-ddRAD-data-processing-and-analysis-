#!/usr/bin/env bash


#3 Genotype Calling using GATK

set -euo pipefail

# run_GATK_parallel_genotyping.sh
# Produces per-sample GVCFs in parallel via GATK HaplotypeCaller, then joint-genotypes to raw_joint.vcf.gz

REF="${REF:-/home/manas/Downloads/Vultures_popGen_analysis/results_modular/gatk_pipeline_20250924_2254/ref.fa}"
BAM_DIR="${BAM_DIR:-$(pwd)}"                   
GVCF_DIR="${GVCF_DIR:-${BAM_DIR}/gvcfs}"
LOG_DIR="${LOG_DIR:-${BAM_DIR}/logs_hc}"
JOINT_DIR="${JOINT_DIR:-${BAM_DIR}/joint_vcf}"
OUT_JOINT_VCF="${OUT_JOINT_VCF:-${JOINT_DIR}/raw_joint.vcf.gz}"

GATK_TMPDIR_BASE="${GATK_TMPDIR_BASE:-${TMPDIR:-/tmp}/gatk_tmp_$(id -u)}"

# Defaults 
THREADS_TOTAL="${THREADS_TOTAL:-$(nproc --all 2>/dev/null || echo 16)}"
JOBS="${JOBS:-3}"     
JAVA_XMX_PREF="${JAVA_XMX_PREF:-24}"  

# detect available CPUs and total memory (KB)
CPU_ON_SYSTEM="$(nproc --all || echo 16)"
MEM_KB="$(awk '/MemTotal/ {print $2}' /proc/meminfo || echo 0)"
MEM_GB=$(( MEM_KB / 1024 / 1024 ))

# normalize THREADS_TOTAL
if [[ "$CPU_ON_SYSTEM" -lt "$THREADS_TOTAL" ]]; then
  echo "Detected only $CPU_ON_SYSTEM CPUs but THREADS_TOTAL=$THREADS_TOTAL requested."
  THREADS_TOTAL="$CPU_ON_SYSTEM"
fi

# ensure JOBS sane
if [[ "$JOBS" -lt 1 ]]; then JOBS=1; fi
if [[ "$JOBS" -gt "$THREADS_TOTAL" ]]; then
  echo "Reducing JOBS from $JOBS to THREADS_TOTAL ($THREADS_TOTAL) to avoid oversubscription."
  JOBS="$THREADS_TOTAL"
fi

HMM_THREADS=$(( THREADS_TOTAL / JOBS ))
if [[ "$HMM_THREADS" -lt 1 ]]; then HMM_THREADS=1; fi

# compute per-job memory Xmx (GB) logic
PER_JOB_MEM_GB=$(( MEM_GB / JOBS ))
if [[ $PER_JOB_MEM_GB -ge $JAVA_XMX_PREF ]]; then
  XMX_G="${JAVA_XMX_PREF}"
elif [[ $PER_JOB_MEM_GB -ge 12 ]]; then
  XMX_G=12
elif [[ $PER_JOB_MEM_GB -ge 8 ]]; then
  XMX_G=8
else
  XMX_G=4
fi

echo "System cpus detected: ${CPU_ON_SYSTEM}"
echo "System RAM (GB): ${MEM_GB}"
echo "THREADS_TOTAL set to: ${THREADS_TOTAL}"
echo "Concurrent JOBS set to: ${JOBS}"
echo "Per-job HMM threads (--native-pair-hmm-threads): ${HMM_THREADS}"
echo "Per-job Java Xmx selected: ${XMX_G}g"
echo "GATK tmpdir base: ${GATK_TMPDIR_BASE}"
echo


command -v gatk >/dev/null 2>&1 || { echo "ERROR: gatk not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH"; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "WARNING: tabix not found - indexing may fail"; }

mkdir -p "$GVCF_DIR" "$LOG_DIR" "$JOINT_DIR" "$GATK_TMPDIR_BASE"

# check reference exists and index (.fai) exists (create if missing)
if [[ ! -f "$REF" ]]; then
  echo "ERROR: reference not found at $REF" >&2
  exit 1
fi
if [[ ! -f "${REF}.fai" ]]; then
  echo "Indexing reference fasta (samtools faidx) ..."
  samtools faidx "$REF"
fi
if [[ ! -f "${REF%.fa}.dict" && ! -f "${REF%.fasta}.dict" ]]; then
  echo "Creating reference dict (gatk CreateSequenceDictionary) ..."
  gatk CreateSequenceDictionary -R "$REF" -O "${REF%.fa}.dict" || true
fi

# Build sample list
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name '*.sorted.bam' -print | sort || true)
if [[ "${#BAM_FILES[@]}" -eq 0 ]]; then
  echo "No *.sorted.bam files found in ${BAM_DIR}. Aborting." >&2
  exit 1
fi

# helper function to create per-sample GVCF command
make_hc_cmd() {
  local bam="$1"
  local sample
  sample=$(basename "$bam" .sorted.bam)
  local out_gvcf="${GVCF_DIR}/${sample}.g.vcf.gz"
  local log="${LOG_DIR}/${sample}.hc.log"
  local tmpdir="${GATK_TMPDIR_BASE}/tmp_${sample}"

  
  printf '%s\n' "mkdir -p '${tmpdir}' && chmod 700 '${tmpdir}' && \
gatk --java-options \"-Xmx${XMX_G}g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true\" HaplotypeCaller \
-R \"${REF}\" \
-I \"${bam}\" \
-O \"${out_gvcf}\" \
-ERC GVCF \
--native-pair-hmm-threads ${HMM_THREADS} \
--tmp-dir \"${tmpdir}\" \
> \"${log}\" 2>&1"
}

# Create a commands file with one command per line for parallel execution 
CMDFILE="${LOG_DIR}/hc_commands.txt"
: > "$CMDFILE"

for bam in "${BAM_FILES[@]}"; do
  sample=$(basename "$bam" .sorted.bam)
  out_gvcf="${GVCF_DIR}/${sample}.g.vcf.gz"
  # skip if already exists
  if [[ -f "$out_gvcf" ]]; then
    echo "Skipping existing gVCF for ${sample}"
    continue
  fi
  make_hc_cmd "$bam" >> "$CMDFILE"
done

NUM_JOBS=$(wc -l < "$CMDFILE" || echo 0)
if [[ "$NUM_JOBS" -eq 0 ]]; then
  echo "No new HaplotypeCaller jobs to run."
else
  echo "Will run ${NUM_JOBS} HaplotypeCaller jobs using up to ${JOBS} concurrent jobs"
  echo


  if command -v parallel >/dev/null 2>&1; then
    echo "Using GNU parallel to run HaplotypeCaller jobs..."
    # run parallel with JOBS concurrency, read commands from CMDFILE
    parallel --joblog "${LOG_DIR}/parallel_joblog.txt" -j "${JOBS}" --bar < "$CMDFILE"
  else
    echo "GNU parallel not found; using xargs fallback"

    # Using --max-procs (or -P) to control concurrency
    cat "$CMDFILE" | xargs -I CMD -P "${JOBS}" /bin/bash -lc CMD
  fi
fi

echo
echo "Verifying generated GVCFs..."
GVCF_LIST=()
for bam in "${BAM_FILES[@]}"; do
  sample=$(basename "$bam" .sorted.bam)
  g="${GVCF_DIR}/${sample}.g.vcf.gz"
  if [[ -f "$g" ]]; then
    GVCF_LIST+=("$g")
    # index if missing
    if [[ ! -f "${g}.tbi" ]]; then
      if command -v tabix >/dev/null 2>&1; then
        echo "Indexing ${g}..."
        tabix -p vcf "$g" || echo "Warning: tabix failed for ${g}"
      else
        echo "tabix not available to index ${g}"
      fi
    fi
  else
    echo "Warning: expected gVCF not found for ${sample}: ${g}"
  fi
done

if [[ "${#GVCF_LIST[@]}" -eq 0 ]]; then
  echo "ERROR: no gVCFs found. Aborting joint genotyping." >&2
  exit 1
fi

# Combine or import and genotype 
echo
echo "Combining GVCFs and running joint genotyping..."

COMBINED_GVCF="${JOINT_DIR}/cohort.combined.g.vcf.gz"

# build CombineGVCFs
COMBINE_CMD=(gatk --java-options "-Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" CombineGVCFs -R "${REF}" -O "${COMBINED_GVCF}")
for g in "${GVCF_LIST[@]}"; do
  COMBINE_CMD+=(--variant "$g")
done

echo "Running CombineGVCFs ..."
"${COMBINE_CMD[@]}"

# index combined gvcf
if [[ -f "${COMBINED_GVCF}" ]]; then
  if command -v tabix >/dev/null 2>&1; then
    tabix -p vcf "${COMBINED_GVCF}" || echo "Warning: tabix failed for ${COMBINED_GVCF}"
  fi
else
  echo "ERROR: CombineGVCFs failed to produce ${COMBINED_GVCF}" >&2
  exit 1
fi

# GenotypeGVCFs to produce joint raw VCF
echo "Running GenotypeGVCFs to produce joint VCF..."
gatk --java-options "-Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
  -R "${REF}" \
  -V "${COMBINED_GVCF}" \
  -O "${OUT_JOINT_VCF}"

# index final joint VCF
if [[ -f "${OUT_JOINT_VCF}" ]]; then
  if command -v tabix >/dev/null 2>&1; then
    tabix -p vcf "${OUT_JOINT_VCF}" || echo "Warning: tabix failed to index ${OUT_JOINT_VCF}"
  fi
  echo "Joint VCF written to: ${OUT_JOINT_VCF}"
else
  echo "ERROR: GenotypeGVCFs failed to create ${OUT_JOINT_VCF}" >&2
  exit 1
fi

if command -v bcftools >/dev/null 2>&1; then
  echo "Checking that PL fields are present in joint VCF (sample check)..."
  # sample a single position from joint VCF header
  first_pos=$(bcftools query -f '%CHROM\t%POS\n' -H "${OUT_JOINT_VCF}" | head -n1 || true)
  if [[ -n "$first_pos" ]]; then
    chrom=$(awk '{print $1}' <<<"$first_pos")
    pos=$(awk '{print $2}' <<<"$first_pos")
    
    bcftools query -f '%GT[\t%PL]\n' -r "${chrom}:${pos}-${pos}" "${OUT_JOINT_VCF}" >/dev/null 2>&1 && echo "PL fields appear present (bcftools query succeeded)."
  fi
fi

echo "Done. Logs are in: ${LOG_DIR}"
echo "Per-sample GVCFs are in: ${GVCF_DIR}"
echo "Joint VCF is: ${OUT_JOINT_VCF}"
