#!/usr/bin/env bash


# 11 Tajima's D (sliding windows)
# sliding_window_tajimasD.sh
# Make sliding windows and compute Tajima's D per-window by extracting regions via bcftools

set -euo pipefail
IFS=$'\n\t'

VCF="Path to vcf file"
POPMAP="Path to popmap.txt file"
OUTDIR="./tajima_sliding_results"
WINDOW_BP="${WINDOW_BP:-1432000}"
STEP_BP="${STEP_BP:-716000}"  
TMPDIR="${OUTDIR}/tmp"
VERBOSE=1

BCFTOOLS="$(command -v bcftools || true)"
TABIX="$(command -v tabix || true)"
VCFTOOLS="$(command -v vcftools || true)"
BEDTOOLS="$(command -v bedtools || true)"

log(){ if [[ "$VERBOSE" -eq 1 ]]; then printf '[%s] %s\n' "$(date '+%F %T')" "$*"; else printf '%s\n' "$*"; fi }

for cmd in "$BCFTOOLS" "$TABIX" "$VCFTOOLS" "$BEDTOOLS"; do
  if [[ -z "$cmd" ]]; then log "ERROR: required tool missing (bcftools, tabix, vcftools, bedtools)"; exit 1; fi
done

mkdir -p "$OUTDIR" "$TMPDIR" "${OUTDIR}/logs"

log "VCF: $VCF"
log "WINDOW_BP: $WINDOW_BP  STEP_BP: $STEP_BP"
log "OUTDIR: $OUTDIR"

# 1) produce contigs file
CONTIGS="${OUTDIR}/contigs.tsv"
$BCFTOOLS view -h "$VCF" | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*(length=([0-9]+)).*/\1\t\3/' > "${CONTIGS}" || true
if [[ ! -s "${CONTIGS}" ]]; then
  log "WARN: contigs not found in VCF header; attempting to infer contigs (may be slower)."
  # fallback: list unique contigs from VCF (no lengths)
  $BCFTOOLS query -f '%CHROM\n' "$VCF" | sort | uniq | awk '{print $1"\t0"}' > "${CONTIGS}"
fi

# 2) make sliding windows across contigs (bedtools makewindows)
WIN_BED="${OUTDIR}/windows.w${WINDOW_BP}.s${STEP_BP}.bed"
log "Generating sliding windows -> $WIN_BED ..."
# bedtools makewindows expects a genome file with contig lengths; our CONTIGS may have 0 for lengths if missing.
# If contig lengths are zero (fallback), bedtools will fail; in that case create windows by scanning positions (not ideal).
if awk '{if($2==0){exit 1}}' "${CONTIGS}" 2>/dev/null; then
  # If any length is zero this command returns 1
  log "Some contig lengths unknown (0). Building windows per-contig from VCF positions (fallback)..."
  rm -f "$WIN_BED"
  while read -r ct len; do
    if [[ "$len" -gt 0 ]]; then
      # use bedtools for this contig
      echo -n "" >> "$WIN_BED"
      bedtools makewindows -g <(echo -e "${ct}\t${len}") -w "$WINDOW_BP" -s "$STEP_BP" >> "$WIN_BED"
    else
      # fallback: get max position from VCF for this contig - FIXED: avoid using END as variable name
      maxpos=$($BCFTOOLS query -r "${ct}" -f '%POS\n' "$VCF" 2>/dev/null | sort -n | tail -1)
      if [[ -z "$maxpos" || "$maxpos" == "0" ]]; then
        log "  No variants for contig ${ct} in VCF; skipping windows for this contig."
        continue
      fi
      # create windows for contig using awk arithmetic
      start=0
      while [[ $start -lt $maxpos ]]; do
        end=$((start + WINDOW_BP - 1))
        if [[ $end -gt $maxpos ]]; then end=$maxpos; fi
        printf "%s\t%s\t%s\n" "$ct" "$start" "$end" >> "$WIN_BED"
        start=$((start + STEP_BP))
      done
    fi
  done < "${CONTIGS}"
else
  # normal case: contig lengths present and >0
  bedtools makewindows -g "${CONTIGS}" -w "$WINDOW_BP" -s "$STEP_BP" > "$WIN_BED"
fi

# 3) prep poplists
awk 'NF>=2{print $1"\t"$2}' "$POPMAP" > "${OUTDIR}/popmap.cleaned.txt"
mapfile -t POPS < <(awk 'NF>=2{print $2}' "${OUTDIR}/popmap.cleaned.txt" | awk '!seen[$0]++')

# 4) iterate windows and compute TajimaD for each window
OUT_MASTER="${OUTDIR}/tajima_sliding_allsets_combined.csv"
echo "chrom,start,end,n_variants,tajimaD,set" > "$OUT_MASTER"

i=0
# Read the BED (tab-delimited); bedtools uses 0-based starts; convert to 1-based for bcftools
while IFS=$'\t' read -r CHR START END_POS; do 
  i=$((i+1))
  # convert 0-based START to 1-based START1 for bcftools
  START1=$((START + 1))
  if [[ $START1 -lt 1 ]]; then START1=1; fi
  region="${CHR}:${START1}-${END_POS}"
  region_prefix="${OUTDIR}/window_${i}_${CHR}_${START}_${END_POS}"
  log "[$i] region ${region} (orig 0-based ${CHR}:${START}-${END_POS})"

  # extract region (fast) - bcftools returns non-zero if no records; capture rc
  set +e
  $BCFTOOLS view -r "${CHR}:${START1}-${END_POS}" -O z -o "${region_prefix}.vcf.gz" "$VCF" 2> "${OUTDIR}/logs/window_${i}.bcftools.stderr"
  rc_bcf=$?
  set -e

  if [[ $rc_bcf -ne 0 || ! -s "${region_prefix}.vcf.gz" ]]; then
    # no variants in region or bcftools failed to extract — produce NA rows for ALL and pops
    echo "${CHR},${START1},${END_POS},0,NA,ALL" >> "$OUT_MASTER"
    for p in "${POPS[@]}"; do
      echo "${CHR},${START1},${END_POS},0,NA,${p}" >> "$OUT_MASTER"
    done
    # cleanup any partial files
    rm -f "${region_prefix}.vcf.gz" "${region_prefix}.vcf.gz.tbi" || true
    continue
  fi

  # create index if missing
  set +e
  $TABIX -f -p vcf "${region_prefix}.vcf.gz" 2>/dev/null || true
  set -e

  # region length for vcftools Tajima (pass window size equal to region length)
  reg_len=$((END_POS - START1 + 1))
  if [[ $reg_len -le 0 ]]; then
    # weird edge case, skip
    echo "${CHR},${START1},${END_POS},0,NA,ALL" >> "$OUT_MASTER"
    for p in "${POPS[@]}"; do
      echo "${CHR},${START1},${END_POS},0,NA,${p}" >> "$OUT_MASTER"
    done
    rm -f "${region_prefix}.vcf.gz" "${region_prefix}.vcf.gz.tbi" || true
    continue
  fi

  # compute TajimaD for ALL in this region
  set +e
  $VCFTOOLS --gzvcf "${region_prefix}.vcf.gz" --TajimaD "${reg_len}" --out "${region_prefix}.tajima" > "${OUTDIR}/logs/window_${i}.vcftools.stdout" 2> "${OUTDIR}/logs/window_${i}.vcftools.stderr"
  rc_all=$?
  set -e

  if [[ $rc_all -eq 0 && -s "${region_prefix}.tajima.Tajima.D" ]]; then
    awk -v CHR="$CHR" -v START1="$START1" -v END_POS="$END_POS" -v S="ALL" '
    BEGIN {OFS=","}
    NR>1 {
      n_variants = $3
      tajimaD = $4

      # Handle missing/invalid tajimaD values
      if (tajimaD == "nan" || tajimaD == "NaN" || tajimaD == "-nan") {
        tajimaD = "NA"
      }

      print CHR, START1, END_POS, n_variants, tajimaD, S
    }' "${region_prefix}.tajima.Tajima.D" >> "$OUT_MASTER"
  else
    echo "${CHR},${START1},${END_POS},0,NA,ALL" >> "$OUT_MASTER"
  fi

  # per-population
  for p in "${POPS[@]}"; do
    poplist="${OUTDIR}/${p}.samples.txt"
    # if a poplist doesn't exist, make it
    if [[ ! -f "$poplist" ]]; then awk -v P="$p" '$2==P{print $1}' "${OUTDIR}/popmap.cleaned.txt" > "$poplist"; fi

    set +e
    $VCFTOOLS --gzvcf "${region_prefix}.vcf.gz" --keep "$poplist" --TajimaD "${reg_len}" --out "${region_prefix}.${p}.tajima" > /dev/null 2> "${OUTDIR}/logs/window_${i}.${p}.vcftools.stderr"
    rcp=$?
    set -e

    if [[ $rcp -eq 0 && -s "${region_prefix}.${p}.tajima.Tajima.D" ]]; then
      # FIXED: Properly parse vcftools output for populations
      awk -v CHR="$CHR" -v START1="$START1" -v END_POS="$END_POS" -v S="$p" '
      BEGIN {OFS=","}
      NR>1 {
        n_variants = $3
        tajimaD = $4

        # Handle missing/invalid tajimaD values
        if (tajimaD == "nan" || tajimaD == "NaN" || tajimaD == "-nan") {
          tajimaD = "NA"
        }

        print CHR, START1, END_POS, n_variants, tajimaD, S
      }' "${region_prefix}.${p}.tajima.Tajima.D" >> "$OUT_MASTER"
    else
      echo "${CHR},${START1},${END_POS},0,NA,${p}" >> "$OUT_MASTER"
    fi
  done

  # cleanup region files to save space
  rm -f "${region_prefix}.vcf.gz" "${region_prefix}.vcf.gz.tbi" "${region_prefix}.tajima.Tajima.D" "${region_prefix}."*".Tajima.D" || true

done < "$WIN_BED"

log "Sliding TajimaD completed. Combined results: $OUT_MASTER"
log "Notes: each row = chrom,start,end,n_variants,tajimaD,set (set is ALL or population code)."
