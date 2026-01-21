#!/usr/bin/env bash

#9 Nucleotide diversity (pi) 
# Compute windowed and site nucleotide diversity (pi) per population using vcftools

set -euo pipefail
IFS=$'\n\t'

VCF_DEFAULT="filtered_final_snps.vcf.reheadered_uncompressed.uniqueIDs.vcf.gz"

VCF="${1:-$VCF_DEFAULT}"
OUTDIR="${2:-pi_results}"
TARGET_SNPS="${3:-10}"   # target number of SNPs per window (default 10)
MIN_WINDOW_BP=5000       # lower bound for window (5 kb)
MAX_WINDOW_BP=5000000    # upper bound (5 Mb)

mkdir -p "$OUTDIR"
echo "[INFO] VCF: $VCF"
echo "[INFO] OUTDIR: $OUTDIR"
echo "[INFO] TARGET SNPs/window: $TARGET_SNPS"
echo

# Check VCF exists
if [[ ! -s "$VCF" ]]; then
  echo "[ERROR] VCF not found or empty: $VCF" >&2
  exit 1
fi

# 1) Total variants: try using .bim if present, otherwise count from VCF index
if [[ -f "renamed_data_plink.bim" ]]; then
  TOTAL_VARIANTS=$(wc -l < renamed_data_plink.bim)
  echo "[INFO] total variant sites from renamed_data_plink.bim: $TOTAL_VARIANTS"
else

  TOTAL_VARIANTS=$(zcat "$VCF" | grep -v '^#' | wc -l)
  echo "[INFO] total variant sites counted from VCF: $TOTAL_VARIANTS"
fi

# 2) genome total length from VCF header (sum contig lengths)
GENOME_LENGTH_BP=0
# parse contig lengths if present in header
if zcat "$VCF" | grep -m1 '^##contig' >/dev/null 2>&1; then
  GENOME_LENGTH_BP=$(zcat "$VCF" | grep '^##contig' | sed -n 's/.*length=\([0-9]*\).*/\1/p' | awk '{s+=$1} END{print s}')
fi

if [[ -z "$GENOME_LENGTH_BP" || "$GENOME_LENGTH_BP" -le 0 ]]; then
  echo "[WARN] Could not parse contig lengths from VCF header. Trying per_contig_variant_counts.tsv..."
  if [[ -s "per_contig_variant_counts.tsv" ]]; then

    echo "[WARN] No genome length; average bp per SNP won't be computed; using conservative defaults."
    GENOME_LENGTH_BP=0
  fi
fi

# If per_contig_variant_counts.tsv is availble, total variants can be obtained as cross-check
if [[ -s "per_contig_variant_counts.tsv" ]]; then
  VCNT_SUM=$(awk '{s+=$2} END{print s}' per_contig_variant_counts.tsv)
  if [[ -n "$VCNT_SUM" && "$VCNT_SUM" -gt 0 ]]; then
    echo "[INFO] total variants from per_contig_variant_counts.tsv: $VCNT_SUM"
    if [[ "$VCNT_SUM" -lt "$TOTAL_VARIANTS" ]]; then
      # do not override if inconsistent; but we still print both
      :
    fi
  fi
fi

# compute avg bp per SNP
if [[ "$GENOME_LENGTH_BP" -gt 0 && "$TOTAL_VARIANTS" -gt 0 ]]; then
  AVG_BP_PER_SNP=$(awk -v g="$GENOME_LENGTH_BP" -v n="$TOTAL_VARIANTS" 'BEGIN{printf("%.2f", g/n)}')
  echo "[RESULT] genome total length from VCF header (sum contig lengths): $GENOME_LENGTH_BP"
  echo "[RESULT] total variant sites in VCF: $TOTAL_VARIANTS"
  echo "[RESULT] average bp per SNP: $AVG_BP_PER_SNP"
else
  echo "[WARN] Could not compute average bp per SNP; using fallback avg_bp_per_snp = 14322 (from your earlier run)."
  AVG_BP_PER_SNP=14322
fi

# compute recommended window (bp) to capture TARGET_SNPS
RECOMMENDED_WINDOW_BP=$(awk -v a="$AVG_BP_PER_SNP" -v t="$TARGET_SNPS" 'BEGIN{w=int(a*t+0.5); if(w<1) w=1; print w}')
# clamp
if [[ "$RECOMMENDED_WINDOW_BP" -lt "$MIN_WINDOW_BP" ]]; then RECOMMENDED_WINDOW_BP=$MIN_WINDOW_BP; fi
if [[ "$RECOMMENDED_WINDOW_BP" -gt "$MAX_WINDOW_BP" ]]; then RECOMMENDED_WINDOW_BP=$MAX_WINDOW_BP; fi

WINDOW_BP="${RECOMMENDED_WINDOW_BP}"
STEP_BP=$(( WINDOW_BP / 2 ))

echo "[INFO] Using window size (bp): $WINDOW_BP  (target ~${TARGET_SNPS} SNPs/window)"
echo "[INFO] Using step size (bp):   $STEP_BP"
echo

# sample lists expected
BKN_LIST="BKN_samples_onecol.keep"
HYD_LIST="HYD_samples_onecol.keep"
EG35_LIST="EG35_samples_onecol.keep"

declare -A POP_SAM
if [[ -f "$BKN_LIST" ]]; then POP_SAM["BKN"]="$BKN_LIST"; fi
if [[ -f "$HYD_LIST" ]]; then POP_SAM["HYD"]="$HYD_LIST"; fi
if [[ -f "$EG35_LIST" ]]; then POP_SAM["EG35"]="$EG35_LIST"; fi

# if not all present, try to construct from renamed_data_plink.fam or fallback to ALL
if [[ ${#POP_SAM[@]} -lt 3 ]]; then
  if [[ -f "renamed_data_plink.fam" ]]; then
    echo "[INFO] Not all population sample lists found; creating ALL from renamed_data_plink.fam"
    awk '{print $2}' renamed_data_plink.fam > "${OUTDIR}/ALL_samples_from_fam.keep"
    POP_SAM["ALL"]="${OUTDIR}/ALL_samples_from_fam.keep"
  else
    echo "[ERROR] Population sample lists missing and renamed_data_plink.fam not available." >&2
    exit 1
  fi
else
  # also create an ALL list
  cat "${BKN_LIST}" "${HYD_LIST}" "${EG35_LIST}" | sort | uniq > "${OUTDIR}/ALL_samples.keep"
  POP_SAM["ALL"]="${OUTDIR}/ALL_samples.keep"
fi

# Run vcftools for each population
echo "[INFO] Running vcftools for populations: ${!POP_SAM[@]}"
if ! command -v vcftools >/dev/null 2>&1 ; then
  echo "[ERROR] vcftools not found in PATH. Install vcftools and re-run." >&2
  exit 1
fi

for pop in "${!POP_SAM[@]}"; do
  sampfile="${POP_SAM[$pop]}"
  if [[ ! -s "$sampfile" ]]; then
    echo "[WARN] Sample list $sampfile for pop $pop is empty or missing. Skipping." >&2
    continue
  fi
  outbase="${OUTDIR}/${pop}"
  echo "[RUN] vcftools --gzvcf $VCF --keep $sampfile --window-pi $WINDOW_BP --window-pi-step $STEP_BP --out $outbase"
  vcftools --gzvcf "$VCF" --keep "$sampfile" --window-pi "$WINDOW_BP" --window-pi-step "$STEP_BP" --out "$outbase" >/dev/null
  echo "[RUN] vcftools --gzvcf $VCF --keep $sampfile --site-pi --out $outbase"
  vcftools --gzvcf "$VCF" --keep "$sampfile" --site-pi --out "$outbase" >/dev/null
  echo "[INFO] Done $pop: ${outbase}.windowed.pi  and  ${outbase}.sites.pi"
done

# Merge per-population windowed pi files
MERGED="${OUTDIR}/merged_window_pi.tsv"
echo -e "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI\tPOP" > "$MERGED"
for pop in $(printf "%s\n" "${!POP_SAM[@]}" | sort); do
  f="${OUTDIR}/${pop}.windowed.pi"
  if [[ -s "$f" ]]; then
    tail -n +2 "$f" | awk -v pop="$pop" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,pop}' >> "$MERGED"
  fi
done
echo "[RESULT] Merged windowed pi written to: $MERGED"

# Per-population summary CSV
SUMMARY_CSV="${OUTDIR}/pi_window_summary_by_pop.csv"
echo -e "POP,N_windows,mean_pi,median_pi,sd_pi" > "$SUMMARY_CSV"
for pop in $(printf "%s\n" "${!POP_SAM[@]}" | sort); do
  f="${OUTDIR}/${pop}.windowed.pi"
  if [[ -s "$f" ]]; then
    stats=$(tail -n +2 "$f" | awk '{if($5!="nan"){a[count++]=$5}} END{if(count==0){print "0,NA,NA,NA"} else {sum=0; for(i=0;i<count;i++) sum+=a[i]; mean=sum/count; asort(a); if(count%2==1) med=a[int((count+1)/2)-1]; else med=(a[int(count/2)-1]+a[int(count/2)])/2; sumsq=0; for(i=0;i<count;i++){sumsq+=(a[i]-mean)^2}; sd=(count>1)?sqrt(sumsq/(count-1)):0; printf("%d,%.8f,%.8f,%.8f", count, mean, med, sd)}}')
    echo -e "${pop},${stats}" >> "$SUMMARY_CSV"
  else
    echo -e "${pop},0,NA,NA,NA" >> "$SUMMARY_CSV"
  fi
done
echo "[RESULT] Summary CSV: $SUMMARY_CSV"

# call R plotting script if present
if [[ -f "plot_pi.R" ]]; then
  echo "[INFO] Running R plotting script (plot_pi_v2.R)..."
  Rscript plot_pi.R "$MERGED" "$OUTDIR" "$SUMMARY_CSV" || echo "[WARN] R plotting script failed. Check R packages."
else
  echo "[WARN] plot_pi_v2.R not found in working directory. Skip plotting."
fi

echo "[DONE] All pi computations complete. Outputs in: $OUTDIR"
