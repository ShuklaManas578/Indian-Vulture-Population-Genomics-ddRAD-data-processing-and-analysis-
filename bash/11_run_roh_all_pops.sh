#!/usr/bin/env bash


#7 Runs of Homozygosity Quantification
# run_roh_all_pops.sh
# Run PLINK ROH (--homozyg) for populations and summarize results.

set -euo pipefail
IFS=$'\n\t'

WORKDIR=${1:-$(pwd)}
cd "$WORKDIR"

# PLINK homozygosity parameters
HOM_PARAMS=(
  "--allow-extra-chr"
  "--homozyg"
  "--homozyg-density" "50"
  "--homozyg-gap" "1000"
  "--homozyg-het" "1"
  "--homozyg-kb" "1000"
  "--homozyg-snp" "50"
  "--homozyg-window-het" "1"
  "--homozyg-window-snp" "50"
  "--homozyg-window-missing" "5"
  "--homozyg-window-threshold" "0.05"
)

# populations to process
POPS=("BKN" "HYD" "BKN_35")

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="roh_results_${TIMESTAMP}"
mkdir -p "$OUTDIR"

echo "Working dir: $WORKDIR"
echo "Output dir:  $OUTDIR"
echo "PLINK homozyg params: ${HOM_PARAMS[*]}"
echo

# find bfile prefix helper
find_bfile_prefix() {
  local pop="$1"
  local candidates=(
    "renamed_data_plink_${pop}"
    "renamed_data_plink_${pop//_/}"
    "renamed_data_plink_${pop^^}"
  )
  for cand in "${candidates[@]}"; do
    if [[ -f "${cand}.bed" && -f "${cand}.bim" && -f "${cand}.fam" ]]; then
      echo "$cand"
      return 0
    fi
  done
  return 1
}

COMBINED_SUMMARY="$OUTDIR/combined_roh_summary.tsv"
echo -e "POP\tFID\tIID\tn_roh_segments\ttotal_roh_kb\tmean_roh_kb_per_segment" > "$COMBINED_SUMMARY"

for pop in "${POPS[@]}"; do
  echo "=== POP: $pop ==="
  bprefix=$(find_bfile_prefix "$pop" || true)
  if [[ -z "$bprefix" ]]; then
    echo "WARNING: could not find PLINK bfile for pop '$pop. Skipping."
    echo "POP $pop: MISSING_BFILE" >> "$OUTDIR/roh_run_status.txt"
    continue
  fi
  echo "Using bfile prefix: $bprefix"

  out_pref="${OUTDIR}/${pop}_roh"
  log_file="${out_pref}.plink.log"

  echo "Running PLINK for $pop ..."
  plink --bfile "$bprefix" "${HOM_PARAMS[@]}" --out "$out_pref" > "$log_file" 2>&1 || {
    echo "ERROR: PLINK failed for $pop. See $log_file"
    echo "POP $pop: PLINK_FAILED" >> "$OUTDIR/roh_run_status.txt"
    continue
  }
  echo "PLINK finished for $pop. Log: $log_file"

  hom_segments="${out_pref}.hom"
  hom_indiv="${out_pref}.hom.indiv"

  if [[ -f "$hom_indiv" ]]; then
    echo "Found per-individual summary: $hom_indiv"
    # detect NSEG-like column and KB-like column
    nseg_i=$(awk 'NR==1{for(i=1;i<=NF;i++){h=tolower($i); if(h ~ /nseg|n_segments|n.segments|nsegments|nseg/i){print i; exit}}}' "$hom_indiv" || true)
    kb_i=$(awk 'NR==1{for(i=1;i<=NF;i++){h=tolower($i); if(h ~ /kb|length|kb_len|length_kb|sum_kb|total_kb/i){print i; exit}}}' "$hom_indiv" || true)

    if [[ -n "$nseg_i" && -n "$kb_i" && "$nseg_i" -gt 0 && "$kb_i" -gt 0 ]]; then
      # append to combined summary: read each data line and print POP, FID, IID, nseg, kb, mean
      awk -v POP="$pop" -v NIDX="$nseg_i" -v KIDX="$kb_i" 'NR>1{
         fid=$1; iid=$2;
         nseg = (NIDX<=NF ? $(NIDX) : 0);
         kb   = (KIDX<=NF ? $(KIDX) : 0);
         if(kb=="" || kb=="NA") kb=0;
         if(nseg=="" || nseg=="NA") nseg=0;
         mean_k = (nseg>0 ? kb / nseg : 0);
         # ensure numeric formatting
         printf("%s\t%s\t%s\t%d\t%.2f\t%.4f\n", POP, fid, iid, nseg+0, kb+0, mean_k+0);
      }' "$hom_indiv" >> "$COMBINED_SUMMARY"

      echo "Appended per-individual RoH summary from $hom_indiv"
    else
      echo "Could not auto-detect NSEG/KB columns in $hom_indiv; skipping that parsing."
    fi
  fi

  # 2) If .hom (segments) exists, compute per-sample totals from segments file
  if [[ -f "$hom_segments" ]]; then
    echo "Found segment file: $hom_segments"
    kb_col=$(awk 'NR==1{for(i=1;i<=NF;i++){h=tolower($i); if(h ~ /kb|length|kb_len|length_kb|seglen|segment_kb/) {print i; exit}}; print 0}' "$hom_segments")
    if [[ "$kb_col" -gt 0 ]]; then
      echo "Detected KB column at index: $kb_col"
      # compute per-sample summaries and append
      awk -v POP="$pop" -v KBCOL="$kb_col" 'NR>1{
         fid=$1; iid=$2;
         kb=$(KBCOL<=NF ? $(KBCOL) : 0);
         if(kb=="" || kb=="NA") kb=0;
         key = fid "\t" iid;
         segs[key] += 1;
         tot[key] += kb;
      }
      END{
         OFS="\t";
         for(k in segs){
           n = segs[k];
           total = tot[k];
           mean = (n>0 ? total / n : 0);
           split(k, parts, "\t");
           printf("%s\t%s\t%s\t%d\t%.2f\t%.4f\n", POP, parts[1], parts[2], n, total, mean);
         }
      }' "$hom_segments" >> "$COMBINED_SUMMARY"
      echo "Appended per-sample RoH summary derived from segments for $pop."
    else
      echo "Could not detect KB column in $hom_segments; skipping per-segment summarization."
    fi
  else
    echo "No .hom segments file found for $pop (expected ${out_pref}.hom)."
  fi

  echo "POP $pop: DONE" >> "$OUTDIR/roh_run_status.txt"
  echo
done

# create a small report
report="$OUTDIR/roh_summary_report.txt"
{
  echo "RoH run summary: $(date)"
  echo "Working dir: $WORKDIR"
  echo "Output dir: $OUTDIR"
  echo
  echo "Per-pop RoH summary (first 200 lines):"
  if [[ -s "$COMBINED_SUMMARY" ]]; then
    head -n 200 "$COMBINED_SUMMARY" | column -t -s $'\t'
  else
    echo "No combined summary produced."
  fi
  echo
  echo "PLINK individual run statuses (roh_run_status.txt):"
  if [[ -f "$OUTDIR/roh_run_status.txt" ]]; then
    cat "$OUTDIR/roh_run_status.txt"
  else
    echo "No status file found."
  fi
} > "$report"

echo "Done. Results and logs are in: $OUTDIR"
echo "Combined per-sample RoH TSV: $COMBINED_SUMMARY"
echo "Human-readable report: $report"
