#!/usr/bin/env bash


# 12 Tajima's D (non-overlapping)
# compute_tajimas_vcftools.sh

set -euo pipefail
IFS=$'\n\t'

VCF="Path to vcf file"
POPMAP="Path to popmap.txt file"
OUTDIR="./tajima_vcftools_results"
WINDOW_BP="${WINDOW_BP:-1432000}"   # recommended from SNP density (~1.432 Mb for ~100 SNPs)
VERBOSE=1
VCFTOOLS="$(command -v vcftools || true)"
BCFTOOLS="$(command -v bcftools || true)"
TABIX="$(command -v tabix || true)"

mkdir -p "$OUTDIR"/logs
log(){ if [[ "$VERBOSE" -eq 1 ]]; then printf '[%s] %s\n' "$(date '+%F %T')" "$*"; else printf '%s\n' "$*"; fi }

if [[ ! -x "$VCFTOOLS" ]]; then log "ERROR: vcftools not found in PATH"; exit 1; fi
if [[ ! -x "$BCFTOOLS" ]]; then log "ERROR: bcftools not found in PATH"; exit 1; fi

log "VCF: $VCF"
log "POPMAP: $POPMAP"
log "OUTDIR: $OUTDIR"
log "Window (bp): $WINDOW_BP"

# 1) prepare poplists from popmap (sample TAB pop)
if [[ ! -f "$POPMAP" ]]; then log "ERROR: popmap not found: $POPMAP"; exit 1; fi
awk 'NF>=2{print $1"\t"$2}' "$POPMAP" > "${OUTDIR}/popmap.cleaned.txt"
mapfile -t POPS < <(awk 'NF>=2{print $2}' "${OUTDIR}/popmap.cleaned.txt" | awk '!seen[$0]++')
for p in "${POPS[@]}"; do
  awk -v P="$p" '$2==P{print $1}' "${OUTDIR}/popmap.cleaned.txt" > "${OUTDIR}/${p}.samples.txt"
  log "Prepared poplist for $p : $(wc -l < "${OUTDIR}/${p}.samples.txt" | tr -d ' ') samples"
done

# 2) run vcftools TajimaD for full set (non-overlapping windows)
log "Running vcftools TajimaD for full (non-overlapping) windows (window = ${WINDOW_BP})..."
set +e
"$VCFTOOLS" --gzvcf "$VCF" --TajimaD "$WINDOW_BP" --out "${OUTDIR}/tajima_all" > "${OUTDIR}/logs/tajima_all.stdout" 2> "${OUTDIR}/logs/tajima_all.stderr"
rc=$?
set -e
if [[ $rc -ne 0 ]]; then log "WARNING: vcftools TajimaD for full dataset failed (check ${OUTDIR}/logs/tajima_all.stderr)"; fi

# 3) run per-pop TajimaD (non-overlapping)
for p in "${POPS[@]}"; do
  POPLIST="${OUTDIR}/${p}.samples.txt"
  outpref="${OUTDIR}/tajima_${p}"
  log "Running TajimaD for pop $p (window=${WINDOW_BP})"
  set +e
  "$VCFTOOLS" --gzvcf "$VCF" --keep "$POPLIST" --TajimaD "$WINDOW_BP" --out "$outpref" > "${OUTDIR}/logs/tajima_${p}.stdout" 2> "${OUTDIR}/logs/tajima_${p}.stderr"
  rc=$?
  set -e
  if [[ $rc -ne 0 ]]; then log "WARNING: vcftools failed for pop $p (see ${OUTDIR}/logs/tajima_${p}.stderr)"; fi
done

# 4) convert vcftools output (.Tajima.D) into tidy CSVs
to_csv(){
  infile="$1"; setname="$2"; outfile="$3"; window_size="$4"
  if [[ ! -f "$infile" ]]; then
    # write empty file w/header if missing
    echo "chrom,start,end,n_variants,tajimaD,set" > "$outfile"
    return
  fi
  # Fixed awk command - properly handle all fields and compute end position
  awk -v S="$setname" -v win="$window_size" '
    BEGIN {
      OFS=","
      print "chrom","start","end","n_variants","tajimaD","set"
    }
    NR>1 {
      chrom = $1
      start = $2
      # Calculate end position: start + window_size - 1
      end = start + win - 1
      n_variants = $3
      tajimaD = $4

      # Handle missing/invalid tajimaD values
      if (tajimaD == "nan" || tajimaD == "NaN" || tajimaD == "-nan") {
        tajimaD = "NA"
      }

      print chrom, start, end, n_variants, tajimaD, S
    }' "$infile" > "$outfile"
}

to_csv "${OUTDIR}/tajima_all.Tajima.D" "ALL" "${OUTDIR}/tajima_all.csv" "$WINDOW_BP"
for p in "${POPS[@]}"; do
  to_csv "${OUTDIR}/tajima_${p}.Tajima.D" "${p}" "${OUTDIR}/tajima_${p}.csv" "$WINDOW_BP"
done

# 5) combine all into master CSV
echo "chrom,start,end,n_variants,tajimaD,set" > "${OUTDIR}/tajima_allsets_combined.csv"
for f in "${OUTDIR}"/tajima_*.csv; do
  # skip the master
  if [[ "$f" == "${OUTDIR}/tajima_allsets_combined.csv" ]]; then continue; fi
  if [[ -f "$f" ]]; then
    tail -n +2 "$f" >> "${OUTDIR}/tajima_allsets_combined.csv"
  fi
done

log "Done. Outputs (per-pop and combined) in: $OUTDIR"
log "Non-overlapping (vcftools) TajimaD files: ${OUTDIR}/*.Tajima.D  CSVs: ${OUTDIR}/*.csv"
