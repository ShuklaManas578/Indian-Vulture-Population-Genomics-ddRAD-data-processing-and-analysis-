#!/usr/bin/env bash


# 17 Pairwise Weir and Cockerham Fst
set -euo pipefail
IFS=$'\n\t'

# fst_from_pruned_bfiles.sh


PRUNED_DIR="${PRUNED_DIR:-Path/to/pruned/directory}"
POPMAP="${POPMAP:-path/to/popmap.txt}"
OUTDIR="${OUTDIR:-${PRUNED_DIR}/fst_from_pruned_sets}"
THREADS="${THREADS:-8}"
VERBOSE=1

BCFTOOLS="$(command -v bcftools || true)"
PLINK1="$(command -v plink || true)"
PLINK2="$(command -v plink2 || true)"
VCFTOOLS="$(command -v vcftools || true)"
TABIX="$(command -v tabix || true)"

log(){
  if [[ "$VERBOSE" -eq 1 ]]; then
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
  else
    printf '%s\n' "$*"
  fi
}

# basic checks
for cmd in "$BCFTOOLS" "$VCFTOOLS" "$TABIX"; do
  if [[ -z "$cmd" ]]; then
    echo "ERROR: required tool missing (bcftools, vcftools, tabix)"; exit 1
  fi
done

mkdir -p "$OUTDIR"/logs
log "PRUNED_DIR = $PRUNED_DIR"
log "POPMAP = $POPMAP"
log "OUTDIR = $OUTDIR"

# Read popmap and build poplists
if [[ ! -f "$POPMAP" ]]; then
  echo "ERROR: POPMAP not found at $POPMAP" >&2
  exit 1
fi

# Clean popmap (keep first two cols)
awk 'NF>=2{print $1"\t"$2}' "$POPMAP" > "${OUTDIR}/popmap.cleaned.txt"

# Unique populations (preserve order)
mapfile -t POPS < <(awk 'NF>=2{print $2}' "${OUTDIR}/popmap.cleaned.txt" | awk '!seen[$0]++')
if [[ ${#POPS[@]} -lt 2 ]]; then
  echo "ERROR: Need at least 2 populations in popmap" >&2
  exit 1
fi
log "Populations found: ${POPS[*]}"

mkdir -p "${OUTDIR}/poplists"
for p in "${POPS[@]}"; do
  awk -v P="$p" '$2==P{print $1}' "${OUTDIR}/popmap.cleaned.txt" > "${OUTDIR}/poplists/${p}.samples.txt"
  n=$(wc -l < "${OUTDIR}/poplists/${p}.samples.txt" | tr -d ' ')
  log "pop ${p}: ${n} samples"
done

# find pruned prefixes (detect pruned*.bim files)
mapfile -t BIM_FILES < <(find "$PRUNED_DIR" -maxdepth 1 -type f -name "pruned*.bim" -print | sort)
if [[ ${#BIM_FILES[@]} -eq 0 ]]; then
  echo "ERROR: No pruned .bim files found in $PRUNED_DIR (expected pruned_relaxed.bim / pruned_stringent.bim)" >&2
  exit 1
fi

# CSV summary header
SUMMARY_CSV="${OUTDIR}/fst_summary_by_pair.csv"
mkdir -p "$(dirname "$SUMMARY_CSV")"
echo "set,label,vcf_sites,total_samples,paired_popA,paired_popB,popA_samples,popB_samples,vcftools_sites,mean_fst,vcftools_status,vcftools_stderr_snippet" > "$SUMMARY_CSV"

# helper: compute mean fst from vcftools .weir.fst (last numeric column)
compute_mean_fst() {
  local fstfile="$1"
  if [[ ! -f "$fstfile" ]]; then echo "NA"; return; fi
  awk '{
    if($0 ~ /^#/) next
    f=$NF
    gsub("NaN","nan",f)
    if(f=="nan" || f=="NaN" || f=="") next
    if(f ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) { sum += f; n++ }
  } END{ if(n>0) printf("%.6g", sum/n); else print "NA" }' "$fstfile"
}

# process each pruned prefix
for bim in "${BIM_FILES[@]}"; do
  prefix="${bim%.bim}"           # e.g. .../pruned_relaxed
  label="$(basename "$prefix")"  # pruned_relaxed
  set_outdir="${OUTDIR}/${label}"
  mkdir -p "$set_outdir" "${set_outdir}/logs"
  log "------------------------------------------------------------"
  log "Processing pruned set label: $label"
  log "PLINK prefix: $prefix"

  existing_vcf=""
  if [[ -f "${prefix}.vcf.gz" ]]; then existing_vcf="${prefix}.vcf.gz"; fi
  if [[ -z "$existing_vcf" && -f "${prefix}.vcf" ]]; then existing_vcf="${prefix}.vcf"; fi
  out_vcf="${set_outdir}/${label}.pruned.vcf.gz"

  need_make_vcf=1
  if [[ -n "$existing_vcf" ]]; then
    log "Found existing VCF for $label: $existing_vcf -- copying to set output dir"
    cp -f "$existing_vcf" "$out_vcf"
    # create index if missing
    if [[ ! -f "${out_vcf}.tbi" ]]; then
      set +e
      "$TABIX" -f -p vcf "$out_vcf" >/dev/null 2>&1
      set -e
    fi
    need_make_vcf=0
  else
    if [[ -f "$out_vcf" && -s "$out_vcf" ]]; then
      sites=$(bcftools view -H "$out_vcf" 2>/dev/null | wc -l || echo 0)
      if [[ "$sites" -gt 0 ]]; then
        log "Existing VCF found in output: $out_vcf (sites: $sites). Will still ensure sample names / header"
        need_make_vcf=0
      fi
    fi
  fi

  # If necessary, try to create VCF from PLINK prefix
  if [[ $need_make_vcf -eq 1 ]]; then
    log "No existing pruned VCF. Attempting to create from PLINK prefix: $prefix -> $out_vcf"
    created=0

    # try plink1 first 
    if [[ -n "$PLINK1" && -x "$PLINK1" ]]; then
      log "Trying plink1 recode -> vcf bgz"
      set +e
      "${PLINK1}" --allow-extra-chr --bfile "$prefix" --out "${set_outdir}/${label}.plink_recoded" --recode vcf bgz > "${set_outdir}/logs/${label}.plink_recoded.stdout" 2> "${set_outdir}/logs/${label}.plink_recoded.stderr"
      rc=$?
      set -e
      if [[ $rc -eq 0 && -f "${set_outdir}/${label}.plink_recoded.vcf.gz" ]]; then
        mv -f "${set_outdir}/${label}.plink_recoded.vcf.gz" "$out_vcf"
        created=1
        log "Created pruned VCF using plink1: $out_vcf"
      else
        log "plink1 failed for $label; rc=$rc"
      fi
    fi

    # plink2 fallback
    if [[ $created -eq 0 && -n "$PLINK2" && -x "$PLINK2" ]]; then
      log "Trying plink2 recode -> vcf bgz"
      set +e
      "${PLINK2}" --allow-extra-chr --bfile "$prefix" --out "${set_outdir}/${label}.plink2_recoded" --recode vcf bgz > "${set_outdir}/logs/${label}.plink2_recoded.stdout" 2> "${set_outdir}/logs/${label}.plink2_recoded.stderr"
      rc=$?
      set -e
      if [[ $rc -eq 0 && -f "${set_outdir}/${label}.plink2_recoded.vcf.gz" ]]; then
        mv -f "${set_outdir}/${label}.plink2_recoded.vcf.gz" "$out_vcf"
        created=1
        log "Created pruned VCF using plink2: $out_vcf"
      else
        log "plink2 failed for $label; rc=$rc"
      fi
    fi

    if [[ $created -eq 0 ]]; then
      log "[ERROR] Failed to create VCF from PLINK prefix: $prefix (see logs ${set_outdir}/logs)"
      echo "[WARN] Failed to produce VCF from PLINK for $prefix. Skipping this set." > "${set_outdir}/logs/${label}.error"
      continue
    fi
  fi

  # check VCF exists
  if [[ ! -f "$out_vcf" ]]; then
    log "[WARN] Expected VCF missing: $out_vcf. Skipping."
    continue
  fi

  # ensure VCF is indexed (tabix)
  set +e
  "$TABIX" -f -p vcf "$out_vcf" >/dev/null 2>&1
  set -e

  # reheader samples to IID (use .fam)
  famfile="${prefix}.fam"
  if [[ -f "$famfile" ]]; then
    reheader_names="${set_outdir}/logs/${label}.new_iid_names.txt"
    awk '{print $2}' "$famfile" > "$reheader_names"
    vcf_count=$(bcftools query -l "$out_vcf" | wc -l | tr -d ' ')
    fam_count=$(wc -l < "$reheader_names" | tr -d ' ')
    if [[ "$vcf_count" -ne "$fam_count" ]]; then
      log "[WARN] sample count mismatch between VCF ($vcf_count) and .fam ($fam_count). Proceeding but check names."
    fi
    log "[INFO] Reheadering VCF sample names using IID from $famfile"
    set +e
    "$BCFTOOLS" reheader -s "$reheader_names" -o "${out_vcf}.reheadered.vcf.gz" "$out_vcf" 2> "${set_outdir}/logs/${label}.bcf_reheader.stderr"
    rc=$?
    set -e
    if [[ $rc -eq 0 && -f "${out_vcf}.reheadered.vcf.gz" ]]; then
      mv -f "${out_vcf}.reheadered.vcf.gz" "$out_vcf"
      "$TABIX" -f -p vcf "$out_vcf"
      log "[INFO] Reheadered VCF (sample names replaced with IID)."
    else
      log "[WARN] bcftools reheader failed; see ${set_outdir}/logs/${label}.bcf_reheader.stderr"
    fi
  else
    log "[WARN] No .fam file for $prefix; cannot reheader to IID."
  fi

  # if header fileformat=VCFv4.3, downgrade to 4.2 for vcftools 
  header_fmt="$($BCFTOOLS view -h "$out_vcf" | grep -i '^fileformat' | sed 's/.*=//I' || true)"
  if [[ -n "$header_fmt" && "$header_fmt" =~ ^VCFv4\.3$ ]]; then
    log "[INFO] VCF header says VCFv4.3; downgrading to VCFv4.2 so vcftools accepts it."
    hdrtmp="${set_outdir}/logs/${label}.header.tmp.txt"
    "$BCFTOOLS" view -h "$out_vcf" > "$hdrtmp"
    awk 'BEGIN{done=0} { if(!done && tolower($0) ~ /^fileformat=/) { sub(/=.*/,"=VCFv4.2",$0); done=1 } print }' "$hdrtmp" > "${hdrtmp}.mod"
    set +e
    "$BCFTOOLS" reheader -h "${hdrtmp}.mod" -o "${out_vcf}.fmtfix.vcf.gz" "$out_vcf" 2> "${set_outdir}/logs/${label}.bcf_reheader_fmt.stderr"
    rc=$?
    set -e
    if [[ $rc -eq 0 && -f "${out_vcf}.fmtfix.vcf.gz" ]]; then
      mv -f "${out_vcf}.fmtfix.vcf.gz" "$out_vcf"
      "$TABIX" -f -p vcf "$out_vcf"
      log "[INFO] Downgraded fileformat to VCFv4.2"
    else
      log "[WARN] Could not change header fileformat; vcftools may fail (see ${set_outdir}/logs/${label}.bcf_reheader_fmt.stderr)"
    fi
  fi

 
  if ! $BCFTOOLS view -h "$out_vcf" | grep -q '^##FORMAT=<ID=GT'; then
    log "[WARN] VCF appears to lack GT FORMAT field; vcftools needs genotypes. Will continue but vcftools may fail."
  fi

  # count sites and total samples
  sites_total=$(bcftools view -H "$out_vcf" | wc -l | tr -d ' ')
  total_samples=$(bcftools query -l "$out_vcf" | wc -l | tr -d ' ')
  log "[INFO] fileformat header: $($BCFTOOLS view -h "$out_vcf" | grep -i '^fileformat' || true)"
  log "[INFO] sites in VCF to be used for vcftools: $sites_total (file: $out_vcf)"
  log "[INFO] total samples in VCF: $total_samples"

  # prepare sample intersection per-pop and run pairwise vcftools
  for ((i=0;i<${#POPS[@]};i++)); do
    for ((j=i+1;j<${#POPS[@]};j++)); do
      A=${POPS[i]}; B=${POPS[j]}
      popA="${OUTDIR}/poplists/${A}.samples.txt"
      popB="${OUTDIR}/poplists/${B}.samples.txt"
      outpref="${set_outdir}/fst_${A}_${B}"
      log "[INFO] vcftools: ${A} vs ${B} -> ${outpref}.weir.fst"

      # check overlap between vcf samples and poplist
      mapfile -t vcf_samples < <($BCFTOOLS query -l "$out_vcf")
      mapfile -t popA_samples < <(awk 'NF{print $1}' "$popA")
      mapfile -t popB_samples < <(awk 'NF{print $1}' "$popB")
      # compute intersection counts
      tmp_vcf_sorted=$(mktemp)
      printf "%s\n" "${vcf_samples[@]}" | sort > "${tmp_vcf_sorted}"
      tmp_A_sorted=$(mktemp); printf "%s\n" "${popA_samples[@]}" | sort > "${tmp_A_sorted}"
      tmp_B_sorted=$(mktemp); printf "%s\n" "${popB_samples[@]}" | sort > "${tmp_B_sorted}"
      nA=$(comm -12 "${tmp_vcf_sorted}" "${tmp_A_sorted}" | wc -l | tr -d ' ')
      nB=$(comm -12 "${tmp_vcf_sorted}" "${tmp_B_sorted}" | wc -l | tr -d ' ')
      rm -f "${tmp_vcf_sorted}" "${tmp_A_sorted}" "${tmp_B_sorted}"
      log "[INFO] pop ${A}: $nA samples available in VCF; pop ${B}: $nB samples available in VCF"

      # run vcftools
      set +e
      "$VCFTOOLS" --gzvcf "$out_vcf" --weir-fst-pop "$popA" --weir-fst-pop "$popB" --out "$outpref" > "${outpref}.vcftools.stdout" 2> "${outpref}.vcftools.stderr"
      rc=$?
      set -e

      vcftools_status="OK"
      stderr_snip=""
      if [[ $rc -ne 0 ]]; then
        vcftools_status="FAIL"
        stderr_snip="$(sed -n '1,6p' "${outpref}.vcftools.stderr" | tr '\n' ' ' | sed 's/"/\\"/g')"
        log "[WARN] vcftools failed for ${A} vs ${B}; snippet: ${stderr_snip}"
      fi

      fstfile="${outpref}.weir.fst"
      if [[ -f "$fstfile" ]]; then
        fst_sites=$(wc -l < "$fstfile" | tr -d ' ')
        meanfst="$(compute_mean_fst "$fstfile")"
      else
        fst_sites=0
        meanfst="NA"
      fi

      # append to summary CSV
      echo "${label},${label},${sites_total},${total_samples},${A},${B},${nA},${nB},${fst_sites},${meanfst},\"${vcftools_status}\",\"${stderr_snip}\"" >> "$SUMMARY_CSV"
    done
  done

  log "[INFO] Done with set $label. Outputs in: $set_outdir"

done

log "ALL DONE. Summary written to: $SUMMARY_CSV"
