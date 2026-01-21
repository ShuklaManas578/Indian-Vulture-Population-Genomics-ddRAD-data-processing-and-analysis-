#!/usr/bin/env bash


# 18 Genetic Ancestry (ADMIXTURE) 
set -euo pipefail
IFS=$'\n\t'

# run_admixture.sh

WORKDIR="${WORKDIR:-$PWD}"

LD_OUT="${LD_OUT:-${WORKDIR}/path/to/ld_prune_results}"
OUTDIR="${OUTDIR:-${WORKDIR}/admix_runs}"

PRUNE_LABELS=( "relaxed" "stringent" )
ADMIX_CMD="${ADMIX_CMD:-admixture}"      # admixture binary 
KS=(1 2 3 4 5 6)                           # K values
NREPS=${NREPS:-5}
THREADS=${THREADS:-8}
SEED_BASE=${SEED_BASE:-13345}
TMPROOT="${TMPROOT:-/tmp/admixture_work}"   
COPY_FAM=true  

mkdir -p "$OUTDIR"
mkdir -p "$TMPROOT"
log(){ printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

if ! command -v "${ADMIX_CMD%% *}" >/dev/null 2>&1; then
  log "ERROR: admixture executable not found at '${ADMIX_CMD}'. Adjust ADMIX_CMD or install admixture."
  exit 1
fi

log "START ADMIXTURE with enforced integer CHR"
log "LD pruned dir: $LD_OUT"
log "Output dir: $OUTDIR"
log "Ks: ${KS[*]}  reps: $NREPS  threads: $THREADS"


MAP_AND_RUNS="${OUTDIR}/map_and_runs"
mkdir -p "$MAP_AND_RUNS"

# helper: produce an integer-mapped .bim from any .bim 

make_integer_bim() {
  local orig_bim="$1"; local out_int_bim="$2"; local out_map="$3"

  if [[ ! -f "$orig_bim" ]]; then
    echo "ERR: orig .bim not found: $orig_bim" >&2
    return 2
  fi

  # Build contig ordering (unique in order of first occurrence)
  awk '{ if(!seen[$1]++){ print $1 } }' "$orig_bim" > "${out_map}.contigs"

  # write mapping: original -> integer index 
  awk 'NR==FNR{cont[++n]=$1; next} END{ for(i=1;i<=n;i++) print cont[i]"\t"i }' "${out_map}.contigs" "${out_map}.contigs" > "$out_map"

  # Create integer-mapped bim (replace col1 using map)
  awk -v OFS="\t" 'NR==FNR{map[$1]=$2; next} { if($1 in map) $1=map[$1]; print }' "$out_map" "$orig_bim" > "$out_int_bim"


  if ! awk '{ if ($1 !~ /^[0-9]+$/) { exit 1 } } END { exit 0 }' "$out_int_bim"; then
    echo "ERR: integer-mapping failed: $out_int_bim still has non-integer CHR codes" >&2
    return 3
  fi

  return 0
}

# helper: create integer-copy of a PLINK bfile (bim->int.bim; copy bed/fam)
make_integer_bfile_copy() {
  local orig_pref="$1" out_pref="$2" out_map="$3"
  local orig_bim="${orig_pref}.bim"
  local orig_bed="${orig_pref}.bed"
  local orig_fam="${orig_pref}.fam"

  if [[ ! -f "$orig_bim" || ! -f "$orig_bed" || ! -f "$orig_fam" ]]; then
    echo "ERR: missing original bfile components for $orig_pref (checked ${orig_bim}, ${orig_bed}, ${orig_fam})" >&2
    return 2
  fi

  make_integer_bim "$orig_bim" "${out_pref}.bim" "$out_map" || return $?


  cp -f "$orig_bed" "${out_pref}.bed"
  cp -f "$orig_fam" "${out_pref}.fam"

  return 0
}

# Run admixture on a given integer bfile prefix
run_admixture_on_bfile() {
  local int_bprefix="$1"    
  local label="$2"          
  local mapfile="$3"        

  if [[ ! -f "${int_bprefix}.bed" ]]; then
    log "SKIP: integer bfile missing for ${label}: ${int_bprefix}.*"
    return
  fi

  local nvar nsmpl
  nvar=$(wc -l < "${int_bprefix}.bim" | tr -d ' ')
  nsmpl=$(wc -l < "${int_bprefix}.fam" | tr -d ' ')
  log "Preparing ADMIXTURE for ${label}: variants=${nvar}, samples=${nsmpl}"

  # run K/reps
  for K in "${KS[@]}"; do
    for ((rep=1; rep<=NREPS; rep++)); do
      local seed=$((SEED_BASE + K*100 + rep))
      local repdir="${OUTDIR}/${label}/K${K}/rep${rep}"
      mkdir -p "$repdir"
      local tmpd
      tmpd="$(mktemp -d "${TMPROOT}/${label}.K${K}.rep${rep}.XXXX")"

      # copy integer bfile to tmpd as work.bed/bim/fam
      cp -f "${int_bprefix}.bed" "${tmpd}/work.bed"
      cp -f "${int_bprefix}.bim" "${tmpd}/work.bim"
      cp -f "${int_bprefix}.fam" "${tmpd}/work.fam"

      # verify integer CHR right before admixture
      if ! awk 'NR>0{ if($1 !~ /^[0-9]+$/){ exit 1 } } END{ exit 0 }' "${tmpd}/work.bim"; then
        log "ERROR: tmp work.bim contains non-integer CHR codes - aborting this run (label=${label} K=${K} rep=${rep})"
        echo "tmp work.bim head:" > "${repdir}/work.bim.head.txt"
        head -n 20 "${tmpd}/work.bim" >> "${repdir}/work.bim.head.txt" || true
        rm -rf "$tmpd"
        continue
      fi

      log "  Running ADMIXTURE: ${label} K=${K} rep=${rep} seed=${seed}  (tmp: ${tmpd})"
      log "    tmp variants = $(wc -l < "${tmpd}/work.bim" | tr -d ' ')"

      pushd "$tmpd" >/dev/null
      set +e
      "${ADMIX_CMD}" --cv=10 --seed="$seed" -j"$THREADS" "work.bed" "$K" > admixture.stdout 2> admixture.stderr
      rc=$?
      set -e
      # move outputs regardless of rc
      mv -f work.* "${repdir}/" 2>/dev/null || true
      [[ -f admixture.stdout ]] && mv -f admixture.stdout "${repdir}/"
      [[ -f admixture.stderr ]] && mv -f admixture.stderr "${repdir}/"

      
      cp -f "$mapfile" "${repdir}/" 2>/dev/null || true
      if [[ "$COPY_FAM" = true ]]; then
        cp -f "${int_bprefix}.fam" "${repdir}/" 2>/dev/null || true
      fi

      # check created P/Q
      if [[ -f "${repdir}/work.P" && -f "${repdir}/work.Q" ]]; then
        log "    ADMIXTURE OK: produced work.P & work.Q -> ${repdir}"
      else
        log "    ADMIXTURE FAILED or produced no P/Q (rc=${rc}). See ${repdir}/admixture.stderr (first 40 lines):"
        if [[ -f "${repdir}/admixture.stderr" ]]; then
          sed -n '1,40p' "${repdir}/admixture.stderr" | sed 's/^/      /' || true
        else
          log "      (no admixture.stderr found)"
        fi
      fi

      popd >/dev/null
      rm -rf "$tmpd"
    done
  done
}

for lab in "${PRUNE_LABELS[@]}"; do
  # candidate prefixes to try 
  cand1="${LD_OUT}/pruned_${lab}"
  cand2="${LD_OUT}/pruned_${lab}.bfile"
  cand3="${LD_OUT}/filtered_final_snps_biallelic_pruned_${lab}.bfile"
  cand4="${LD_OUT}/filtered_final_snps_biallelic_pruned_${lab}.bfile.bfile"
  cand5="${LD_OUT}/pruned_${lab}.bfile.bfile"   
  orig_pref=""
  for c in "${cand1}" "${cand2}" "${cand3}" "${cand4}" "${cand5}"; do
    if [[ -f "${c}.bim" && -f "${c}.bed" && -f "${c}.fam" ]]; then
      orig_pref="$c"
      break
    fi
  done

  if [[ -z "$orig_pref" ]]; then
    log "SKIP: pruned bfile missing for label '${lab}'. Looked for candidates: ${cand1}, ${cand2}, ${cand3}"
    continue
  fi

  # prepare out prefixes
  int_pref="${MAP_AND_RUNS}/pruned_${lab}.bfile.int"
  mapfile="${MAP_AND_RUNS}/contig_map.pruned_${lab}.tsv"
  mkdir -p "$(dirname "$int_pref")"

  log "Processing label: pruned_${lab}"
  log "  Original bfile: ${orig_pref}.*"
  log "  Will write integer copy: ${int_pref}.* (map -> ${mapfile})"


  make_integer_bfile_copy "${orig_pref}" "${int_pref}" "${mapfile}" || {
    log "  ERROR: failed to produce integer-bfile for ${orig_pref}. Check ${mapfile}*.contigs"
    continue
  }

  # final sanity check
  if ! awk 'NR>0{ if($1 !~ /^[0-9]+$/){ exit 1 } } END{ exit 0 }' "${int_pref}.bim"; then
    log "  FATAL: integer mapping did not produce numeric CHR column in ${int_pref}.bim - aborting label ${lab}"
    sed -n '1,40p' "${int_pref}.bim" > "${OUTDIR}/debug_${lab}.bim.head.txt" || true
    continue
  fi

  run_admixture_on_bfile "${int_pref}" "pruned_${lab}" "${mapfile}"

  noev_pref="${LD_OUT}/pruned_${lab}.noEV_24"
  if [[ -f "${noev_pref}.bim" && -f "${noev_pref}.bed" && -f "${noev_pref}.fam" ]]; then
    int_pref2="${MAP_AND_RUNS}/pruned_${lab}.bfile.noEV_24.int"
    mapfile2="${MAP_AND_RUNS}/contig_map.pruned_${lab}.noEV_24.tsv"
    make_integer_bfile_copy "${noev_pref}" "${int_pref2}" "${mapfile2}" || {
      log "  WARNING: failed integer copy for ${noev_pref}; skipping"
    } || true
    run_admixture_on_bfile "${int_pref2}" "pruned_${lab}.bfile_noEV_24" "${mapfile2}"
  else
    log "  noEV_24 subset not present for ${lab} -> skip"
  fi
done

log "DONE. Results & integer copies under: $OUTDIR"
log "If any run failed, inspect the corresponding admixture.stderr under each rep dir."
