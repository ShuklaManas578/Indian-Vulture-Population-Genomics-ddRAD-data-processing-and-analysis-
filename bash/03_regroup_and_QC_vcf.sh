#!/usr/bin/env bash


#5 Per-pop PLINK conversion and heterozygosity check
# regroup_and_QC_vcf.sh

set -euo pipefail
IFS=$'\n\t'

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 <input.vcf.gz> <output_dir> [reference.fa]"
  exit 1
fi

IN_VCF="$1"
OUT_DIR="$2"
REF_FA=${3-}

mkdir -p "$OUT_DIR"

# resolve paths
if command -v readlink &>/dev/null; then
  OUT_DIR=$(readlink -m "$OUT_DIR")
else
  OUT_DIR="$(cd "$OUT_DIR" && pwd)"
fi
IN_VCF=$(realpath "$IN_VCF")
if [[ -n "$REF_FA" ]]; then REF_FA=$(realpath "$REF_FA"); fi

TMPDIR=$(mktemp -d -p "$OUT_DIR" tmp.rename.XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

echo "Input VCF: $IN_VCF"
echo "Output dir: $OUT_DIR"
if [[ -n "$REF_FA" ]]; then echo "Reference fasta: $REF_FA"; fi
echo

for prog in bcftools bgzip tabix plink awk sed grep column realpath; do
  if ! command -v "$prog" &>/dev/null; then
    echo "Error: required program '$prog' not found in PATH." >&2
    exit 2
  fi
done

if [[ -n "$REF_FA" ]]; then
  if [[ ! -f "$REF_FA" ]]; then
    echo "Error: reference fasta '$REF_FA' not found." >&2
    exit 3
  fi
  if [[ ! -f "${REF_FA}.fai" ]]; then
    if command -v samtools &>/dev/null; then
      echo "Indexing reference with samtools faidx"
      samtools faidx "$REF_FA"
    else
      echo "Warning: ${REF_FA}.fai missing and samtools not available" >&2
    fi
  fi
fi

# 1) original sample list
ORIG_SAMPLE_LIST="$TMPDIR/orig_samples.txt"
bcftools query -l "$IN_VCF" > "$ORIG_SAMPLE_LIST"
if [[ ! -s "$ORIG_SAMPLE_LIST" ]]; then
  echo "No samples found in VCF or failed to read samples." >&2
  exit 4
fi

# 2) build mapping
MAPPING_TSV="$OUT_DIR/sample_mapping_old_to_new.tsv"
NEW_NAMES="$TMPDIR/new_names.txt"
: > "$MAPPING_TSV"
: > "$NEW_NAMES"

while read -r s; do
  new="$s"
  if [[ "$s" == "EG_35" ]]; then
    new="EG_35"
  elif [[ "$s" == WR* ]]; then
    new="${s/#WR/HYD}"
  elif [[ "$s" == EG* ]]; then
    new="${s/#EG/BKN}"
  elif [[ "$s" == EV* ]]; then
    new="${s/#EV/BKN}"
  else
    new="$s"
  fi
  printf "%s\t%s\n" "$s" "$new" >> "$MAPPING_TSV"
  printf "%s\n" "$new" >> "$NEW_NAMES"
done < "$ORIG_SAMPLE_LIST"

echo "Sample mapping written to: $MAPPING_TSV"
echo "---- sample mapping preview ----"
column -t "$MAPPING_TSV" | sed -n '1,200p'
echo

# 3) reheader -> force ASCII VCF
TMP_VCF="$OUT_DIR/$(basename "${IN_VCF%.*}").reheadered_uncompressed.vcf"
echo "Creating uncompressed reheadered VCF: $TMP_VCF"
if ! bcftools reheader -s "$NEW_NAMES" "$IN_VCF" | bcftools view -Ov -o "$TMP_VCF" - ; then
  echo "Error: bcftools reheader|view failed. Inspect $TMP_VCF (if present) and $MAPPING_TSV." >&2
  exit 5
fi

if ! head -n 10 "$TMP_VCF" | grep -q '^##fileformat=VCF'; then
  echo "Error: The reheadered file $TMP_VCF does not appear to be a valid VCF." >&2
  head -n 40 "$TMP_VCF" || true
  exit 6
fi

# 4) bgzip + tabix
OUT_VCF="$OUT_DIR/$(basename "${TMP_VCF}").gz"
echo "Compressing to bgzip: $OUT_VCF"
bgzip -c "$TMP_VCF" > "$OUT_VCF"
echo "Indexing with tabix..."
if ! tabix -p vcf "$OUT_VCF"; then
  echo "tabix failed on $OUT_VCF. Print first 80 lines of uncompressed file for debug:" >&2
  head -n 80 "$TMP_VCF" >&2 || true
  exit 7
fi

# 5) verify sample list
NEW_SAMPLE_LIST="$TMPDIR/new_samples_extracted.txt"
bcftools query -l "$OUT_VCF" > "$NEW_SAMPLE_LIST"
if diff -u "$NEW_NAMES" "$NEW_SAMPLE_LIST" > "$TMPDIR/name_diff.txt"; then
  echo "Renaming verification: OK"
else
  echo "Warning: mismatch between expected names and actual names in renamed VCF. See $TMPDIR/name_diff.txt" >&2
  sed -n '1,200p' "$TMPDIR/name_diff.txt"
fi

# 6) bcftools stats
BCF_STATS="$OUT_DIR/$(basename "${OUT_VCF}").bcftools.stats.txt"
BCF_STATS_SUMMARY="$OUT_DIR/$(basename "${OUT_VCF}").bcftools.stats.summary.txt"
echo "Running bcftools stats..."
if [[ -n "$REF_FA" ]]; then
  bcftools stats -F "$REF_FA" -s - "$OUT_VCF" > "$BCF_STATS"
else
  bcftools stats -s - "$OUT_VCF" > "$BCF_STATS"
fi
grep -E '^TST|^PSC' -n "$BCF_STATS" | sed -n '1,200p' > "$BCF_STATS_SUMMARY" || true
echo "bcftools stats saved to: $BCF_STATS"
echo

# 7) Convert to PLINK using the uncompressed TMP_VCF
PLINK_PREFIX="$OUT_DIR/renamed_data_plink"
echo "Converting uncompressed reheadered VCF to PLINK: prefix=$PLINK_PREFIX"
plink --vcf "$TMP_VCF" --double-id --allow-extra-chr --make-bed --out "$PLINK_PREFIX"

echo "Computing per-sample heterozygosity/inbreeding and missingness for ALL samples..."
plink --bfile "$PLINK_PREFIX" --allow-extra-chr --het --out "${PLINK_PREFIX}_het"
plink --bfile "$PLINK_PREFIX" --allow-extra-chr --missing --out "${PLINK_PREFIX}_missing"

awk '{print $2}' "$MAPPING_TSV" > "$TMPDIR/all_new_names_onecol.txt"
awk '{if($2 ~ /^BKN/) print $2}' "$MAPPING_TSV" > "$OUT_DIR/BKN_samples_onecol.txt"
awk '{if($2 ~ /^HYD/) print $2}' "$MAPPING_TSV" > "$OUT_DIR/HYD_samples_onecol.txt"
awk '{if($2 ~ /^EG_35$/) print $2}' "$MAPPING_TSV" > "$OUT_DIR/EG35_samples_onecol.txt"

# create PLINK keep files with FID IID (double-id -> both identical)
for f in "$OUT_DIR"/{BKN,HYD,EG35}_samples_onecol.txt; do
  if [[ -s "$f" ]]; then
    awk '{print $1, $1}' "$f" > "${f%.txt}.keep"
  fi
done

# 8) per-pop PLINKs + het/missing (all with allow-extra-chr)
for pop in BKN HYD EG35; do
  POP_KEEP="$OUT_DIR/${pop}_samples_onecol.keep"
  if [[ -s "$POP_KEEP" ]]; then
    POP_PREFIX="$OUT_DIR/renamed_data_plink_${pop}"
    echo "Creating PLINK dataset for $pop (keep=$POP_KEEP) -> $POP_PREFIX"
    plink --bfile "$PLINK_PREFIX" --allow-extra-chr --keep "$POP_KEEP" --make-bed --out "$POP_PREFIX" || true

    if [[ -f "${POP_PREFIX}.bed" ]]; then
      # compute het on the per-pop PLINK (from the per-pop bed)
      plink --bfile "$POP_PREFIX" --allow-extra-chr --het --out "${POP_PREFIX}_het" || true
      plink --bfile "$POP_PREFIX" --allow-extra-chr --missing --out "${POP_PREFIX}_missing" || true

      # ALSO compute het for this pop FROM THE ALL PLINK using --keep (ensures same variant set as ALL)
      plink --bfile "$PLINK_PREFIX" --allow-extra-chr --keep "$POP_KEEP" --het --out "${PLINK_PREFIX}_${pop}_het_fromALL" || true
    fi
  else
    echo "No samples for $pop (skip)."
  fi
done

# 9) summarize het ALL robustly
HET_TSV="$OUT_DIR/heterozygosity_inbreeding_summary_ALL.tsv"
if [[ -f "${PLINK_PREFIX}_het.het" ]]; then
  # robust awk: accept headers "N(NM)" or "N_SITES" or others; print FID IID N F
  awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1{
    nidx = h["N(NM)"]; if(!nidx) nidx = h["N_SITES"]; if(!nidx) nidx = h["N_SITES(NM)"];
    # fallback: if we do not find the header, assume column 5 is N and 6 is F (common PLINK order)
    if(!nidx) nidx = 5;
    fidx = h["F"]; if(!fidx) fidx = 6;
    print $h["FID"], $h["IID"], ($(nidx) ? $(nidx) : "NA"), ($(fidx) ? $(fidx) : "NA")
  }' "${PLINK_PREFIX}_het.het" > "$HET_TSV"
  echo "Wrote het summary to: $HET_TSV"
else
  echo "Warning: expected het file not found: ${PLINK_PREFIX}_het.het" >&2
fi

# 10) create comparable per-pop het summaries (from ALL and from per-pop PLINKs) and compare
for pop in BKN HYD EG35; do
  POP_KEEP="$OUT_DIR/${pop}_samples_onecol.keep"
  POP_PREFIX="$OUT_DIR/renamed_data_plink_${pop}"
  ALL_HET_FROMALL="${PLINK_PREFIX}_${pop}_het_fromALL.het"
  PERPOP_HET="${POP_PREFIX}_het.het"

  OUT_SUM_ALL="$OUT_DIR/${pop}_het_fromALL_summary.tsv"
  OUT_SUM_PERPOP="$OUT_DIR/${pop}_het_perpop_summary.tsv"
  COMPARE_OUT="$OUT_DIR/${pop}_het_compare.diff.tsv"

  if [[ -f "$ALL_HET_FROMALL" ]]; then
    awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1{
      nidx = h["N(NM)"]; if(!nidx) nidx=h["N_SITES"]; if(!nidx) nidx=5;
      fidx = h["F"]; if(!fidx) fidx=6;
      print $h["FID"],$h["IID"],$(nidx),$(fidx)
    }' "$ALL_HET_FROMALL" > "$OUT_SUM_ALL"
  fi

  if [[ -f "$PERPOP_HET" ]]; then
    awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++){h[$i]=i}} NR>1{
      nidx = h["N(NM)"]; if(!nidx) nidx=h["N_SITES"]; if(!nidx) nidx=5;
      fidx = h["F"]; if(!fidx) fidx=6;
      print $h["FID"],$h["IID"],$(nidx),$(fidx)
    }' "$PERPOP_HET" > "$OUT_SUM_PERPOP"
  fi

  if [[ -s "$OUT_SUM_ALL" && -s "$OUT_SUM_PERPOP" ]]; then
    # sort and join (both files have "FID IID N F")
    sort -k1,1 -k2,2 "$OUT_SUM_ALL" > "$OUT_SUM_ALL.sorted"
    sort -k1,1 -k2,2 "$OUT_SUM_PERPOP" > "$OUT_SUM_PERPOP.sorted"
    join -1 1 -2 1 -o '1.1 1.2 1.3 1.4 2.3 2.4' -t $'\t' <(awk '{print $1,$2,$3,$4}' "$OUT_SUM_ALL.sorted") <(awk '{print $1,$2,$3,$4}' "$OUT_SUM_PERPOP.sorted") \
      > "$COMPARE_OUT" || true
    if [[ -s "$COMPARE_OUT" ]]; then
      echo -e "Differences for $pop written to: $COMPARE_OUT\nColumns: FID IID N_fromALL F_fromALL N_fromPERPOP F_fromPERPOP"
      # print any sample where N differs
      awk -F'\t' 'NR>0 && ($3 != $5 || $4 != $6){print $0}' "$COMPARE_OUT" | sed -n '1,200p' || true
    else
      echo "No differences between per-pop (derived) and per-pop (from ALL) for $pop (or no data to compare)."
    fi
  fi
done

# README
cat > "$OUT_DIR/README_rename_qc.txt" <<EOF
Files produced:
- ${OUT_VCF}: renamed, bgzipped VCF (and .tbi)
- ${BCF_STATS}: bcftools stats
- ${BCF_STATS_SUMMARY}: quick bcftools TST/PSC summary
- ${PLINK_PREFIX}.bed/.bim/.fam : PLINK binary files (all)
- ${PLINK_PREFIX}_het.het : per-sample het/inbreeding (ALL)
- ${PLINK_PREFIX}_missing.imiss : per-sample missingness (ALL)
- per-pop PLINKs: renamed_data_plink_BKN.*, renamed_data_plink_HYD.*, renamed_data_plink_EG35.*
- per-pop het computed FROM per-pop PLINKs: renamed_data_plink_<POP>_het.het
- per-pop het computed FROM ALL using --keep: renamed_data_plink_<POP>_het_fromALL.het
- per-pop compare diffs: <POP>_het_compare.diff.tsv (columns: FID IID N_fromALL F_fromALL N_fromPERPOP F_fromPERPOP)
- sample_mapping_old_to_new.tsv
- heterozygosity_inbreeding_summary_ALL.tsv
EOF

echo
echo "Done. All outputs in: $OUT_DIR"
