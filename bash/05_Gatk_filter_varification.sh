#!/usr/bin/env bash


#4 Verifying filters for VCF and genotype summary generation

set -euo pipefail

MIN_DP=10
MIN_GQ=20
QUAL_MIN=30
QD_MIN=2.0
MQ_MIN=40.0
MAF_THRESHOLD=0.03
MIN_CALL_RATE=0.70
MAX_DP_CAP=2000
# -----------------------------------

WORKDIR="$(pwd)"
JOINT_OUT="${WORKDIR}/joint_out"
RAW_VCF="${JOINT_OUT}/raw_joint.vcf.gz"
FILTERED_VCF="${JOINT_OUT}/filtered_final.vcf.gz"
REPORT_DIR="${JOINT_OUT}/report_filters_check"
mkdir -p "$REPORT_DIR"

BCFTOOLS="${BCFTOOLS:-bcftools}"
TABIX="${TABIX:-tabix}"
SAMTOOLS="${SAMTOOLS:-samtools}"
PYTHON="${PYTHON:-python3}"

for t in "$BCFTOOLS" "$TABIX" "$SAMTOOLS" "$PYTHON"; do
  if ! command -v "$t" >/dev/null 2>&1; then
    echo "ERROR: required tool not found: $t" >&2
    exit 1
  fi
done

# files
if [[ ! -f "$RAW_VCF" ]]; then echo "ERROR: raw vcf not found: $RAW_VCF" >&2; exit 1; fi
if [[ ! -f "$FILTERED_VCF" ]]; then echo "ERROR: filtered vcf not found: $FILTERED_VCF" >&2; exit 1; fi

# ensure indexed
for f in "$RAW_VCF" "$FILTERED_VCF"; do
  if [[ ! -f "${f}.tbi" && ! -f "${f}.csi" ]]; then
    echo "Indexing $f ..."
    $TABIX -p vcf "$f"
  fi
done

echo "Writing reports to: $REPORT_DIR"
echo

# Basic counts
RAW_TOTAL=$($BCFTOOLS view -H "$RAW_VCF" | wc -l)
FILTERED_TOTAL=$($BCFTOOLS view -H "$FILTERED_VCF" | wc -l)
RAW_SAMPLES=$($BCFTOOLS query -l "$RAW_VCF" | wc -l)

cat > "${REPORT_DIR}/summary_basic.txt" <<EOF
Raw joint VCF: $RAW_VCF
Filtered joint VCF: $FILTERED_VCF

Raw variant count (lines): $RAW_TOTAL
Filtered variant count (lines): $FILTERED_TOTAL
Variants removed by filtering: $(( RAW_TOTAL - FILTERED_TOTAL ))

Raw sample count: $RAW_SAMPLES
EOF

echo "Basic summary written to ${REPORT_DIR}/summary_basic.txt"

# Biallelic check on filtered
FILTERED_MULTI=$($BCFTOOLS view -H "$FILTERED_VCF" | awk -F'\t' '{if(gsub(",","",$5)>0) print $1":"$2}' | wc -l)
FILTERED_BIAL=$($BCFTOOLS view -H -m2 -M2 "$FILTERED_VCF" | wc -l)
cat >> "${REPORT_DIR}/summary_basic.txt" <<EOF

Filtered VCF multi-allelic ALT count: $FILTERED_MULTI
Filtered VCF biallelic site count (m2 M2): $FILTERED_BIAL
EOF

echo "Evaluating site-level filters on RAW joint VCF..."

SORT_CMD="LC_ALL=C sort -k1,1V -k2,2n"

# (1) QUAL: missing or < QUAL_MIN
$BCFTOOLS view -H "$RAW_VCF" | awk -F'\t' -v qmin="$QUAL_MIN" '{
  qual=$6;
  if(qual=="."){ print $1"\t"$2 }
  else if((qual+0) < qmin){ print $1"\t"$2 }
}' | eval "$SORT_CMD" > "${REPORT_DIR}/fail_QUAL.list"

# (2) QD (INFO/QD)
$BCFTOOLS query -f '%CHROM\t%POS\t%INFO/QD\n' "$RAW_VCF" | awk -v qdmin="$QD_MIN" '{
  if($3=="." || $3==""){ print $1"\t"$2 }
  else if(($3+0) < qdmin){ print $1"\t"$2 }
}' | eval "$SORT_CMD" > "${REPORT_DIR}/fail_QD.list"

# (3) MQ (INFO/MQ)
$BCFTOOLS query -f '%CHROM\t%POS\t%INFO/MQ\n' "$RAW_VCF" | awk -v mqmin="$MQ_MIN" '{
  if($3=="." || $3==""){ print $1"\t"$2 }
  else if(($3+0) < mqmin){ print $1"\t"$2 }
}' | eval "$SORT_CMD" > "${REPORT_DIR}/fail_MQ.list"

# (4) INFO/DP missing or outside [MIN_DP,MAX_DP_CAP]
$BCFTOOLS query -f '%CHROM\t%POS\t%INFO/DP\n' "$RAW_VCF" | awk -v mind="$MIN_DP" -v maxd="$MAX_DP_CAP" '{
  if($3=="." || $3==""){ print $1"\t"$2 }
  else if(($3+0) < mind || ($3+0) > maxd){ print $1"\t"$2 }
}' | eval "$SORT_CMD" > "${REPORT_DIR}/fail_INFO_DP.list"

# (5) MAF (INFO/AF or AC/AN) - missing counts as fail
$BCFTOOLS query -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\n' "$RAW_VCF" | \
awk -v mafth="$MAF_THRESHOLD" '{
  af=$3; ac=$4; an=$5;
  if(af=="." || af==""){
    if(ac=="" || ac=="." || an=="" || an=="."){ print $1"\t"$2 }
    else {
      split(ac,AC,","); a=AC[1]+0; n=an+0;
      if(n==0){ print $1"\t"$2 }
      else { afcalc=a/n; maf = (afcalc<=0.5?afcalc:1-afcalc); if(maf < mafth) print $1"\t"$2 }
    }
  } else {
    afv=af+0; maf=(afv<=0.5?afv:1-afv);
    if(maf < mafth) print $1"\t"$2;
  }
}' | eval "$SORT_CMD" > "${REPORT_DIR}/fail_MAF.list"

# (6) Call-rate via INFO/AN (AN/(2*N) < MIN_CALL_RATE) - missing counts as fail
NUM_SAMPLES=$RAW_SAMPLES
$BCFTOOLS query -f '%CHROM\t%POS\t%INFO/AN\n' "$RAW_VCF" | awk -v ns="$NUM_SAMPLES" -v callmin="$MIN_CALL_RATE" '{
  an=$3;
  if(an=="" || an=="." || an==0){ print $1"\t"$2 }
  else {
    frac = an/(2*ns);
    if(frac < callmin) print $1"\t"$2
  }
}' | eval "$SORT_CMD" > "${REPORT_DIR}/fail_CALLRATE.list"

# Summarize counts
QUAL_FAIL_COUNT=$(wc -l < "${REPORT_DIR}/fail_QUAL.list")
QD_FAIL_COUNT=$(wc -l < "${REPORT_DIR}/fail_QD.list")
MQ_FAIL_COUNT=$(wc -l < "${REPORT_DIR}/fail_MQ.list")
INFO_DP_FAIL_COUNT=$(wc -l < "${REPORT_DIR}/fail_INFO_DP.list")
MAF_FAIL_COUNT=$(wc -l < "${REPORT_DIR}/fail_MAF.list")
CALLRATE_FAIL_COUNT=$(wc -l < "${REPORT_DIR}/fail_CALLRATE.list")

cat > "${REPORT_DIR}/site_filter_fail_summary.txt" <<EOF
Site-filter failure counts on RAW joint VCF

Total raw sites: $RAW_TOTAL

QUAL fail (missing or < $QUAL_MIN): $QUAL_FAIL_COUNT
QD fail (missing or < $QD_MIN): $QD_FAIL_COUNT
MQ fail (missing or < $MQ_MIN): $MQ_FAIL_COUNT
INFO/DP fail (missing or outside [$MIN_DP,$MAX_DP_CAP]): $INFO_DP_FAIL_COUNT
MAF fail (maf < $MAF_THRESHOLD or AF/AC+AN missing): $MAF_FAIL_COUNT
Call-rate fail (INFO/AN/(2*N) < $MIN_CALL_RATE or AN missing): $CALLRATE_FAIL_COUNT

EOF

echo "Site-level fail summary written to ${REPORT_DIR}/site_filter_fail_summary.txt"

# Compare raw vs filtered
echo "Listing removed sites (present in RAW but absent in FILTERED) ..."

$BCFTOOLS query -f '%CHROM\t%POS\n' "$FILTERED_VCF" > "${REPORT_DIR}/filtered_sites.unsorted.txt"
$BCFTOOLS query -f '%CHROM\t%POS\n' "$RAW_VCF" > "${REPORT_DIR}/raw_sites.unsorted.txt"

awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR { k[$1":"$2]=1; next } { if(!(($1":"$2) in k)) print $1,$2 }' \
    "${REPORT_DIR}/filtered_sites.unsorted.txt" "${REPORT_DIR}/raw_sites.unsorted.txt" \
    > "${REPORT_DIR}/sites_removed_by_filtering.txt"

# sorted version for readability
LC_ALL=C sort -k1,1V -k2,2n "${REPORT_DIR}/sites_removed_by_filtering.txt" -o "${REPORT_DIR}/sites_removed_by_filtering.sorted.txt"

REMOVED_COUNT=$(wc -l < "${REPORT_DIR}/sites_removed_by_filtering.txt")
echo "Sites removed by filtering: $REMOVED_COUNT (list: ${REPORT_DIR}/sites_removed_by_filtering.txt). Sorted copy: sites_removed_by_filtering.sorted.txt"

# For each per-filter list, compute overlap with removed sites
echo "Intersecting removed sites with per-filter fail lists..."
for f in QUAL QD MQ INFO_DP MAF CALLRATE; do
  src="${REPORT_DIR}/fail_${f}.list"
  out="${REPORT_DIR}/removed_by_${f}.list"
  awk 'BEGIN{FS="\t"; OFS="\t"} NR==FNR { rem[$1":"$2]=1; next } { if(($1":"$2) in rem) print $1,$2 }' \
      "${REPORT_DIR}/sites_removed_by_filtering.txt" "$src" > "$out" || true
  LC_ALL=C sort -k1,1V -k2,2n "$out" -o "$out"
  cnt=$(wc -l < "$out" 2>/dev/null || echo 0)
  echo "$f: $cnt removed sites overlap" >> "${REPORT_DIR}/site_filter_fail_summary.txt"
done

echo "Appended per-filter overlap counts to site_filter_fail_summary.txt"

# Genotype-level masking: extract GT/DP/GQ matrices for RAW and FILTERED
echo "Extracting GT/DP/GQ matrices for RAW and FILTERED (for masking check)..."
$BCFTOOLS query -f '%CHROM\t%POS[\t%GT]\n' "$RAW_VCF" | gzip -c > "${REPORT_DIR}/genotypes_raw.GT.tsv.gz"
$BCFTOOLS query -f '%CHROM\t%POS[\t%GT]\n' "$FILTERED_VCF" | gzip -c > "${REPORT_DIR}/genotypes_filtered.GT.tsv.gz"
$BCFTOOLS query -f '%CHROM\t%POS[\t%DP]\n' "$RAW_VCF" | gzip -c > "${REPORT_DIR}/genotypes_raw.DP.tsv.gz"
$BCFTOOLS query -f '%CHROM\t%POS[\t%GQ]\n' "$RAW_VCF" | gzip -c > "${REPORT_DIR}/genotypes_raw.GQ.tsv.gz"

echo "Genotype matrices written to ${REPORT_DIR}"

# Final report pointer
cat > "${REPORT_DIR}/final_report.txt" <<EOF
FILTER VERIFICATION REPORT
See files in ${REPORT_DIR}:
 - summary_basic.txt
 - site_filter_fail_summary.txt
 - fail_*.list (per-filter lists)
 - sites_removed_by_filtering.txt and sites_removed_by_filtering.sorted.txt
 - removed_by_*.list (overlap of removed sites with per-filter fails)
 - genotypes_raw.GT.tsv.gz, genotypes_filtered.GT.tsv.gz
 - genotypes_raw.DP.tsv.gz, genotypes_raw.GQ.tsv.gz
