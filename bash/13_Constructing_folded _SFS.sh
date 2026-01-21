#!/usr/bin/env bash


# 15 Constructing folded SFS
set -euo pipefail
# foldedSFS_angsd.sh  -- robust folded SFS pipeline using ANGSD + realSFS

# USER-CONFIGURATION
THREADS=8
MINMAPQ=30
MINQ=20
GL=1                       # 2 = GATK genotype likelihood model
minIndBKN=5                 
minIndHYD=4              
REF_FA="/home/manas/Downloads/Vultures_popGen_analysis/results_modular/gatk_pipeline_20250924_2254/ref.fa"
OUTDIR="angsd_folded_results_GL${GL}_gatk"
COMMON_OPTS="-minMapQ ${MINMAPQ} -minQ ${MINQ} -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -baq 1 -C 50"
REN_DIR="renamed_bams"


echo "[$(date +\"%A %d %B %Y %r %Z\")] Starting folded-SFS pipeline -> ${OUTDIR}"

# check required programs
for prog in angsd realSFS samtools python3; do
  if ! command -v ${prog} >/dev/null 2>&1; then
    echo "ERROR: ${prog} not found on PATH. Activate angsd_env and ensure samtools & python3 in PATH." >&2
    exit 1
  fi
done

# check REF exists
if [[ ! -f "${REF_FA}" ]]; then
  echo "ERROR: REF fasta not found at: ${REF_FA}" >&2
  exit 1
fi

# ensure fasta index exists
if [[ ! -f "${REF_FA}.fai" ]]; then
  echo "[$(date)] Creating fasta index (${REF_FA}.fai) with samtools faidx..."
  samtools faidx "${REF_FA}"
fi

# prepare directories
mkdir -p "${OUTDIR}"
mkdir -p "${REN_DIR}"
rm -f "${REN_DIR}"/* || true

# Create renamed symlinks (WR -> HYD, EV/EG -> BKN) and skip EG_35
echo "[$(date)] Creating symlinks (WR->HYD, EV/EG->BKN) excluding EG_35..."
shopt -s nullglob
for f in *.bam; do
  base="$(basename "$f")"
  if [[ "$base" == EG_35* ]]; then
    echo "Skipping excluded sample: $base"
    continue
  fi
  if [[ "$base" == WR_* ]]; then
    new="HYD_${base#WR_}"
  elif [[ "$base" == EV_* || "$base" == EG_* ]]; then
    new="BKN_${base#*_}"
  else
    echo "Skipping non-matching file: $base"
    continue
  fi
  ln -sf "$(pwd)/${base}" "${REN_DIR}/${new}"
done

# ensure .bai for each renamed bam: link if found, otherwise create index
pushd "${REN_DIR}" >/dev/null
for f in *.bam; do
  if [[ ! -f "${f}.bai" ]]; then
    suffix="${f#*_}"
    found=0
    for cand in ../*"${suffix}"; do
      if [[ -f "${cand}.bai" ]]; then
        ln -sf "${cand}.bai" "${f}.bai"
        found=1
        break
      fi
    done
    if [[ $found -eq 0 ]]; then
      echo "[$(date)] Indexing ${f} -> creating ${f}.bai with samtools index"
      samtools index "${f}"
    fi
  fi
done
popd >/dev/null

# write absolute-path bam lists
BAMLIST_BKN="${OUTDIR}/bamlist_BKN.txt"
BAMLIST_HYD="${OUTDIR}/bamlist_HYD.txt"
ls -1 "$(pwd)/${REN_DIR}"/BKN_*.bam 2>/dev/null | sort > "${BAMLIST_BKN}" || true
ls -1 "$(pwd)/${REN_DIR}"/HYD_*.bam 2>/dev/null | sort > "${BAMLIST_HYD}" || true

if [[ ! -s "${BAMLIST_BKN}" || ! -s "${BAMLIST_HYD}" ]]; then
  echo "ERROR: One of the bamlist files is empty. Check ${REN_DIR} for symlinks." >&2
  ls -l "${REN_DIR}" || true
  exit 1
fi

nBKN=$(wc -l < "${BAMLIST_BKN}")
nHYD=$(wc -l < "${BAMLIST_HYD}")
echo "[$(date)] Found nBKN=${nBKN} nHYD=${nHYD}"

# IMPORTANT: some ANGSD builds require -anc when running -doSaf.
# supply the reference as -anc to allow ANGSD to produce .saf.idx;
ANC_OPT="-anc ${REF_FA}"

OUT_BKN="${OUTDIR}/BKN"
OUT_HYD="${OUTDIR}/HYD"

echo "[$(date)] Running ANGSD (BKN)..."
angsd -b "${BAMLIST_BKN}" -GL ${GL} -doSaf 1 ${ANC_OPT} -ref "${REF_FA}" ${COMMON_OPTS} -out "${OUT_BKN}" -minInd ${minIndBKN} -P ${THREADS} > "${OUTDIR}/angsd_BKN.log" 2>&1 || {
  echo "ERROR: angsd failed for BKN. See ${OUTDIR}/angsd_BKN.log" >&2
  tail -n 200 "${OUTDIR}/angsd_BKN.log"
  exit 1
}
if [[ ! -s "${OUT_BKN}.saf.idx" ]]; then
  echo "ERROR: ANGSD did not produce ${OUT_BKN}.saf.idx. See ${OUTDIR}/angsd_BKN.log" >&2
  tail -n 200 "${OUTDIR}/angsd_BKN.log"
  exit 1
fi

echo "[$(date)] Running ANGSD (HYD)..."
angsd -b "${BAMLIST_HYD}" -GL ${GL} -doSaf 1 ${ANC_OPT} -ref "${REF_FA}" ${COMMON_OPTS} -out "${OUT_HYD}" -minInd ${minIndHYD} -P ${THREADS} > "${OUTDIR}/angsd_HYD.log" 2>&1 || {
  echo "ERROR: angsd failed for HYD. See ${OUTDIR}/angsd_HYD.log" >&2
  tail -n 200 "${OUTDIR}/angsd_HYD.log"
  exit 1
}
if [[ ! -s "${OUT_HYD}.saf.idx" ]]; then
  echo "ERROR: ANGSD did not produce ${OUT_HYD}.saf.idx. See ${OUTDIR}/angsd_HYD.log" >&2
  tail -n 200 "${OUTDIR}/angsd_HYD.log"
  exit 1
fi

# realSFS 1D and fold
echo "[$(date)] Running realSFS for 1D..."
realSFS "${OUT_BKN}.saf.idx" -P ${THREADS} > "${OUTDIR}/BKN.raw.sfs"
realSFS "${OUT_HYD}.saf.idx" -P ${THREADS} > "${OUTDIR}/HYD.raw.sfs"
echo "[$(date)] Folding 1D SFSs..."
realSFS fold "${OUTDIR}/BKN.raw.sfs" > "${OUTDIR}/BKN.folded.sfs"
realSFS fold "${OUTDIR}/HYD.raw.sfs" > "${OUTDIR}/HYD.folded.sfs"

# realSFS 2D joint then fold
echo "[$(date)] Running realSFS for 2D joint..."
realSFS "${OUT_BKN}.saf.idx" "${OUT_HYD}.saf.idx" -P ${THREADS} > "${OUTDIR}/BKN_HYD.raw.2dsfs"
echo "[$(date)] Folding joint 2D..."
realSFS fold "${OUTDIR}/BKN_HYD.raw.2dsfs" > "${OUTDIR}/BKN_HYD.folded.2dsfs"

# convert flattened folded 2D to CSV matrix
catsBKN=$(( nBKN * 2 + 1 ))
catsHYD=$(( nHYD * 2 + 1 ))

python3 - <<PY
import numpy as np, sys, os
fn = "${OUTDIR}/BKN_HYD.folded.2dsfs"
if not os.path.exists(fn):
    print("ERROR: expected 2D file not found:", fn)
    sys.exit(1)
with open(fn) as fh:
    toks = fh.read().strip().split()
vec = np.array([float(x) for x in toks]) if toks else np.array([])
expected = ${catsBKN} * ${catsHYD}
if vec.size != expected:
    print("ERROR: unexpected length of 2D file: got", vec.size, "expected", expected)
    sys.exit(1)
mat = vec.reshape((${catsBKN}, ${catsHYD}))
np.savetxt("${OUTDIR}/BKN_HYD.folded.2dsfs.matrix.csv", mat, delimiter=",", fmt="%.6f")
print("WROTE:", "${OUTDIR}/BKN_HYD.folded.2dsfs.matrix.csv")
PY

echo "[$(date)] Folded SFS pipeline finished. Output directory: ${OUTDIR}"
ls -lh "${OUTDIR}" | sed -n '1,200p'
