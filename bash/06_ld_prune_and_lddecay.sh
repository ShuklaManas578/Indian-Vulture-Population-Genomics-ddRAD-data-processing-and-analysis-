#!/usr/bin/env bash

# 16 LD_prune and Decay
# ld_prune_and_lddecay.sh

set -euo pipefail
IFS=$'\n\t'

REL_WINDOW=50      # kb
REL_STEP=5
REL_R2=0.40

STR_WINDOW=100
STR_STEP=10
STR_R2=0.30

ALLOW_EXTRA_CHR="--allow-extra-chr"
PLINK_PREFIX="renamed_data_plink"
OUTDIR="ld_prune_results"
mkdir -p "$OUTDIR"

# pairwise LD params
LD_WINDOW_KB=1000
LD_WINDOW=99999
LD_WINDOW_R2=0

# binning (bp)
BIN_SIZE_BP=1000
MAX_PLOT_DIST_BP=$((LD_WINDOW_KB * 1000))

REL_PREFIX="${OUTDIR}/pruned_relaxed"
STR_PREFIX="${OUTDIR}/pruned_stringent"

ALL_R2="${OUTDIR}/ld_pairs_all.ld"
REL_R2_FILE="${OUTDIR}/ld_pairs_relaxed.ld"
STR_R2_FILE="${OUTDIR}/ld_pairs_stringent.ld"

REL_PRUNE_IN="${OUTDIR}/relaxed.prune.in"
STR_PRUNE_IN="${OUTDIR}/stringent.prune.in"

echo "Starting LD pipeline; outputs -> ${OUTDIR}"
echo

# 1) RELAXED prune
echo "1) RELAXED prune: ${REL_WINDOW}kb ${REL_STEP} ${REL_R2}"
plink --bfile "${PLINK_PREFIX}" ${ALLOW_EXTRA_CHR} \
  --indep-pairwise ${REL_WINDOW} ${REL_STEP} ${REL_R2} \
  --out "${OUTDIR}/relaxed" --threads 1

if [[ -f "${OUTDIR}/relaxed.prune.in" ]]; then
  if [[ "$(realpath "${OUTDIR}/relaxed.prune.in")" != "$(realpath "${REL_PRUNE_IN}")" ]]; then
    cp -f "${OUTDIR}/relaxed.prune.in" "${REL_PRUNE_IN}"
  else
    echo "Relaxed prune file already in target location: ${REL_PRUNE_IN}"
  fi
else
  echo "ERROR: relaxed.prune.in not found" >&2
  exit 1
fi

# 2) STRINGENT prune
echo "2) STRINGENT prune: ${STR_WINDOW}kb ${STR_STEP} ${STR_R2}"
plink --bfile "${PLINK_PREFIX}" ${ALLOW_EXTRA_CHR} \
  --indep-pairwise ${STR_WINDOW} ${STR_STEP} ${STR_R2} \
  --out "${OUTDIR}/stringent" --threads 1

if [[ -f "${OUTDIR}/stringent.prune.in" ]]; then
  if [[ "$(realpath "${OUTDIR}/stringent.prune.in")" != "$(realpath "${STR_PRUNE_IN}")" ]]; then
    cp -f "${OUTDIR}/stringent.prune.in" "${STR_PRUNE_IN}"
  else
    echo "Stringent prune file already in target location: ${STR_PRUNE_IN}"
  fi
else
  echo "ERROR: stringent.prune.in not found" >&2
  exit 1
fi

# 3) create pruned PLINK datasets
echo "3) Creating pruned PLINK datasets..."
plink --bfile "${PLINK_PREFIX}" ${ALLOW_EXTRA_CHR} --extract "${REL_PRUNE_IN}" --make-bed --out "${REL_PREFIX}" --threads 1
plink --bfile "${PLINK_PREFIX}" ${ALLOW_EXTRA_CHR} --extract "${STR_PRUNE_IN}" --make-bed --out "${STR_PREFIX}" --threads 1

# 4) create bgzipped VCFs
echo "4) Writing VCFs (bgzip + tabix)..."
plink --bfile "${REL_PREFIX}" --recode vcf --out "${REL_PREFIX}" ${ALLOW_EXTRA_CHR}
bgzip -c "${REL_PREFIX}.vcf" > "${REL_PREFIX}.vcf.gz"
tabix -p vcf "${REL_PREFIX}.vcf.gz"
rm -f "${REL_PREFIX}.vcf"

plink --bfile "${STR_PREFIX}" --recode vcf --out "${STR_PREFIX}" ${ALLOW_EXTRA_CHR}
bgzip -c "${STR_PREFIX}.vcf" > "${STR_PREFIX}.vcf.gz"
tabix -p vcf "${STR_PREFIX}.vcf.gz"
rm -f "${STR_PREFIX}.vcf"


cp -v "${STR_PREFIX}.vcf.gz" "${OUTDIR}/ld_pruned.vcf.gz" || true
if [[ -f "${STR_PREFIX}.vcf.gz.tbi" ]]; then
  cp -v "${STR_PREFIX}.vcf.gz.tbi" "${OUTDIR}/ld_pruned.vcf.gz.tbi" || true
fi

# 5) counts & intersections
echo "5) marker counts & intersections..."
awk '{print $2}' "${PLINK_PREFIX}.bim" > "${OUTDIR}/all_variant_ids.txt"
awk '{print $2}' "${REL_PREFIX}.bim" > "${OUTDIR}/relaxed_markers.txt"
awk '{print $2}' "${STR_PREFIX}.bim" > "${OUTDIR}/stringent_markers.txt"

ALL_N=$(wc -l < "${OUTDIR}/all_variant_ids.txt")
REL_N=$(wc -l < "${OUTDIR}/relaxed_markers.txt")
STR_N=$(wc -l < "${OUTDIR}/stringent_markers.txt")

sort "${OUTDIR}/relaxed_markers.txt" > "${OUTDIR}/relaxed_markers.sorted"
sort "${OUTDIR}/stringent_markers.txt" > "${OUTDIR}/stringent_markers.sorted"
sort "${OUTDIR}/all_variant_ids.txt" > "${OUTDIR}/all_variant_ids.sorted"

REL_STR_INTER_COUNT=$(comm -12 "${OUTDIR}/relaxed_markers.sorted" "${OUTDIR}/stringent_markers.sorted" | wc -l)
REL_ALL_INTER_COUNT=$(comm -12 "${OUTDIR}/relaxed_markers.sorted" "${OUTDIR}/all_variant_ids.sorted" | wc -l)
STR_ALL_INTER_COUNT=$(comm -12 "${OUTDIR}/stringent_markers.sorted" "${OUTDIR}/all_variant_ids.sorted" | wc -l)

cat > "${OUTDIR}/prune_counts_summary.txt" <<EOF
LD pruning summary
-------------------
Original markers (ALL): ${ALL_N}
Relaxed prune retained: ${REL_N}
Stringent prune retained: ${STR_N}

Intersections:
Relaxed ∩ Stringent: ${REL_STR_INTER_COUNT}
Relaxed ∩ All: ${REL_ALL_INTER_COUNT}
Stringent ∩ All: ${STR_ALL_INTER_COUNT}

Files:
ALL marker list: ${OUTDIR}/all_variant_ids.txt
RELAXED marker list: ${OUTDIR}/relaxed_markers.txt
STRINGENT marker list: ${OUTDIR}/stringent_markers.txt
EOF

cat "${OUTDIR}/prune_counts_summary.txt"
echo

# 6) compute pairwise LD (plink --r2)
echo "6) computing pairwise LD..."

# ALL
plink --bfile "${PLINK_PREFIX}" ${ALLOW_EXTRA_CHR} \
  --r2 --ld-window-kb "${LD_WINDOW_KB}" --ld-window "${LD_WINDOW}" --ld-window-r2 "${LD_WINDOW_R2}" \
  --out "${OUTDIR}/ld_pairs_all" --threads 1
if [[ -f "${OUTDIR}/ld_pairs_all.ld" ]]; then
  mv -f "${OUTDIR}/ld_pairs_all.ld" "${ALL_R2}"
fi

# RELAXED
plink --bfile "${REL_PREFIX}" ${ALLOW_EXTRA_CHR} \
  --r2 --ld-window-kb "${LD_WINDOW_KB}" --ld-window "${LD_WINDOW}" --ld-window-r2 "${LD_WINDOW_R2}" \
  --out "${OUTDIR}/ld_pairs_relaxed" --threads 1
if [[ -f "${OUTDIR}/ld_pairs_relaxed.ld" ]]; then
  mv -f "${OUTDIR}/ld_pairs_relaxed.ld" "${REL_R2_FILE}"
fi

# STRINGENT
plink --bfile "${STR_PREFIX}" ${ALLOW_EXTRA_CHR} \
  --r2 --ld-window-kb "${LD_WINDOW_KB}" --ld-window "${LD_WINDOW}" --ld-window-r2 "${LD_WINDOW_R2}" \
  --out "${OUTDIR}/ld_pairs_stringent" --threads 1
if [[ -f "${OUTDIR}/ld_pairs_stringent.ld" ]]; then
  mv -f "${OUTDIR}/ld_pairs_stringent.ld" "${STR_R2_FILE}"
fi

# 7) prepare R inputs
cat > "${OUTDIR}/ld_plot_inputs.tsv" <<EOF
label	file
ALL	${ALL_R2}
RELAXED	${REL_R2_FILE}
STRINGENT	${STR_R2_FILE}
EOF

cat > "${OUTDIR}/ld_plot_params.sh" <<EOF
OUTDIR='$(realpath "${OUTDIR}")'
BIN_SIZE_BP=${BIN_SIZE_BP}
MAX_PLOT_DIST_BP=${MAX_PLOT_DIST_BP}
DPI=1200
EOF

# 8) write the R plotting script
cat > "${OUTDIR}/plot_ld_decay.R" <<'R_SCRIPT'
#!/usr/bin/env Rscript
# plot_ld_decay.R
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: plot_ld_decay.R <inputs.tsv> <params.sh>")
inputs_tsv <- args[1]
params_sh <- args[2]
source(params_sh)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(viridis)
})

inp <- read.table(inputs_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
all_results <- list()
breaks <- seq(0, MAX_PLOT_DIST_BP, by = BIN_SIZE_BP)

for (i in seq_len(nrow(inp))) {
  label <- inp$label[i]
  file <- inp$file[i]
  if (!file.exists(file)) {
    warning("LD file not found: ", file); next
  }
  dat <- tryCatch(read.table(file, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(dat)) next
  cn <- names(dat)
  # detect BP columns and R2
  bpAcol <- grep("BP", cn, ignore.case = TRUE, value = TRUE)[1]
  bpBcol <- grep("BP", cn, ignore.case = TRUE, value = TRUE)[2]
  r2col  <- grep("^R2$|^R.2$|r2|R2\\.", cn, ignore.case = TRUE, value = TRUE)[1]
  if (is.na(bpAcol) || is.na(bpBcol) || is.na(r2col)) {
    if (ncol(dat) >= 7) {
      bpAcol <- cn[2]; bpBcol <- cn[5]; r2col <- cn[7]
    } else stop("Unable to parse LD file columns in ", file)
  }
  dat$BP_A <- as.numeric(dat[[bpAcol]])
  dat$BP_B <- as.numeric(dat[[bpBcol]])
  dat$R2   <- as.numeric(dat[[r2col]])
  dat <- dat %>% filter(!is.na(BP_A) & !is.na(BP_B) & !is.na(R2))
  dat$dist_bp <- abs(dat$BP_A - dat$BP_B)
  dat <- dat %>% filter(dist_bp > 0 & dist_bp <= MAX_PLOT_DIST_BP)
  if (nrow(dat) == 0) next
  idx <- findInterval(dat$dist_bp, breaks, rightmost.closed = TRUE)
  valid <- which(idx > 0 & idx < length(breaks))
  if (length(valid) == 0) next
  dat2 <- dat[valid, , drop = FALSE]
  idx2 <- idx[valid]
  mids <- (breaks[idx2] + breaks[pmin(idx2+1, length(breaks))]) / 2
  summary_by_bin <- data.frame(bin_mid = tapply(mids, idx2, mean),
                               mean_r2 = tapply(dat2$R2, idx2, mean),
                               n_pairs = tapply(dat2$R2, idx2, length),
                               stringsAsFactors = FALSE)
  summary_by_bin$label <- label
  all_results[[label]] <- summary_by_bin
}

if (length(all_results) == 0) stop("No LD data to plot.")
df <- bind_rows(all_results)
df$dist_kb <- df$bin_mid / 1000

# annotation from prune counts
prune_file <- file.path(OUTDIR, "prune_counts_summary.txt")
anno_text <- NULL
if (file.exists(prune_file)) {
  txt <- readLines(prune_file)
  lines <- grep("Original markers|Relaxed prune retained|Stringent prune retained", txt, value = TRUE)
  anno_text <- paste(lines, collapse = "\n")
}

palette_lines <- viridis::viridis(n = length(unique(df$label)))
p <- ggplot(df, aes(x = dist_kb, y = mean_r2, color = label)) +
  geom_line(size = 1) + geom_point(size = 0.8, alpha = 0.9) +
  labs(x = "Distance (kb)", y = expression(mean~r^2),
       title = "LD decay (mean r^2 by distance bin)", color = "Dataset") +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = palette_lines) +
  theme(plot.title = element_text(face = "bold"))

if (!is.null(anno_text)) {
  p <- p + annotate("text", x = Inf, y = Inf, label = anno_text, hjust = 1.02, vjust = 1.2, size = 3.6)
}

outfile <- file.path(OUTDIR, "ld_decay_mean_r2.png")
ggsave(outfile, p, width = 10, height = 6, dpi = DPI)
cat("Wrote LD decay plot to:", outfile, "\n")
R_SCRIPT

echo "9) Running R plotting script..."
Rscript "${OUTDIR}/plot_ld_decay.R" "${OUTDIR}/ld_plot_inputs.tsv" "${OUTDIR}/ld_plot_params.sh" || {
  echo "Warning: R plotting failed; check R packages and input files." >&2
}

echo "All done. Key outputs in ${OUTDIR}:"
ls -1 "${OUTDIR}" | sed -n '1,200p'
echo "Done."
