#!/usr/bin/env Rscript


#6 Heterozygosity plots (Rscript)

# Reads a PLINK .het and the matching .imiss file (if available), computes/uses missingness,


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(viridis)
})

het_file <- "/home/manas/Downloads/Vultures_popGen_analysis/raw_data/mapped_bam/gvcfs/joint_out/rename_qc_allvariants/renamed_data_plink_het_allvariants.het"
imiss_file <- "/home/manas/Downloads/Vultures_popGen_analysis/raw_data/mapped_bam/gvcfs/joint_out/rename_qc_allvariants/renamed_data_plink_missing.imiss"
out_base <- file.path(dirname(het_file), "plots", "combined")
dpi <- 1200
label_f_thresh <- 0.1
label_miss_thresh <- 0.05


dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

safe_read_table <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  res <- tryCatch({
    read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "", fill = TRUE)
  }, error = function(e) {
    readr::read_table2(path, col_names = TRUE)
  })
  # normalize column names
  colnames(res) <- trimws(colnames(res))
  return(res)
}

find_col_name <- function(cols, patterns) {
  cols_trim <- trimws(cols)
  for (p in patterns) {
    idx <- grep(p, cols_trim, ignore.case = TRUE)
    if (length(idx) > 0) return(cols[idx[1]])
  }
  return(NA_character_)
}

safe_num <- function(x) {
  if (is.null(x)) return(NA_real_)
  out <- suppressWarnings(as.numeric(as.character(x)))
  return(out)
}

all_na <- function(vec) {
  if (length(vec) == 0) return(TRUE)
  all(is.na(vec))
}

save_figure <- function(plot_obj, png_name, width = 10, height = 7, dpi_val = dpi) {
  outpng <- file.path(out_base, png_name)
  outsvg <- sub("\\.png$", ".svg", outpng)
  ggsave(filename = outpng, plot = plot_obj, width = width, height = height, dpi = dpi_val)
  # write SVG without forcing dpi
  ggsave(filename = outsvg, plot = plot_obj, width = width, height = height)
  message("Saved: ", basename(outpng), " and ", basename(outsvg))
}

# Read .het
if (!file.exists(het_file)) stop("het file not found: ", het_file)
message("[INFO] Reading .het file: ", het_file)
het_raw <- safe_read_table(het_file)
orig_cols <- colnames(het_raw)

iid_col <- find_col_name(orig_cols, c("^IID$","^i$","^id$","^sample$","^s$"))
fid_col <- find_col_name(orig_cols, c("^FID$","^family$","^f$","^family.id$"))
o_col <- find_col_name(orig_cols, c("^O\\W*HOM$","^O\\(HOM\\)$","^OBS\\.HOM$","^OBS_HOM$","^O$","^OBS$"))
e_col <- find_col_name(orig_cols, c("^E\\W*HOM$","^E\\(HOM\\)$","^EXP\\.HOM$","^EXP_HOM$","^E$","^EXP$"))
n_col <- find_col_name(orig_cols, c("^OBS_CT$","^OBS\\.CT$","^N\\(.+\\)$","^N\\(NM\\)$","^N$","^COUNT$","^NOBS$","^NCOUNT$"))
f_col <- find_col_name(orig_cols, c("^F$","^F\\w*HAT$","^INBREEDING$","^FPLINK$","^F_HAT$"))
pop_col <- find_col_name(orig_cols, c("^POP$","^POPULATION$","^GROUP$","^POP_ID$"))

if (is.na(iid_col)) {
  if (ncol(het_raw) >= 2) iid_col <- orig_cols[2] else iid_col <- orig_cols[1]
  message("[WARN] Could not find IID by name; falling back to column: ", iid_col)
}

df <- het_raw
df$sample <- trimws(as.character(df[[iid_col]]))
df$fid <- if (!is.na(fid_col)) trimws(as.character(df[[fid_col]])) else NA_character_
df$O_HOM <- if (!is.na(o_col)) safe_num(df[[o_col]]) else NA_real_
df$E_HOM <- if (!is.na(e_col)) safe_num(df[[e_col]]) else NA_real_
df$N_COUNT <- if (!is.na(n_col)) safe_num(df[[n_col]]) else NA_real_
df$F_plink <- if (!is.na(f_col)) safe_num(df[[f_col]]) else NA_real_

# missingness: prefer imiss file if present, otherwise try columns in het or fallback compute
df$missing_prop <- NA_real_

miss_col_in_het <- find_col_name(orig_cols, c("^MISSING$","^MISSING_PROP$","^MISSINGNESS$","^MISS$","^F_MISS$"))
if (!is.na(miss_col_in_het)) {
  df$missing_prop <- safe_num(df[[miss_col_in_het]])
  message("[INFO] Using missingness column in .het: ", miss_col_in_het)
}

# Read imiss if exists and use it
if (file.exists(imiss_file)) {
  message("[INFO] Reading imiss file: ", imiss_file)
  im <- safe_read_table(imiss_file)
  im_cols <- colnames(im)
  iid_im_col <- find_col_name(im_cols, c("^IID$","^i$","^id$","^sample$"))
  fmiss_col <- find_col_name(im_cols, c("^F_MISS$","^F_MISS\\b","^F_MISS$","^F_MISS"))
  nmiss_col <- find_col_name(im_cols, c("^N_MISS$","^N_MISS\\b","^N_MISS"))
  ngen_col <- find_col_name(im_cols, c("^N_GENO$","^N_GENO","^N_GENO\\b"))
  if (is.na(iid_im_col)) iid_im_col <- im_cols[1]
  
  im <- im %>% mutate(iid_im = as.character(.data[[iid_im_col]]))
  
  # if F_MISS is present, use as missing_prop; otherwise compute from N_MISS / N_GENO
  if (!is.na(fmiss_col)) {
    im <- im %>% mutate(missing_prop_im = safe_num(.data[[fmiss_col]]))
  } else if (!is.na(nmiss_col) && !is.na(ngen_col)) {
    im <- im %>% mutate(missing_prop_im = safe_num(.data[[nmiss_col]]) / safe_num(.data[[ngen_col]]))
  } else if (!is.na(nmiss_col)) {
    im <- im %>% mutate(missing_prop_im = safe_num(.data[[nmiss_col]]))
  } else {
    im$missing_prop_im <- NA_real_
  }
  
  # join to main df by sample/IID
  df <- df %>% left_join(im %>% select(iid_im, missing_prop_im), by = c("sample" = "iid_im"))
  be_used <- !is.na(df$missing_prop_im)
  if (any(be_used)) {
    df$missing_prop[be_used] <- df$missing_prop_im[be_used]
    message("[INFO] Applied missing_prop from imiss file for ", sum(be_used), " samples")
  } else {
    message("[INFO] imiss file read but no usable missingness found")
  }
} else {
  message("[INFO] imiss file not found at path: ", imiss_file, "  -> will attempt fallback missingness.")
}

# Fallback
if (all(is.na(df$missing_prop)) && !all(is.na(df$N_COUNT))) {
  maxNc <- max(df$N_COUNT, na.rm = TRUE)
  if (is.finite(maxNc) && maxNc > 0) {
    df$missing_prop <- 1 - (df$N_COUNT / maxNc)
    message("[INFO] Computed missing_prop fallback from N_COUNT for all samples")
  }
}

# If some samples still missing missing_prop, compute per-sample fallback using available N_COUNT
na_idx <- which(is.na(df$missing_prop) & !is.na(df$N_COUNT))
if (length(na_idx) > 0 && is.finite(max(df$N_COUNT, na.rm = TRUE))) {
  maxNc <- max(df$N_COUNT, na.rm = TRUE)
  df$missing_prop[na_idx] <- 1 - (df$N_COUNT[na_idx] / maxNc)
  message("[INFO] Computed missing_prop for ", length(na_idx), " additional samples using N_COUNT fallback")
}

if (!is.na(pop_col) && pop_col %in% colnames(het_raw)) {
  df$pop_raw <- as.character(het_raw[[pop_col]])
} else {
  df$pop_raw <- NA_character_
  df$pop_raw[grepl("^(EG|EG_|EG-)", df$sample, ignore.case = TRUE)] <- "EG"
  df$pop_raw[grepl("^(BKN|BKN_)", df$sample, ignore.case = TRUE)] <- "BKN"
  df$pop_raw[grepl("^(HYD|HYD_)", df$sample, ignore.case = TRUE)] <- "HYD"
  na_idx <- which(is.na(df$pop_raw))
  if (length(na_idx) > 0) {
    df$pop_raw[na_idx] <- sapply(df$sample[na_idx], function(s) {
      t <- strsplit(s, "_|-")[[1]][1]
      if (t %in% c("EG","BKN","HYD")) return(t)
      return(NA_character_)
    })
  }
}

safe_div <- function(a,b) ifelse(is.na(a) | is.na(b) | b == 0, NA_real_, a/b)
df <- df %>% mutate(
  prop_obs_hom = safe_div(O_HOM, N_COUNT),
  prop_exp_hom = safe_div(E_HOM, N_COUNT),
  obs_het = ifelse(!is.na(prop_obs_hom), 1 - prop_obs_hom, NA_real_),
  exp_het = ifelse(!is.na(prop_exp_hom), 1 - prop_exp_hom, NA_real_),
  F_calc = ifelse(!is.na(F_plink), F_plink,
                  ifelse(!is.na(O_HOM) & !is.na(E_HOM) & !is.na(N_COUNT) & (N_COUNT - E_HOM) != 0,
                         (O_HOM - E_HOM) / (N_COUNT - E_HOM), NA_real_)),
  pop = pop_raw
)

# Ensure unique by sample
df$sample <- trimws(as.character(df$sample))
df <- df %>% distinct(sample, .keep_all = TRUE)

# Main filtered dataset: BKN & HYD only, excluding EG_35 (for all main plots)
df_filtered <- df %>% filter(sample != "EG_35", pop %in% c("BKN", "HYD"))

df_remap <- df

df_remap$pop <- as.character(df_remap$pop)
idx_eg35 <- which(df_remap$sample == "EG_35")
if (length(idx_eg35) > 0) {
  df_remap$sample[idx_eg35] <- "BKN_UV_o_01"
  df_remap$pop[idx_eg35] <- "BKN_o"
  message("[INFO] Remapped EG_35 -> BKN_UV_o_01 with pop = BKN_o for remapped dataset")
}
df_remap$pop <- factor(df_remap$pop, levels = unique(na.omit(df_remap$pop)))

# -------------- Colors --------------
col_BKN <- "#008B8B"
col_HYD <- "#DE3163"
col_BKN_o <- "#B5A642" 

pal_filtered <- c("BKN" = col_BKN, "HYD" = col_HYD)
pal_with_bkno <- pal_filtered
pal_with_bkno["BKN_o"] <- col_BKN_o

extra_filt <- setdiff(unique(na.omit(df_filtered$pop)), names(pal_filtered))
if (length(extra_filt) > 0) {
  addcols <- viridis::cividis(length(extra_filt))
  names(addcols) <- extra_filt
  pal_filtered <- c(pal_filtered, addcols)
}
extra_remap <- setdiff(unique(na.omit(df_remap$pop)), names(pal_with_bkno))
if (length(extra_remap) > 0) {
  addcols2 <- viridis::cividis(length(extra_remap))
  names(addcols2) <- extra_remap
  pal_with_bkno <- c(pal_with_bkno, addcols2)
}

write.csv(df, file.path(out_base, "combined_heterozygosity_per_sample_full_input.csv"), row.names = FALSE)
write.csv(df_filtered, file.path(out_base, "combined_heterozygosity_per_sample_BKN_HYD_only.csv"), row.names = FALSE)
message("[INFO] Wrote combined CSVs to: ", out_base)

# Per-pop summary for BKN & HYD
pop_summary_filtered <- df_filtered %>%
  group_by(pop) %>%
  summarise(
    n = n(),
    mean_prop_obs_hom = mean(prop_obs_hom, na.rm = TRUE),
    mean_prop_exp_hom = mean(prop_exp_hom, na.rm = TRUE),
    mean_obs_het = mean(obs_het, na.rm = TRUE),
    mean_exp_het = mean(exp_het, na.rm = TRUE),
    mean_F = mean(F_calc, na.rm = TRUE),
    mean_missing = mean(missing_prop, na.rm = TRUE)
  ) %>% arrange(pop)

write.csv(pop_summary_filtered, file.path(out_base, "combined_heterozygosity_summary_by_pop_BKN_HYD.csv"), row.names = FALSE)
message("[INFO] Wrote per-pop summary (BKN & HYD): ", file.path(out_base, "combined_heterozygosity_summary_by_pop_BKN_HYD.csv"))

# Produce main plots (BKN & HYD only)
# Observed homozygosity
if (!all_na(df_filtered$prop_obs_hom) && nrow(df_filtered) > 0) {
  p_obs_hom <- ggplot(df_filtered, aes(x = pop, y = prop_obs_hom, fill = pop)) +
    geom_violin(width = 0.9, alpha = 0.25, trim = TRUE, na.rm = TRUE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.8, na.rm = TRUE) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 3.5, color = "black", na.rm = TRUE) +
    scale_fill_manual(values = pal_filtered) +
    labs(x = "Population", y = "Observed homozygosity (proportion)", title = "Observed homozygosity") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  save_figure(p_obs_hom, "combined_box_obs_homozygosity_BKN_HYD_only.png")
} else message("[SKIP] Observed homozygosity (no usable data for BKN/HYD)")

# Expected homozygosity
if (!all_na(df_filtered$prop_exp_hom) && nrow(df_filtered) > 0) {
  p_exp_hom <- ggplot(df_filtered, aes(x = pop, y = prop_exp_hom, fill = pop)) +
    geom_violin(width = 0.9, alpha = 0.25, trim = TRUE, na.rm = TRUE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.8, na.rm = TRUE) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 3.5, color = "black", na.rm = TRUE) +
    scale_fill_manual(values = pal_filtered) +
    labs(x = "Population", y = "Expected homozygosity (proportion)", title = "Expected homozygosity") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  save_figure(p_exp_hom, "combined_box_exp_homozygosity_BKN_HYD_only.png")
} else message("[SKIP] Expected homozygosity (no usable data for BKN/HYD)")

# Observed heterozygosity
if (!all_na(df_filtered$obs_het) && nrow(df_filtered) > 0) {
  p_obs_het <- ggplot(df_filtered, aes(x = pop, y = obs_het, fill = pop)) +
    geom_violin(width = 0.9, alpha = 0.25, trim = TRUE, na.rm = TRUE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.8, na.rm = TRUE) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 3.5, color = "black", na.rm = TRUE) +
    scale_fill_manual(values = pal_filtered) +
    labs(x = "Population", y = "Observed heterozygosity (per-sample)", title = "Observed heterozygosity") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  save_figure(p_obs_het, "combined_box_obs_heterozygosity_BKN_HYD_only.png")
} else message("[SKIP] Observed heterozygosity (no usable data for BKN/HYD)")

# Expected heterozygosity
if (!all_na(df_filtered$exp_het) && nrow(df_filtered) > 0) {
  p_exp_het <- ggplot(df_filtered, aes(x = pop, y = exp_het, fill = pop)) +
    geom_violin(width = 0.9, alpha = 0.25, trim = TRUE, na.rm = TRUE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.8, na.rm = TRUE) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean, geom = "point", shape = 3, size = 3.5, color = "black", na.rm = TRUE) +
    scale_fill_manual(values = pal_filtered) +
    labs(x = "Population", y = "Expected heterozygosity (per-sample)", title = "Expected heterozygosity") +
    theme_minimal(base_size = 12) + theme(legend.position = "none")
  save_figure(p_exp_het, "combined_box_exp_heterozygosity_BKN_HYD_only.png")
} else message("[SKIP] Expected heterozygosity (no usable data for BKN/HYD)")

# TWO F vs missingness scatter plots
# Plot A: BKN & HYD only (EG_35 excluded)
df_bkn_hyd <- df %>% filter(pop %in% c("BKN", "HYD"), sample != "EG_35")
if (nrow(df_bkn_hyd) > 0 && !all_na(df_bkn_hyd$F_calc)) {
  # If missing_prop is NA for some rows, keep them out of scatter (geom_point drops NA)
  df_bkn_hyd <- df_bkn_hyd %>% mutate(label_flag = (abs(F_calc) > label_f_thresh) | (!is.na(missing_prop) & missing_prop > label_miss_thresh))
  p_f_both <- ggplot(df_bkn_hyd, aes(x = missing_prop, y = F_calc, color = pop, label = sample)) +
    geom_point(size = 2.8, alpha = 0.9, na.rm = TRUE) +
    scale_color_manual(values = pal_filtered, na.value = "#999999") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    theme_minimal(base_size = 12) +
    labs(x = "Per-sample missingness (proportion)", y = "Inbreeding coefficient F", title = "F vs missingness") +
    theme(legend.position = "right")
  if (any(df_bkn_hyd$label_flag, na.rm = TRUE)) {
    p_f_both <- p_f_both + geom_text_repel(data = df_bkn_hyd %>% filter(label_flag), size = 3.2, max.overlaps = 200)
  }
  save_figure(p_f_both, "combined_scatter_F_vs_missing_BKN_HYD_only.png")
} else {
  message("[SKIP] F vs missingness (BKN & HYD) -> no usable rows or all F missing")
}

# Plot B: include BKN, HYD and BKN_o (remapped EG_35 -> BKN_UV_o_01)
df_with_bkno <- df_remap %>% filter(pop %in% c("BKN", "HYD", "BKN_o"))
if (nrow(df_with_bkno) > 0 && !all_na(df_with_bkno$F_calc)) {
  df_with_bkno <- df_with_bkno %>% mutate(label_flag = (abs(F_calc) > label_f_thresh) | (!is.na(missing_prop) & missing_prop > label_miss_thresh))
  p_f_bkno <- ggplot(df_with_bkno, aes(x = missing_prop, y = F_calc, color = pop, label = sample)) +
    geom_point(size = 2.8, alpha = 0.9, na.rm = TRUE) +
    scale_color_manual(values = pal_with_bkno, na.value = "#999999") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    theme_minimal(base_size = 12) +
    labs(x = "Per-sample missingness (proportion)", y = "Inbreeding coefficient F", title = "F vs missingness") +
    theme(legend.position = "right")
  if (any(df_with_bkno$label_flag, na.rm = TRUE)) {
    p_f_bkno <- p_f_bkno + geom_text_repel(data = df_with_bkno %>% filter(label_flag), size = 3.2, max.overlaps = 200)
  }
  save_figure(p_f_bkno, "combined_scatter_F_vs_missing_BKN_HYD_BKN_o_included.png")
} else {
  message("[SKIP] F vs missingness (BKN+HYD+BKN_o) -> no usable rows or all F missing")
}

# Density of F (BKN & HYD only)
if (!all_na(df_filtered$F_calc) && nrow(df_filtered) > 0) {
  p_density <- ggplot(df_filtered, aes(x = F_calc, fill = pop)) +
    geom_density(alpha = 0.45, na.rm = TRUE) +
    scale_fill_manual(values = pal_filtered) +
    theme_minimal(base_size = 12) +
    labs(title = "Density of inbreeding coefficient F", x = "F", y = "Density")
  save_figure(p_density, "combined_density_F_by_pop_BKN_HYD_only.png")
} else message("[SKIP] Density of F -> no usable BKN/HYD data")


summary_txt <- file.path(out_base, "combined_heterozygosity_text_summary_BKN_HYD_only.txt")
sink(summary_txt)
cat("Combined heterozygosity summary (BKN & HYD only)\n")
cat("====\n\n")
cat("Input .het:", het_file, "\n")
cat("Input imiss (if used):", ifelse(file.exists(imiss_file), imiss_file, "NOT FOUND"), "\n\n")
print(pop_summary_filtered)
cat("\nSamples with extreme F or missingness (thresholds):\n")
extremes <- df_filtered %>% mutate(label_flag = (abs(F_calc) > label_f_thresh) | (!is.na(missing_prop) & missing_prop > label_miss_thresh)) %>% filter(label_flag)
print(extremes %>% select(sample, pop, F_calc, missing_prop, obs_het, exp_het))
cat("\nNotes:\n - EG_35 was excluded from the main BKN/HYD-only plots.\n - For the second F vs missingness plot EG_35 was remapped to BKN_UV_o_01 with pop = BKN_o and included in that plot.\n")
sink()
message("[INFO] Wrote text summary: ", summary_txt)

message("[DONE] All outputs written to: ", out_base)
