#!/usr/bin/env Rscript

# PCA and Plottig (Rscript)
# pca_from_pruned_sets.R
# PCA (SNPRelate) from pruned VCFs; specific outlier labeling.

suppressPackageStartupMessages({
  library(optparse)
  library(SNPRelate)
  library(gdsfmt)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(RColorBrewer)
  library(ggrepel)
})

option_list <- list(
  make_option(c("--pruned_dir"), type = "character",
              default = "Path/to/ld_prune",
              help = "Directory containing pruned_relaxed.vcf.gz and pruned_stringent.vcf.gz"),
  make_option(c("--popmap"), type = "character",
              default = "Path/to/popmap.txt",
              help = "popmap file (sample<TAB>pop)"),
  make_option(c("--outdir"), type = "character", default = NULL,
              help = "output root (default: <pruned_dir>/pca_results)"),
  make_option(c("--dpi"), type = "integer", default = 1200, help = "PNG DPI"),
  make_option(c("--png_w"), type = "numeric", default = 10, help = "PNG width (in)"),
  make_option(c("--png_h"), type = "numeric", default = 8, help = "PNG height (in)"),
  make_option(c("--keep_gds"), type = "logical", default = FALSE, help = "Keep intermediate GDS files")
)
opt <- parse_args(OptionParser(option_list = option_list))

PRUNED_DIR <- normalizePath(opt$pruned_dir)
POPMAP <- opt$popmap
OUTROOT <- if (is.null(opt$outdir)) file.path(PRUNED_DIR, "pca_results") else normalizePath(opt$outdir)
DPI <- opt$dpi; PNG_W <- opt$png_w; PNG_H <- opt$png_h; KEEP_GDS <- opt$keep_gds

dir.create(OUTROOT, recursive = TRUE, showWarnings = FALSE)
message("[INFO] pruned_dir: ", PRUNED_DIR)
message("[INFO] popmap: ", POPMAP)
message("[INFO] outdir: ", OUTROOT)

# expected file names
sets <- c("pruned_relaxed.vcf.gz", "pruned_stringent.vcf.gz")
set_paths <- file.path(PRUNED_DIR, sets)
names(set_paths) <- c("relaxed", "stringent")

ORIG_EG35_NAME <- "EG_35"
RENAMED_EG35 <- "BKN_UV_o_01"
RENAMED_EG35_POP <- "BKN_o"
# -------------------------------

# load popmap
popmap_df <- NULL
if (file.exists(POPMAP)) {
  pm <- tryCatch(read.table(POPMAP, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
                 error = function(e) {
                   message("[ERROR] Failed to read popmap: ", e$message)
                   NULL
                 })
  if (!is.null(pm) && ncol(pm) >= 2) {
    colnames(pm)[1:2] <- c("sample", "pop")
    pm$sample <- trimws(pm$sample)
    pm$pop <- trimws(pm$pop)
    
    pm$sample[pm$sample == ORIG_EG35_NAME] <- RENAMED_EG35
    pm$pop[pm$sample == RENAMED_EG35] <- RENAMED_EG35_POP
    
    if (!any(pm$sample == RENAMED_EG35)) {
      pm <- bind_rows(pm, tibble(sample = RENAMED_EG35, pop = RENAMED_EG35_POP))
      message("[INFO] Added mapping to popmap for renamed sample: ", RENAMED_EG35, " -> ", RENAMED_EG35_POP)
    } else {
      message("[INFO] Ensured popmap contains renamed sample mapping for EG_35 -> ", RENAMED_EG35)
    }
    popmap_df <- pm[, c("sample", "pop")]
    message("[INFO] Loaded popmap: ", nrow(popmap_df), " rows")
    message("[INFO] Unique populations in popmap: ", paste(unique(popmap_df$pop), collapse = ", "))
  } else {
    warning("Could not parse popmap at: ", POPMAP)
  }
} else {
  warning("Popmap not found at: ", POPMAP, "  (script will still try to match by sample name)")
}

# Colors
col_HYD   <- "#DE3163"  
col_BKN   <- "#008B8B"  
col_EG35  <- "#B5A642"  
col_UNK   <- "#7F7F7F"  

match_pop_for_sample_v2 <- function(samp, pm_df) {
  if (is.null(pm_df) || nrow(pm_df) == 0) return(NA_character_)
  samp_clean <- trimws(samp)
  # 1) exact match
  idx <- which(pm_df$sample == samp_clean)
  if (length(idx) > 0) return(pm_df$pop[idx[1]])
  # 2) base duplicate patterns
  parts <- strsplit(samp_clean, "_")[[1]]
  if (length(parts) >= 2) {
    base2 <- paste(parts[1:2], collapse = "_")
    idx <- which(pm_df$sample == base2)
    if (length(idx) > 0) return(pm_df$pop[idx[1]])
    if (length(parts) >= 3) {
      base3 <- paste(parts[1:3], collapse = "_")
      idx <- which(pm_df$sample == base3)
      if (length(idx) > 0) return(pm_df$pop[idx[1]])
    }
  }
  # 3) case-insensitive
  idx <- which(tolower(pm_df$sample) == tolower(samp_clean))
  if (length(idx) > 0) return(pm_df$pop[idx[1]])
  # 4) repeated pattern detection
  if (grepl("_(.+)_\\1", samp_clean)) {
    base_name <- sub("_(.+)_\\1", "\\1", samp_clean)
    idx <- which(pm_df$sample == base_name)
    if (length(idx) > 0) return(pm_df$pop[idx[1]])
  }
  # 5) last token fallback
  last_token <- tail(strsplit(samp_clean, "_")[[1]], 1)
  idx <- which(pm_df$sample == last_token)
  if (length(idx) > 0) return(pm_df$pop[idx[1]])
  return(NA_character_)
}

# Function to extract clean sample name
extract_clean_sample_name <- function(samp) {
  samp_clean <- trimws(samp)
  parts <- strsplit(samp_clean, "_")[[1]]

  if (length(parts) >= 4) {
    mid <- length(parts)/2
    first_half <- paste(parts[1:mid], collapse = "_")
    second_half <- paste(parts[(mid+1):length(parts)], collapse = "_")
    if (first_half == second_half) return(first_half)
  }
  if (length(parts) >= 2) {
    base2 <- paste(parts[1:2], collapse = "_")
    remaining <- paste(parts[3:length(parts)], collapse = "_")
    if (base2 == remaining) return(base2)
  }
  return(samp_clean)
}

# create color mapping for groups 
build_group_color_map <- function(groups) {
  groups_u <- unique(na.omit(groups))
  known_pops <- c("HYD", "BKN", "EG", "UNK", RENAMED_EG35_POP)
  others <- setdiff(sort(groups_u), known_pops)
  ordered <- c(intersect(known_pops, groups_u), others)
  color_map <- character(length(ordered)); names(color_map) <- ordered
  if ("HYD" %in% ordered) color_map["HYD"] <- col_HYD
  if ("BKN" %in% ordered) color_map["BKN"] <- col_BKN
  if ("EG" %in% ordered) color_map["EG"] <- col_EG35
  if ("UNK" %in% ordered) color_map["UNK"] <- col_UNK

  if (RENAMED_EG35_POP %in% ordered) color_map[RENAMED_EG35_POP] <- col_EG35
  remaining <- ordered[which(color_map == "")]
  if (length(remaining) > 0) {
    if (length(remaining) <= 12) base_pal <- brewer.pal(12, "Set3") else base_pal <- colorRampPalette(brewer.pal(12, "Set3"))(length(remaining))
    color_map[remaining] <- base_pal[1:length(remaining)]
  }
  return(color_map)
}

is_bkn25_outlier_in_pc13 <- function(df) {
  bkn25_row <- df[grepl("BKN_25", df$sample), ]
  if (nrow(bkn25_row) == 0) { message("[INFO] BKN_25 not found in dataset"); return(FALSE) }
  bkn_data <- df[df$group == "BKN", ]
  if (nrow(bkn_data) <= 2) { message("[INFO] Not enough BKN samples for outlier detection"); return(FALSE) }
  centroid_pc1 <- mean(bkn_data$PC1, na.rm = TRUE); centroid_pc3 <- mean(bkn_data$PC3, na.rm = TRUE)
  distances <- sqrt((bkn_data$PC1 - centroid_pc1)^2 + (bkn_data$PC3 - centroid_pc3)^2)
  bkn25_distance <- distances[grepl("BKN_25", bkn_data$sample)]
  threshold <- mean(distances, na.rm = TRUE) + 1.5 * sd(distances, na.rm = TRUE)
  is_outlier <- bkn25_distance > threshold
  if (is_outlier) message("[OUTLIER] BKN_25 is an outlier in PC1 vs PC3 (distance: ", round(bkn25_distance,3), ", threshold: ", round(threshold,3), ")")
  else message("[INFO] BKN_25 is not an outlier in PC1 vs PC3 (distance: ", round(bkn25_distance,3), ", threshold: ", round(threshold,3), ")")
  return(is_outlier)
}

# plotting helper
plot_pcs_specific_labeling <- function(df, eigvals, xcol, ycol, out_png, group_colors) {
  idx_x <- as.integer(sub("PC","", xcol)); idx_y <- as.integer(sub("PC","", ycol))
  prop <- eigvals / sum(eigvals)
  xlabel <- paste0(xcol, " (", sprintf("%.2f%%", 100 * prop[idx_x]), ")")
  ylabel <- paste0(ycol, " (", sprintf("%.2f%%", 100 * prop[idx_y]), ")")
  df$group <- factor(df$group, levels = names(group_colors))
  df$clean_name <- sapply(df$sample, extract_clean_sample_name)
  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]], color = group)) +
    geom_hline(yintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "gray90", linewidth = 0.3) +
    geom_point(size = 4, alpha = 0.85) +
    scale_color_manual(values = group_colors, name = "Population", drop = FALSE) +
    labs(x = xlabel, y = ylabel, title = paste0("PCA — ", tools::file_path_sans_ext(basename(out_png)))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "gray95"),
      panel.grid.minor = element_blank()
    )
  if (xcol == "PC1" && ycol == "PC3") {
    if (is_bkn25_outlier_in_pc13(df)) {
      bkn25_data <- df[grepl("BKN_25", df$sample), ]
      if (nrow(bkn25_data) > 0) {
        p <- p + geom_text_repel(data = bkn25_data, aes(label = clean_name), size = 3.5,
                                 box.padding = 0.8, point.padding = 0.5, segment.color = "gray40",
                                 segment.alpha = 0.8, min.segment.length = 0.2, max.overlaps = 10,
                                 show.legend = FALSE, fontface = "bold")
        message("[LABEL] Added BKN_25 label to PC1 vs PC3 plot")
      }
    }
  }
  ggsave(filename = out_png, plot = p, width = PNG_W, height = PNG_H, dpi = DPI)
  message("[INFO] Wrote: ", out_png)
}

# main loop over sets
for (nm in names(set_paths)) {
  vcf_path <- set_paths[[nm]]
  if (!file.exists(vcf_path)) {
    warning("VCF not found for set '", nm, "': ", vcf_path, "  -> skipping")
    next
  }
  message("----------------------------------------------------")
  message("[INFO] Processing set: ", nm)
  outdir_set <- file.path(OUTROOT, nm)
  dir.create(outdir_set, recursive = TRUE, showWarnings = FALSE)
  
  gds_file <- file.path(outdir_set, paste0(nm, ".gds"))
  message("[INFO] Converting VCF -> GDS: ", gds_file)
  
  snpgdsVCF2GDS(vcf.fn = vcf_path, out.fn = gds_file, method = "biallelic.only", verbose = TRUE)
  genofile <- snpgdsOpen(gds_file)
  on.exit({ try(snpgdsClose(genofile), silent = TRUE) }, add = TRUE)
  
  message("[INFO] Running PCA (SNPRelate)...")
  pca <- snpgdsPCA(genofile, autosome.only = FALSE, remove.monosnp = TRUE, verbose = TRUE)
  eigvals <- pca$eigenval
  write.table(eigvals, file = file.path(outdir_set, "eigenvalues.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  write_csv(tibble(eigen_index = seq_along(eigvals), eigenvalue = eigvals, prop_variance = eigvals / sum(eigvals)),
            file.path(outdir_set, "eigenvalues.csv"))
  
  sample_ids <- pca$sample.id
  message("[INFO] VCF samples: ", paste(sample_ids, collapse = ", "))
  
  pcs_df <- data.frame(sample = sample_ids,
                       PC1 = pca$eigenvect[,1],
                       PC2 = pca$eigenvect[,2],
                       PC3 = if (ncol(pca$eigenvect) >= 3) pca$eigenvect[,3] else NA_real_,
                       stringsAsFactors = FALSE)
  
  # map populations with improved matching
  message("[INFO] Matching populations...")
  pcs_df$pop <- vapply(pcs_df$sample, function(s) match_pop_for_sample_v2(s, popmap_df), FUN.VALUE = character(1))
  

  pcs_df$clean_sample <- sapply(pcs_df$sample, extract_clean_sample_name, USE.NAMES = FALSE)
  idx_eg35 <- which(pcs_df$clean_sample == ORIG_EG35_NAME | grepl(paste0("\\b", ORIG_EG35_NAME, "\\b"), pcs_df$sample))
  if (length(idx_eg35) > 0) {
    message("[INFO] Found sample(s) matching EG_35 pattern: ", paste(pcs_df$sample[idx_eg35], collapse = ", "))
    pcs_df$sample[idx_eg35] <- RENAMED_EG35
    pcs_df$pop[idx_eg35] <- RENAMED_EG35_POP
    message("[INFO] Renamed these sample(s) to: ", RENAMED_EG35, " and set pop -> ", RENAMED_EG35_POP)
  }
  
  # Handle NA populations and create group column
  pcs_df$group <- ifelse(is.na(pcs_df$pop) | pcs_df$pop == "", "UNK", pcs_df$pop)
  message("[INFO] Assigned groups: ")
  print(table(pcs_df$group, useNA = "always"))
  
  # build group->color map
  group_colors <- build_group_color_map(pcs_df$group)
  message("[INFO] Color mapping: ")
  print(group_colors)
  
  # save annotated PCA scores & colors: ensure sample column contains renamed sample
  pcs_df_out <- pcs_df %>% mutate(color = unname(group_colors[group]))
  write_csv(pcs_df_out, file.path(outdir_set, "pca_scores_with_group_and_color.csv"))
  write_csv(pcs_df %>% select(sample, pop, group, PC1, PC2, PC3), file.path(outdir_set, "pca_scores.csv"))
  
  # Save color mapping for reference
  color_ref <- data.frame(population = names(group_colors), color = group_colors, stringsAsFactors = FALSE)
  write_csv(color_ref, file.path(outdir_set, "population_colors.csv"))
  
  # plots (PC1 vs PC2 and PC1 vs PC3 if available)
  out_png12 <- file.path(outdir_set, paste0(nm, "_PC1_vs_PC2.png"))
  plot_pcs_specific_labeling(pcs_df, eigvals, "PC1", "PC2", out_png12, group_colors)
  
  if (!all(is.na(pcs_df$PC3))) {
    out_png13 <- file.path(outdir_set, paste0(nm, "_PC1_vs_PC3.png"))
    plot_pcs_specific_labeling(pcs_df, eigvals, "PC1", "PC3", out_png13, group_colors)
  } else {
    message("[WARN] PC3 not available for set: ", nm)
  }
  
  snpgdsClose(genofile)
  if (!KEEP_GDS) {
    try(file.remove(gds_file), silent = TRUE)
    message("[INFO] Removed GDS file: ", gds_file)
  }
  
  message("[DONE] outputs in: ", outdir_set)
}

message("[ALL DONE] PCA outputs in: ", OUTROOT)
