#!/usr/bin/env Rscript


# Fst_plots (Rscript)
# make_all_fst_plots.R


suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(ape)
  library(ggrepel)
  library(vcfR)
})

option_list <- list(
  make_option(c("-i","--indir"), type="character",
              default = "path/to/input/directory",
              help="Directory with per-set fst result folders"),
  make_option(c("-p","--popmap"), type="character",
              default = "popmap_from_fam.txt",
              help="popmap sample<TAB>pop"),
  make_option(c("-o","--outdir"), type="character", default = NULL,
              help="(optional) root output dir; default = --indir"),
  make_option(c("--dpi"), type="integer", default=1200, help="PNG DPI (default 1200)"),
  make_option(c("--min_sites_indfst"), type="integer", default=100,
              help="Min overlapping sites for individual FST (default 100)"),
  make_option(c("--verbose"), type="logical", default=TRUE, help="Verbose output")
)
opt <- parse_args(OptionParser(option_list=option_list))

INDIR <- normalizePath(opt$indir, mustWork = TRUE)
POPFILE_CMDLINE <- opt$popmap
POPFILE_CLEAN <- file.path(INDIR, "popmap.cleaned.txt")
POPFILE <- if (file.exists(POPFILE_CLEAN)) POPFILE_CLEAN else POPFILE_CMDLINE
if (!file.exists(POPFILE)) {
  warning("Provided popmap not found at either: ", POPFILE_CLEAN, " or ", POPFILE_CMDLINE, ". Population annotations will be inferred from sample names.")
  POPFILE <- NULL
}
OUTROOT <- ifelse(is.null(opt$outdir), INDIR, opt$outdir)
DPI <- opt$dpi
MIN_SITES <- opt$min_sites_indfst
VERBOSE <- opt$verbose

logmsg <- function(...) if(VERBOSE) message("[", format(Sys.time(), "%F %T"), "] ", ...)

# helper functions
ensure_dir <- function(p) { if(!dir.exists(p)) dir.create(p, recursive = TRUE); normalizePath(p) }
safe_pcoa <- function(distmat) {
  res <- NULL
  try({ res <- ape::pcoa(as.matrix(distmat), correction = "cailliez") }, silent = TRUE)
  if(is.null(res)) {
    d <- as.dist(distmat)
    km <- cmdscale(d, k = min(3, nrow(distmat)-1), eig = TRUE)
    res <- list(vectors = km$points, values = list(Relative_eig = (km$eig / sum(km$eig, na.rm = TRUE))))
  }
  return(res)
}
read_fst_file <- function(path) {
  df <- tryCatch(read.table(path, header = TRUE, comment.char = "#", stringsAsFactors = FALSE),
                 error = function(e) tryCatch(read.table(path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE),
                                              error = function(e2) stop("Cannot read FST file: ", path)))
  col_idx <- grep("WEIR|WEIR_AND_COCKERHAM|FST", colnames(df), ignore.case = TRUE)
  if(length(col_idx) == 0) {
    if(ncol(df) >= 1 && all(grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$|^NaN$|^nan$", as.character(df[[ncol(df)]]), perl=TRUE), na.rm=TRUE)) {
      col_idx <- ncol(df)
    } else if (ncol(df) >= 4) {
      col_idx <- 4
    } else stop("No FST-like column detected in ", path)
  }
  vec <- suppressWarnings(as.numeric(as.character(df[[col_idx[1]]])))
  vec <- vec[!is.na(vec)]
  return(vec)
}

# clean duplicated sample names
clean_sample_name <- function(s) {
  if (is.na(s)) return(NA_character_)
  s2 <- gsub("^(.+?)_\\1$", "\\1", s, perl = TRUE)
  s3 <- gsub("^(.+?)_\\1$", "\\1", s2, perl = TRUE)
  s3 <- trimws(s3)
  return(s3)
}

ORIG_EG35_NAME <- "EG_35"
RENAMED_EG35 <- "BKN_UV_o_01"
RENAMED_EG35_POP <- "BKN_o"
EXCLUDE_POP <- "EG"           
EXCLUDE_SAMPLE_ORIG <- ORIG_EG35_NAME  
# helper to apply the strict sample mapping
map_sample_label_hard <- function(x) {
  x <- as.character(x)
  x_clean <- sapply(x, clean_sample_name, USE.NAMES = FALSE)
  x_clean[x_clean == ORIG_EG35_NAME] <- RENAMED_EG35
  return(x_clean)
}

popmap_df <- NULL
if (!is.null(POPFILE) && file.exists(POPFILE)) {
  pm_raw <- tryCatch(read.table(POPFILE, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE),
                     error = function(e) read.table(POPFILE, header = FALSE, stringsAsFactors = FALSE, fill = TRUE))
  if (ncol(pm_raw) >= 2) {
    colnames(pm_raw)[1:2] <- c("sample_raw","pop_raw")
    pm_raw$sample <- sapply(pm_raw$sample_raw, clean_sample_name, USE.NAMES = FALSE)
    # Apply renaming: EG_35 -> BKN_UV_o_01
    pm_raw$sample[pm_raw$sample == ORIG_EG35_NAME] <- RENAMED_EG35
    pm_raw$pop_raw[pm_raw$sample == RENAMED_EG35] <- RENAMED_EG35_POP
    popmap_df <- pm_raw %>% select(sample, pop = pop_raw)
    logmsg("Loaded popmap and applied EG_35 -> BKN_UV_o_01 mapping.")
  } else {
    logmsg("Popmap present but could not parse into two columns; ignoring.")
    popmap_df <- NULL
  }
} else {
  logmsg("No popmap found; will infer populations from sample names.")
}

OUTROOT <- ensure_dir(OUTROOT)
logmsg("INDIR =", INDIR)
if (!is.null(popmap_df)) logmsg("Using popmap with mapping applied.")
logmsg("Output root:", OUTROOT)

set_dirs <- list.dirs(INDIR, full.names = TRUE, recursive = FALSE)
if(length(set_dirs) == 0) stop("No subdirectories found in ", INDIR)

for(setdir in set_dirs) {
  setname <- basename(setdir)
  logmsg("Processing set: ", setname)
  outset <- file.path(OUTROOT, setname)
  plots_dir <- file.path(outset, "plots")
  ensure_dir(outset); ensure_dir(plots_dir)
  

  # Population-level mean FST (exclude population EG)

  fst_files <- list.files(setdir, pattern="\\.weir\\.fst$", full.names = TRUE)
  if(length(fst_files) > 0) {
    logmsg(" Found ", length(fst_files), " population FST files.")
    res_list <- vector("list", length = length(fst_files))
    k <- 1
    for(f in fst_files) {
      nm <- basename(f) %>% sub("^fst_","",.) %>% sub("\\.weir\\.fst$","",.) %>% sub("\\.vcftools.*$","",.)
      pops <- strsplit(nm, "_")[[1]]
      if(length(pops) < 2) pops <- strsplit(nm, "\\.")[[1]]
      pa <- pops[1]; pb <- ifelse(length(pops)>=2, pops[2], NA_character_)
      vals <- tryCatch(read_fst_file(f), error = function(e) { logmsg("  Error reading ", f, ": ", e$message); numeric(0) })
      res_list[[k]] <- tibble(file = basename(f), popA = pa, popB = pb,
                              mean_fst = ifelse(length(vals)>0, mean(vals, na.rm=TRUE), NA_real_),
                              median_fst = ifelse(length(vals)>0, median(vals, na.rm=TRUE), NA_real_),
                              n_sites = length(vals),
                              all_vals = list(vals))
      k <- k + 1
    }
    res_df_all <- bind_rows(res_list)
    # Filter out any population pairs involving EXCLUDE_POP
    res_df <- res_df_all %>% filter(!(popA == EXCLUDE_POP | popB == EXCLUDE_POP))
    
    if(nrow(res_df) > 0) {
      pops <- sort(unique(c(res_df$popA, res_df$popB)))
      pops <- pops[!is.na(pops)]
      mat <- matrix(NA_real_, nrow=length(pops), ncol=length(pops), dimnames=list(pops,pops))
      for(i in seq_len(nrow(res_df))) {
        a <- res_df$popA[i]; b <- res_df$popB[i]; v <- res_df$mean_fst[i]
        if(!is.na(a) && !is.na(b)) {
          mat[a,b] <- v; mat[b,a] <- v
        }
      }
      diag(mat) <- 0
      write.csv(as.data.frame(mat), file=file.path(outset,"pairwise_FST_matrix_mean.csv"), row.names=TRUE)
      write_csv(res_df %>% select(-all_vals), file.path(outset,"fst_pairwise_summary_per_file.filtered_popEG.csv"))
      
      # heatmap (population mean) - using filtered pops
      mat_plot <- mat; mat_plot[is.na(mat_plot)] <- 0
      heat_png <- file.path(plots_dir,"FST_heatmap_mean.filtered_popEG.png")
      heat_svg <- file.path(plots_dir,"FST_heatmap_mean.filtered_popEG.svg")
      pal <- colorRampPalette(rev(brewer.pal(9,"YlGnBu")))(200)
      png(heat_png, width=9, height=7, units="in", res=DPI)
      pheatmap(mat_plot, color=pal, cluster_rows=TRUE, cluster_cols=TRUE,
               display_numbers=TRUE, number_format="%.3f",
               main=paste0(setname," — mean pairwise Fst"), fontsize_number=8, fontsize=12)
      dev.off()
      svg(heat_svg, width=9, height=7)
      pheatmap(mat_plot, color=pal, cluster_rows=TRUE, cluster_cols=TRUE,
               display_numbers=TRUE, number_format="%.3f",
               main=paste0(setname," — mean pairwise Fst"), fontsize_number=8, fontsize=12)
      dev.off()
      
      # PCoA from filtered mean pairwise FST
      mat_for_pcoa <- mat
      if(any(is.na(mat_for_pcoa))) {
        logmsg("  Filling NA in mean matrix with column means (for PCoA only).")
        col_means <- apply(mat_for_pcoa, 2, function(x) mean(x, na.rm=TRUE))
        for(r in seq_len(nrow(mat_for_pcoa))) for(c in seq_len(ncol(mat_for_pcoa))) if(is.na(mat_for_pcoa[r,c])) mat_for_pcoa[r,c] <- col_means[c]
        diag(mat_for_pcoa) <- 0
      }
      pcoa_res <- tryCatch(safe_pcoa(mat_for_pcoa), error = function(e) NULL)
      if(!is.null(pcoa_res)) {
        vecs <- pcoa_res$vectors
        rel_eig <- if(!is.null(pcoa_res$values$Relative_eig)) pcoa_res$values$Relative_eig else (pcoa_res$values$Eigenvalues / sum(pcoa_res$values$Eigenvalues, na.rm=TRUE))
        var1 <- round(100 * rel_eig[1], 2); var2 <- ifelse(length(rel_eig)>=2, round(100 * rel_eig[2], 2), 0)
        mds_df <- data.frame(POP = rownames(mat_for_pcoa), PCo1 = vecs[,1], PCo2 = ifelse(ncol(vecs)>=2, vecs[,2], rep(0,nrow(vecs))))
        pcoa_png <- file.path(plots_dir,"FST_PCoA.filtered_popEG.png"); pcoa_svg <- file.path(plots_dir,"FST_PCoA.filtered_popEG.svg")
        png(pcoa_png, width=10, height=8, units="in", res=DPI)
        gg <- ggplot(mds_df, aes(x=PCo1, y=PCo2, label=POP)) +
          geom_point(size=4, color="steelblue") + ggrepel::geom_text_repel(size=4, max.overlaps=50) +
          labs(x=paste0("PCo1 (",var1,"% )"), y=paste0("PCo2 (",var2,"% )"), title=paste0("PCoA from mean pairwise Fst — ", setname, "")) +
          theme_minimal(base_size=14)
        print(gg); dev.off()
        svg(pcoa_svg, width=10, height=8); print(gg); dev.off()
      } else {
        logmsg("  PCoA failed for set ", setname)
      }
      
      # Expand per-site FST values
      expanded <- bind_rows(lapply(seq_len(nrow(res_df)), function(i) {
        row <- res_df[i, ]
        vals <- ifelse(is.null(row$all_vals[[1]]), numeric(0), row$all_vals[[1]])
        tibble(pop_pair = paste0(row$popA,"_",row$popB), popA=row$popA, popB=row$popB, value = vals)
      }))
      if(nrow(expanded) > 0) {
        bp_png <- file.path(plots_dir, "FST_pairwise_values_violin_boxplot.filtered_popEG.png")
        png(bp_png, width=12, height=6, units="in", res=DPI)
        ggplot(expanded, aes(x=pop_pair, y=value)) +
          geom_violin(trim=TRUE) + geom_boxplot(width=0.12, outlier.size=0.8) +
          labs(x="Population pair", y="FST per-site", title=paste0(setname, " per-site FST distribution by pair (EG excluded)")) +
          theme_minimal(base_size=12) + theme(axis.text.x = element_text(angle=45, hjust=1))
        dev.off()
      }
    } else {
      logmsg("  No usable population FST rows (after excluding EG) for set ", setname)
    }
  } else {
    logmsg("  No population-level .weir.fst files in ", setdir)
  }
  
  # Individual-level FST from pruned VCF
  pruned_vcfs <- list.files(setdir, pattern="\\.pruned\\.vcf\\.gz$|pruned_.*\\.vcf\\.gz$|.*\\.pruned\\.vcf\\.gz$", full.names = TRUE)
  if(length(pruned_vcfs) == 0) {
    logmsg("  No pruned VCF found in ", setdir, " — skip individual-level FST.")
    next
  }
  pruned_vcf <- pruned_vcfs[1]
  logmsg("  Reading pruned VCF: ", pruned_vcf)
  vcf <- tryCatch(vcfR::read.vcfR(pruned_vcf, verbose=FALSE), error = function(e) { logmsg("   VCF read error: ", e$message); return(NULL) })
  if(is.null(vcf)) next
  
  gt_mat <- tryCatch(vcfR::extract.gt(vcf, element="GT"), error = function(e) { logmsg("   extract.gt error: ", e$message); return(NULL) })
  if(is.null(gt_mat)) next
  gt_mat <- as.matrix(gt_mat)
  if(is.null(colnames(gt_mat)) && !is.null(rownames(gt_mat))) gt_mat <- t(gt_mat)
  samples_raw <- colnames(gt_mat)
  if(is.null(samples_raw) || length(samples_raw) < 2) {
    logmsg("   Unable to determine sample names from VCF GT matrix; skipping individual-level FST for ", setname)
    next
  }
  # Clean sample names
  orig_cleaned_samples <- sapply(samples_raw, clean_sample_name, USE.NAMES = FALSE)
  new_sample_names <- orig_cleaned_samples
  new_sample_names[new_sample_names == ORIG_EG35_NAME] <- RENAMED_EG35
  colnames(gt_mat) <- new_sample_names
  
  
  sample_map_df <- tibble(original = orig_cleaned_samples, sample = new_sample_names)
  
  # If popmap exists, ensure its sample names are mapped identically and EG_35 row changed
  if(!is.null(popmap_df)) {
    if(!(RENAMED_EG35 %in% popmap_df$sample)) {
      popmap_df <- bind_rows(popmap_df, tibble(sample = RENAMED_EG35, pop = RENAMED_EG35_POP))
      logmsg("  Added mapping for renamed sample to popmap: ", RENAMED_EG35, " -> ", RENAMED_EG35_POP)
    } else {
      popmap_df$pop[popmap_df$sample == RENAMED_EG35] <- RENAMED_EG35_POP
    }
  }
  
  # Convert GT strings to numeric (0,1,2)
  gt_to_num_vec <- function(g) {
    g[is.na(g)] <- "./."
    g2 <- gsub("\\|", "/", g)
    r <- rep(NA_real_, length(g2))
    r[g2 %in% c("0/0","0|0")] <- 0
    r[g2 %in% c("0/1","1/0","0|1","1|0")] <- 1
    r[g2 %in% c("1/1","1|1")] <- 2
    return(r)
  }
  Gnum <- apply(gt_mat, 2, gt_to_num_vec)
  if(!is.matrix(Gnum)) Gnum <- as.matrix(Gnum)
  nvar <- nrow(Gnum)
  logmsg("  Variants in pruned VCF: ", nvar)
  
  compute_pairwise_fst_vec <- function(x, y, min_sites = MIN_SITES) {
    ok <- !is.na(x) & !is.na(y)
    n_ok <- sum(ok)
    if(n_ok < min_sites) return(NA_real_)
    x2 <- x[ok]; y2 <- y[ok]
    p1 <- mean(x2)/2; p2 <- mean(y2)/2; p_avg <- (p1 + p2)/2
    if(is.na(p1) || is.na(p2) || is.na(p_avg)) return(NA_real_)
    if(p_avg == 0 || p_avg == 1) return(0)
    MSG <- (p1*(1-p1) + p2*(1-p2)) / 2
    MSP <- ((p1 - p_avg)^2 + (p2 - p_avg)^2)
    denom <- (MSP + 2*MSG)
    if(denom == 0) return(0)
    FST <- MSP / denom
    if(is.na(FST) || FST < 0) return(0)
    return(FST)
  }
  
  samples <- colnames(Gnum)
  n <- length(samples)
  fst_mat_ind <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(samples, samples))
  for(i in seq_len(n)) {
    fst_mat_ind[i,i] <- 0
    if(i < n) {
      for(j in seq(i+1, n)) {
        v <- compute_pairwise_fst_vec(Gnum[,i], Gnum[,j], min_sites = MIN_SITES)
        fst_mat_ind[i,j] <- v; fst_mat_ind[j,i] <- v
      }
    }
  }
  
  write.csv(as.data.frame(fst_mat_ind), file = file.path(outset, "individual_pairwise_FST_matrix.full.csv"), row.names = TRUE)
  
  # per-sample summary (full)
  fst_summary <- tibble(Sample = samples,
                        Mean_FST = rowMeans(fst_mat_ind, na.rm = TRUE),
                        Min_FST = apply(fst_mat_ind, 1, function(x) if(all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)),
                        Max_FST = apply(fst_mat_ind, 1, function(x) if(all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))

  fst_summary <- left_join(fst_summary, sample_map_df %>% rename(OriginalSample = original, Sample = sample), by = "Sample")

  if(!is.null(popmap_df)) {
    fst_summary <- left_join(fst_summary, popmap_df, by = c("Sample" = "sample"))
    colnames(fst_summary)[which(names(fst_summary) == "pop")] <- "pop"
  } else {
    fst_summary$pop <- sapply(fst_summary$OriginalSample, function(s) {
      if(is.na(s)) return("OTHER")
      if(grepl("^EV", s)) return("EV")
      if(grepl("^EG", s)) return("EG")
      if(grepl("^WR", s)) return("WR")
      return("OTHER")
    })
 
    fst_summary$pop[fst_summary$OriginalSample == ORIG_EG35_NAME] <- RENAMED_EG35_POP
  }
  write_csv(fst_summary, file.path(outset, "individual_FST_summary.full.csv"))
  

  # individual heatmap: INCLUDE ALL samples

  annot_all <- data.frame(Population = fst_summary$pop); rownames(annot_all) <- fst_summary$Sample
  unique_pops_all <- unique(annot_all$Population)
  cols_all <- RColorBrewer::brewer.pal(max(3, length(unique_pops_all)), "Set2")[1:length(unique_pops_all)]
  ann_colors_all <- list(Population = setNames(cols_all, unique_pops_all))
  hm_mat_plot_all <- fst_mat_ind
  hm_mat_plot_all[is.na(hm_mat_plot_all)] <- 0
  png(file.path(plots_dir,"individual_FST_heatmap.allSamples.png"), width=12, height=10, units="in", res=DPI)
  pheatmap(hm_mat_plot_all, color = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(200),
           cluster_rows=TRUE, cluster_cols=TRUE,
           show_rownames=TRUE, show_colnames=TRUE, annotation_row=annot_all, annotation_col=annot_all,
           annotation_colors=ann_colors_all, main=paste0(setname," — individual pairwise FST"),
           fontsize_row = ifelse(n>60,6,10))
  dev.off()
  svg(file.path(plots_dir,"individual_FST_heatmap.allSamples.svg"), width=12, height=10)
  pheatmap(hm_mat_plot_all, color = colorRampPalette(rev(brewer.pal(9,"RdYlBu")))(200),
           cluster_rows=TRUE, cluster_cols=TRUE,
           show_rownames=TRUE, show_colnames=TRUE, annotation_row=annot_all, annotation_col=annot_all,
           annotation_colors=ann_colors_all, main=paste0(setname," — individual pairwise FST"),
           fontsize_row = ifelse(n>60,6,10))
  dev.off()
  

  # MDS (include ALL samples)
  dmat_all <- fst_mat_ind
  if(any(is.na(dmat_all))) {
    col_means <- apply(dmat_all, 2, function(x) mean(x, na.rm = TRUE))
    for(i in seq_len(nrow(dmat_all))) for(j in seq_len(ncol(dmat_all))) if(is.na(dmat_all[i,j])) dmat_all[i,j] <- (col_means[i] + col_means[j]) / 2
  }
  diag(dmat_all) <- 0
  dmat_all <- (dmat_all + t(dmat_all)) / 2
  mds <- tryCatch(cmdscale(as.dist(dmat_all), k = 2, eig = TRUE), error = function(e) NULL)
  if(!is.null(mds)) {
    mdf <- data.frame(Sample = rownames(dmat_all), MDS1 = mds$points[,1], MDS2 = mds$points[,2], pop = fst_summary$pop)
    png(file.path(plots_dir,"individual_MDS.allSamples.png"), width=8, height=6, units="in", res=DPI)
    gg <- ggplot(mdf, aes(x=MDS1, y=MDS2, color = pop, label = Sample)) +
      geom_point(size=3) + ggrepel::geom_text_repel(size=3, max.overlaps = 100) +
      labs(title=paste0(setname," — Individual MDS (from pairwise FST)")) + theme_minimal(base_size = 12)
    print(gg); dev.off()
    svg(file.path(plots_dir,"individual_MDS.allSamples.svg"), width=8, height=6); print(gg); dev.off()
  }
  
 
  # For all other plots (population-level and individual summaries/histograms/boxplots),
  # EXCLUDE: any sample whose OriginalSample == ORIG_EG35_NAME and any sample with pop == EXCLUDE_POP.
-
  excluded_original_samples <- ORIG_EG35_NAME
  excluded_popname <- EXCLUDE_POP
  
  fst_summary_filtered <- fst_summary %>%
    filter(!(OriginalSample %in% excluded_original_samples) & !(pop == excluded_popname))
  write_csv(fst_summary_filtered, file.path(outset, "individual_FST_summary.filtered_excludingEGandEG35.csv"))
  
  # Per-sample mean-FST boxplot (filtered)
  if("pop" %in% colnames(fst_summary_filtered) && nrow(fst_summary_filtered) > 0) {
    png(file.path(plots_dir,"individual_meanFST_by_pop_boxplot.filtered.png"), width=8, height=5, units="in", res=DPI)
    gg <- ggplot(fst_summary_filtered, aes(x = pop, y = Mean_FST, fill = pop)) +
      geom_boxplot() + geom_jitter(width = 0.15, size = 2) +
      labs(title = paste0(setname, " — per-sample mean FST by population"), x = "Population", y = "Mean pairwise FST") +
      theme_minimal() + theme(legend.position = "none")
    print(gg); dev.off()
    svg(file.path(plots_dir,"individual_meanFST_by_pop_boxplot.filtered.svg"), width=8, height=5); print(gg); dev.off()
  }
  
  # Histogram of pairwise FST values (individual) - using filtered matrix
  
  samples_filtered <- fst_summary_filtered$Sample
  if(length(samples_filtered) >= 2) {
    dmat_filtered <- fst_mat_ind[samples_filtered, samples_filtered, drop = FALSE]
    dmat_filtered[is.na(dmat_filtered)] <- 0
    vals <- as.vector(dmat_filtered[upper.tri(dmat_filtered)])
    vals <- vals[!is.na(vals)]
    if(length(vals) > 0) {
      png(file.path(plots_dir,"individual_pairwiseFST_histogram.filtered.png"), width=8, height=4, units="in", res=DPI)
      gg <- ggplot(data.frame(FST = vals), aes(x = FST)) +
        geom_histogram(bins = 40, fill = "steelblue", color = "white") +
        labs(title = paste0(setname, " — histogram of individual pairwise FST values (EG excluded)"), x = "FST", y = "count") +
        theme_minimal()
      print(gg); dev.off()
      svg(file.path(plots_dir,"individual_pairwiseFST_histogram.filtered.svg"), width=8, height=4); print(gg); dev.off()
    }
  } else {
    logmsg("  Not enough samples after filtering for histogram/boxplot (EG excluded).")
  }
  
  
  logmsg("  Wrote outputs for set ", setname, " (filtered population-level and filtered individual-level plots).")
}

logmsg("ALL done. Outputs are under ", OUTROOT)
