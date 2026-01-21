#!/usr/bin/env Rscript

# 20 ADMIXTURE results visuai=lization (Rscript)
# make_admixture_plots_with_cvSD.R
# Builds ancestry stacked plots from ADMIXTURE replicates and prints CV stdev (and n) on multi-K panels.


suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(viridis); library(stringr)
})

option_list <- list(
  make_option(c("-i","--indir"), type="character",
            
              default = "Path/to/admix_runs",
              help="Root directory containing admixture run folders"),
  make_option(c("-o","--outdir"), type="character",
              default = "Path/to/output_dir",
              help="Output directory for CSVs and plots"),
  make_option(c("-p","--popmap"), type="character",
            
              default = "Path/to/popmap.txt",
              help="Optional popmap (sample<TAB>pop)"),
  make_option(c("--topk"), type="character",
              
              default = "2,3,4,5,6",
              help="Comma-separated Ks to stack into a combined multi-K plot, e.g. '2,3' (default: '2,3,4,5,6')"),
  make_option(c("--dpi"), type="integer", default = 1200),
  make_option(c("--png_w"), type="numeric", default = 12),
  make_option(c("--png_h"), type="numeric", default = 6)
)
opt <- parse_args(OptionParser(option_list=option_list))

INDIR <- opt$indir; OUTDIR <- opt$outdir; POPMAP <- opt$popmap
TOPK <- if (!is.null(opt$topk) && nzchar(opt$topk)) as.integer(str_split(opt$topk, ",")[[1]]) else NULL
DPI <- opt$dpi; PNG_W <- opt$png_w; PNG_H <- opt$png_h

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

message("[INFO] Searching for replicates under: ", INDIR)
all_dirs <- list.dirs(INDIR, recursive = TRUE, full.names = TRUE)
rep_dirs <- all_dirs[grepl("/K[0-9]+/rep[0-9]+$", all_dirs)]
if (length(rep_dirs) == 0) stop("No replicate directories found under: ", INDIR, "  (expected pattern: <prefix>/K<k>/rep<r>)")

parse_repdir <- function(d) {
  prefix <- basename(dirname(dirname(d)))
  K <- as.integer(str_match(d, "/K([0-9]+)/rep[0-9]+$")[,2])
  rep <- as.integer(str_match(d, "/K[0-9]+/rep([0-9]+)$")[,2])
  tibble(repdir = d, prefix = prefix, K = K, rep = rep)
}
runs_df <- bind_rows(lapply(rep_dirs, parse_repdir))

extract_cv_seed <- function(repdir) {
  candidates <- list.files(repdir, pattern = "admixture.stderr|admixture.stdout|\\.stderr$|\\.stdout$|admix.*\\.log$|\\.admix\\.log$|\\.log$",
                           full.names = TRUE, ignore.case = TRUE)
  candidates <- unique(c(candidates, list.files(repdir, full.names = TRUE)))
  cv <- NA_real_; seed <- NA_integer_
  for (f in candidates) {
    if (!file.exists(f)) next
    txt <- tryCatch(paste(readLines(f, warn = FALSE), collapse = "\n"), error = function(e) NA_character_)
    if (is.na(txt)) next
    # Try regex-first extraction of CV
    cv_matches <- str_match_all(txt, "(?i)(?:CV\\s*error[^0-9\\-]*|Cross-validation[^0-9\\-]*|CV[^0-9\\-]*[:=\\s]+)([0-9]+\\.[0-9]+)")
    if (length(cv_matches) && length(cv_matches[[1]])>0) {
      last <- tail(cv_matches[[1]][,2], 1)
      if (!is.na(last) && nzchar(last)) cv <- as.numeric(last)
    } else {
      lines_with_cv <- grep("(?i)CV|Cross-validation", str_split(txt, "\n")[[1]], value = TRUE)
      if (length(lines_with_cv)) {
        decs <- unlist(str_extract_all(paste(lines_with_cv, collapse = " "), "[0-9]+\\.[0-9]+"))
        if (length(decs)) cv <- as.numeric(tail(decs,1))
      }
    }
    seed_match <- str_match(txt, "(?i)Random seed[^0-9]*([0-9]+)")
    if (!is.na(seed_match[1,2])) seed <- as.integer(seed_match[1,2])
  }
  list(cv = cv, seed = seed)
}

read_Q_file <- function(repdir) {
  qfiles <- list.files(repdir, pattern = "(?:\\.Q$|\\.[0-9]+\\.Q$)", full.names = TRUE)
  if (length(qfiles) == 0) qfiles <- list.files(repdir, pattern = "\\.Q$", full.names = TRUE)
  if (length(qfiles) == 0) return(NULL)
  qsel <- qfiles[which.max(file.info(qfiles)$size)]
  mat <- tryCatch(as.matrix(read.table(qsel, header = FALSE, stringsAsFactors = FALSE)), error = function(e) NULL)
  if (is.null(mat)) return(NULL)
  list(mat = mat, qfile = qsel)
}

find_fam_for_rep <- function(repdir) {
  fams <- list.files(repdir, pattern="\\.fam$", full.names = TRUE)
  if (length(fams) > 0) return(fams[1])
  fams2 <- list.files(dirname(repdir), pattern="\\.fam$", full.names = TRUE)
  if (length(fams2) > 0) return(fams2[1])
  fams3 <- list.files(dirname(dirname(repdir)), pattern="\\.fam$", full.names = TRUE)
  if (length(fams3) > 0) return(fams3[1])
  return(NA_character_)
}

count_variants_for_rep <- function(repdir) {
  bim <- list.files(repdir, pattern = "work\\.bim$|.*\\.bim$", full.names = TRUE)
  if (length(bim) == 0) return(NA_integer_)
  w <- grep("/work\\.bim$", bim, value = TRUE)
  if (length(w) > 0) bim <- w[1] else bim <- bim[1]
  n <- tryCatch(as.integer(length(readLines(bim, warn = FALSE))), error = function(e) NA_integer_)
  return(n)
}

map_sample_label <- function(labels) {
  if (length(labels) == 0) return(labels)
  labels[labels == "EG_35"] <- "BKN_UV_o_01"
  return(labels)
}

# labels_from_fam: avoid duplicated FID_IID and apply mapping
labels_from_fam <- function(famfile, n_expected = NA_integer_) {
  if (is.na(famfile) || !file.exists(famfile)) {
    if (!is.na(n_expected)) return(make.unique(paste0("ind", seq_len(n_expected))))
    return(character(0))
  }
  fam <- read.table(famfile, stringsAsFactors = FALSE, header = FALSE)
  fid <- as.character(fam[[1]]); iid <- as.character(fam[[2]])
  labels <- ifelse(is.na(fid) | fid == "" | fid == "0" | fid == iid, iid, paste0(fid, "_", iid))
  labels <- make.unique(labels, sep = "_")
  labels <- map_sample_label(labels)
  return(labels)
}

read_Q_matrix_safe <- function(repdir) {
  qinfo <- read_Q_file(repdir)
  if (is.null(qinfo)) return(NULL)
  return(qinfo$mat)
}

meta_list <- lapply(seq_len(nrow(runs_df)), function(i) {
  r <- runs_df[i, ]
  cs <- extract_cv_seed(r$repdir)
  qinfo <- read_Q_file(r$repdir)
  famfile <- find_fam_for_rep(r$repdir)
  nclusters <- if (!is.null(qinfo)) ncol(qinfo$mat) else NA_integer_
  nsamples <- if (!is.null(qinfo)) nrow(qinfo$mat) else NA_integer_
  nvars <- count_variants_for_rep(r$repdir)
  tibble(repdir = r$repdir, prefix = r$prefix, K = r$K, rep = r$rep,
         qfile = if (!is.null(qinfo)) qinfo$qfile else NA_character_,
         q_present = !is.null(qinfo), q_rows = nsamples, q_cols = nclusters,
         fam = if (!is.na(famfile)) famfile else NA_character_, cv = cs$cv, seed = cs$seed,
         n_variants = nvars)
})
meta <- bind_rows(meta_list)
meta_out_csv <- file.path(OUTDIR, "admixture_replicates_all.csv")
write_csv(meta, meta_out_csv)
message("[INFO] Wrote metadata for all replicates: ", meta_out_csv)

# choose best replicate per prefix+K (lowest CV first, NA last)
chosen <- meta %>%
  group_by(prefix, K) %>%
  arrange(is.na(cv), cv, rep) %>%
  slice(1) %>%
  ungroup()
write_csv(chosen, file.path(OUTDIR, "admixture_replicates_chosen.csv"))
message("[INFO] Wrote chosen replicates: ", file.path(OUTDIR, "admixture_replicates_chosen.csv"))

# compute CV statistics (mean, sd, var, n) per prefix+K based on all replicates
cv_stats <- meta %>%
  group_by(prefix, K) %>%
  summarise(cv_n = sum(!is.na(cv)),
            cv_mean = ifelse(cv_n>0, mean(cv, na.rm = TRUE), NA_real_),
            cv_sd = ifelse(cv_n >= 2, sd(cv, na.rm = TRUE), NA_real_),
            cv_var = ifelse(cv_n >= 2, var(cv, na.rm = TRUE), NA_real_),
            .groups = "drop")

# attach cv stats to chosen table for easy access when plotting multi-K
chosen2 <- left_join(chosen, cv_stats, by = c("prefix", "K"))

# optional popmap
popmap_df <- NULL
if (file.exists(POPMAP)) {
  pm <- tryCatch(read.table(POPMAP, header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(pm) && ncol(pm) >= 2) {
    colnames(pm)[1:2] <- c("sample","pop")
    popmap_df <- pm[,c("sample","pop")]

    popmap_df$sample <- map_sample_label(popmap_df$sample)
    message("[INFO] Loaded popmap: ", POPMAP)
  } else message("[WARN] popmap exists but could not parse: ", POPMAP)
} else message("[INFO] popmap not found: ", POPMAP)

make_stacked_plot <- function(qmat, sample_labels, prefix_label, Kval, repnum, cv = NA_real_, seed = NA_integer_,
                              outdir = OUTDIR, dpi = DPI, png_w = PNG_W, png_h = PNG_H, popmap = popmap_df) {
  nK <- ncol(qmat)
  if (is.null(sample_labels) || length(sample_labels) != nrow(qmat)) {
    sample_labels <- make.unique(paste0("ind", seq_len(nrow(qmat))))
  }
  sample_labels <- make.unique(sample_labels, sep = "_")
  # apply hard mapping
  sample_labels <- map_sample_label(sample_labels)
  
  dfQ <- as.data.frame(qmat, stringsAsFactors = FALSE)
  colnames(dfQ) <- paste0("Cluster", seq_len(nK))
  dfQ$sample <- sample_labels
  
  if (!is.null(popmap)) {
    pm_match <- popmap$pop[match(dfQ$sample, popmap$sample)]
    dfQ$pop <- ifelse(is.na(pm_match), "UNK", pm_match)
    dfQ <- dfQ %>% arrange(pop, sample)
  } else dfQ <- dfQ %>% arrange(sample)
  
  df_long <- dfQ %>% pivot_longer(cols = starts_with("Cluster"), names_to = "cluster", values_to = "q")
  df_long$samplef <- factor(df_long$sample, levels = unique(dfQ$sample))
  pal <- viridis::viridis(n = nK, option = "plasma")
  subtitle_text <- paste0("rep=", repnum, if (!is.na(seed)) paste0(" seed=", seed) else "", if (!is.na(cv)) paste0("  CV=", signif(cv, 4)) else "")
  
  p <- ggplot(df_long, aes(x = samplef, y = q, fill = cluster)) +
    geom_bar(stat = "identity", width = 1, color = NA) +
    scale_fill_manual(values = pal, name = "Cluster") +
    labs(x = "", y = "Ancestry proportion", title = paste0(prefix_label, "  K=", Kval), subtitle = subtitle_text) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          panel.grid = element_blank(), legend.position = "top")
  
  if (!is.null(popmap)) {
    ord_samples <- unique(dfQ$sample); ord_pops <- dfQ$pop[match(ord_samples, dfQ$sample)]
    if (length(ord_pops) > 1) {
      boundaries <- which(ord_pops[-1] != ord_pops[-length(ord_pops)]) + 0.5
      if (length(boundaries) > 0) p <- p + geom_vline(xintercept = boundaries, colour = "grey40", size = 0.25)
    }
  }
  
  out_prefix <- file.path(outdir, paste0(gsub("[/ ]", "_", prefix_label), ".K", Kval, ".rep", repnum))
  pngfile <- paste0(out_prefix, ".png"); pdffile <- paste0(out_prefix, ".pdf"); svgfile <- paste0(out_prefix, ".svg")
  ggsave(pngfile, p, width = png_w, height = png_h, dpi = dpi)
  ggsave(pdffile, p, width = png_w, height = png_h)
  ggsave(svgfile, p, width = png_w, height = png_h)
  return(list(png = pngfile, pdf = pdffile, svg = svgfile))
}

# Create individual plots for chosen replicates
plots_made <- list()
for (i in seq_len(nrow(chosen2))) {
  row <- chosen2[i, ]
  if (!row$q_present) {
    message("[WARN] No .Q present for chosen replicate: ", row$repdir, " (skipping)"); next
  }
  qmat <- read_Q_matrix_safe(row$repdir)
  if (is.null(qmat)) { message("[WARN] Could not read Q for: ", row$repdir); next }
  famfile <- if (!is.na(row$fam) && file.exists(row$fam)) row$fam else find_fam_for_rep(row$repdir)
  sample_labels <- labels_from_fam(famfile, n_expected = nrow(qmat))
  sample_labels <- map_sample_label(sample_labels)
  out <- make_stacked_plot(qmat, sample_labels, row$prefix, row$K, row$rep, cv = row$cv, seed = row$seed,
                           outdir = OUTDIR, dpi = DPI, png_w = PNG_W, png_h = PNG_H, popmap = popmap_df)
  message("[INFO] Saved plots for: ", paste0(row$prefix, " K=", row$K, " rep=", row$rep), " -> ", out$png)
  plots_made[[paste0(row$prefix, "_K", row$K)]] <- out
}

# Combined multi-K (use chosen2 which includes cv stats)
if (!is.null(TOPK) && length(TOPK)>0) {
  message("[INFO] Creating combined multi-K plots for Ks: ", paste(TOPK, collapse = ","))
  prefixes <- unique(chosen2$prefix)
  for (pref in prefixes) {
    sel <- chosen2 %>% filter(prefix==pref & K %in% TOPK) %>% arrange(K)
    if (nrow(sel)==0) next
    plot_list <- list()
    for (r in seq_len(nrow(sel))) {
      row <- sel[r,]; qmat <- read_Q_matrix_safe(row$repdir); if (is.null(qmat)) next
      famfile <- if (!is.na(row$fam) && file.exists(row$fam)) row$fam else find_fam_for_rep(row$repdir)
      sample_labels <- labels_from_fam(famfile, n_expected = nrow(qmat))
      sample_labels <- map_sample_label(sample_labels)
      dfQ <- as.data.frame(qmat); colnames(dfQ) <- paste0("Cluster", seq_len(ncol(qmat))); dfQ$sample <- sample_labels
      if (!is.null(popmap_df)) {
        pm_match <- popmap_df$pop[match(dfQ$sample, popmap_df$sample)]
        dfQ$pop <- ifelse(is.na(pm_match), "UNK", pm_match); dfQ <- dfQ %>% arrange(pop, sample)
      } else dfQ <- dfQ %>% arrange(sample)
      df_long <- dfQ %>% pivot_longer(cols = starts_with("Cluster"), names_to = "cluster", values_to = "q")
      df_long$samplef <- factor(df_long$sample, levels = unique(dfQ$sample))
      pal <- viridis::viridis(n = ncol(qmat), option = "plasma")
      subtitle_text <- paste0("rep=", row$rep, if (!is.na(row$seed)) paste0(" seed=", row$seed) else "", if (!is.na(row$cv)) paste0("  CV=", signif(row$cv,4)) else "")
      # CV stats for this K (from chosen2 join)
      cvsd_text <- if (!is.na(row$cv_sd)) paste0(" CVsd=", signif(row$cv_sd, 4), " (n=", row$cv_n, ")") else paste0(" CVsd=NA (n=", row$cv_n, ")")
      subtitle_full <- paste0(subtitle_text, cvsd_text)
      p <- ggplot(df_long, aes(x = samplef, y = q, fill = cluster)) +
        geom_bar(stat="identity", width=1) +
        scale_fill_manual(values=pal, name="Cluster") +
        labs(x = "", y = " ancestry", title = paste0(pref, "  K=", row$K), subtitle = subtitle_full) +
        theme_minimal(base_size=12) +
        theme(axis.text.x=element_text(angle=90,size=6), panel.grid=element_blank(), legend.position="top")
      # also annotate on top-right of panel
      ann_txt <- paste0("CVsd=", ifelse(is.na(row$cv_sd), "NA", signif(row$cv_sd,4)), "  n=", row$cv_n)
      p <- p + annotate("text", x = Inf, y = 1.02, label = ann_txt, hjust = 1.05, vjust = 0, size = 3)
      plot_list[[paste0("K",row$K)]] <- p
    }
    if (length(plot_list)==0) next
    # combine vertically - prefer patchwork, fallback to cowplot
    if (requireNamespace("patchwork", quietly=TRUE)) {
      library(patchwork)
      comb <- wrap_plots(plot_list, ncol=1)
      comb <- comb + plot_annotation(title = paste0(pref, " combined K=", paste(TOPK, collapse=",")))
      out_comb <- file.path(OUTDIR, paste0(pref, ".combined_K", paste(TOPK, collapse="_"), ".png"))
      ggsave(out_comb, comb, width = PNG_W, height = PNG_H*length(plot_list), dpi = DPI)
      message("[INFO] Wrote combined plot via patchwork: ", out_comb)
    } else if (requireNamespace("cowplot", quietly=TRUE)) {
      library(cowplot)
      grobs <- lapply(plot_list, ggplotGrob)
      comb_grob <- plot_grid(plotlist = grobs, ncol = 1, align = "v")
      title <- ggdraw() + draw_label(paste0(pref, " combined K=", paste(TOPK, collapse=",")), fontface='bold')
      final <- plot_grid(title, comb_grob, ncol = 1, rel_heights = c(0.05, 0.95))
      out_comb <- file.path(OUTDIR, paste0(pref, ".combined_K", paste(TOPK, collapse="_"), ".png"))
      ggsave(out_comb, final, width = PNG_W, height = PNG_H*length(plot_list), dpi = DPI)
      message("[INFO] Wrote combined plot via cowplot: ", out_comb)
    } else {
      message("[WARN] Neither patchwork nor cowplot installed; cannot make combined multi-K image. Install patchwork or cowplot.")
    }
  }
}

message("[DONE] Output in: ", OUTDIR)
