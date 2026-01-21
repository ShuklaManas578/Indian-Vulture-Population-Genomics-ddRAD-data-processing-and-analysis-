#!/usr/bin/env Rscript

# 10 Plotting nucleotide diversity (Rscript)
# plot_pi_v2.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(viridis)
  library(scales)
})


args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 3) {
  merged_file <- args[1]
  outdir <- args[2]
  summary_csv <- args[3]
} else {
  wd <- getwd()
  message("[INFO] No command-line args: auto-detecting in ", wd)
  candidates <- list.files(wd, pattern = "(merged_window_pi.*\\.tsv$|.*_window_pi.*\\.tsv$)", full.names = TRUE, ignore.case = TRUE)
  if (length(candidates) == 0) stop("No merged window pi file found. Provide arguments or place file in working dir.")
  merged_file <- candidates[1]; message("[INFO] Using merged file: ", basename(merged_file))
  outdir <- file.path(wd, "pi_plots_BKN_HYD")
  summary_candidates <- list.files(wd, pattern = "(summary.*\\.csv$|.*_summary.*\\.csv$)", full.names = TRUE, ignore.case = TRUE)
  summary_csv <- if (length(summary_candidates) > 0) summary_candidates[1] else ""
  if (nzchar(summary_csv)) message("[INFO] Found summary CSV: ", basename(summary_csv)) else message("[INFO] No summary CSV found; using alphabetical ordering.")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Read input 
message("[INFO] Reading: ", merged_file)
df_all <- read_tsv(merged_file, col_types = cols(
  CHROM = col_character(),
  BIN_START = col_double(),
  BIN_END = col_double(),
  N_VARIANTS = col_double(),
  PI = col_double(),
  POP = col_character()
))

# basic cleaning
df_all <- df_all %>% filter(!is.na(PI) & PI >= 0)

df_all <- df_all %>% mutate(POP = ifelse(POP %in% c("EG35","EG_35","EG-35"), "EG_35", POP))

df_all <- df_all %>% mutate(orig_POP = POP)

# Remap EG_35 into BKN for analysis and tag remapped_sample
df_analysis <- df_all %>%
  mutate(
    remapped_sample = ifelse(orig_POP == "EG_35", "BKN_UV_o_01", NA_character_),
    POP = ifelse(orig_POP == "EG_35", "BKN", POP)
  )

df_plot <- df_analysis %>%
  filter(orig_POP != "EG_35") %>%   # exclude EG_35 rows from plotting
  filter(POP %in% c("BKN", "HYD"))

if (nrow(df_plot) == 0) stop("No data left for plotting after limiting to BKN & HYD (and excluding EG_35). Check POP labels in the merged file.")

# Save the remapped analysis data
write_csv(df_analysis, file.path(outdir, "pi_merged_windowed_for_analysis_including_EG35_remapped_BKN_UV_o_01.csv"))
message("[INFO] Wrote analysis-ready merged data (EG_35 remapped to BKN): ", file.path(outdir, "pi_merged_windowed_for_analysis_including_EG35_remapped_BKN_UV_o_01.csv"))

# Save the plotting-only data (BKN & HYD only, EG_35 excluded)
write_csv(df_plot, file.path(outdir, "pi_merged_windowed_for_plotting_BKN_HYD_only.csv"))
message("[INFO] Wrote plotting-only merged data (BKN & HYD only): ", file.path(outdir, "pi_merged_windowed_for_plotting_BKN_HYD_only.csv"))

# Population ordering: prefer summary_csv if given
if (nzchar(summary_csv) && file.exists(summary_csv)) {
  ssum <- read_csv(summary_csv, show_col_types = FALSE)
  ssum <- ssum %>% mutate(POP = ifelse(POP %in% c("EG35","EG_35","EG-35"), "EG_35", POP))
  pop_levels_plot <- intersect(c("BKN", "HYD"), intersect(ssum$POP, unique(df_plot$POP)))
  if (length(pop_levels_plot) == 0) pop_levels_plot <- unique(df_plot$POP)
  df_plot$POP <- factor(df_plot$POP, levels = pop_levels_plot)
} else {
  df_plot$POP <- factor(df_plot$POP, levels = sort(unique(df_plot$POP)))
}

# Colors
pop_colors <- viridis::cividis(n = length(levels(df_plot$POP)))
names(pop_colors) <- levels(df_plot$POP)

# DPI
DPI <- 1200

use_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

# Summary stats ----
# summary for analysis (includes EG_35 remapped into BKN)
summary_stats_analysis <- df_analysis %>%
  group_by(POP) %>%
  summarise(n = n(), mean_pi = mean(PI), sd_pi = sd(PI),
            se_pi = ifelse(n > 1, sd_pi / sqrt(n), 0), median_pi = median(PI), .groups = "drop")

# summary for plotting (excludes EG_35)
summary_stats_plot <- df_plot %>%
  group_by(POP) %>%
  summarise(n = n(), mean_pi = mean(PI), sd_pi = sd(PI),
            se_pi = ifelse(n > 1, sd_pi / sqrt(n), 0), median_pi = median(PI), .groups = "drop")

# write summaries
write_csv(summary_stats_analysis, file.path(outdir, "pi_window_summary_stats_analysis_including_EG35_remapped.csv"))
write_csv(summary_stats_plot, file.path(outdir, "pi_window_summary_stats_for_plots_BKN_HYD_only.csv"))
message("[INFO] Wrote summary CSVs for analysis and plotting to: ", outdir)

# PLOT 1: Boxplot + beeswarm/jitter (plotting df only)
p1 <- ggplot(df_plot, aes(x = POP, y = PI, fill = POP)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  (if (use_beeswarm)
    ggbeeswarm::geom_beeswarm(aes(color = POP), size = 1.2, alpha = 0.75)
   else
     geom_jitter(aes(color = POP), width = 0.15, alpha = 0.7, size = 1.5)) +
  scale_fill_manual(values = pop_colors) +
  scale_color_manual(values = pop_colors) +
  labs(title = "Windowed nucleotide diversity (π) by population",
       subtitle = "",
       x = "Population", y = expression(pi~"(per window)")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(outdir, "pi_boxplot_windowed.png"), p1, width = 10, height = 6, dpi = DPI)
message("Wrote: ", file.path(outdir, "pi_boxplot_windowed.png"))

# PLOT 2: Violin + mean ± SE (plotting df only)
p2 <- ggplot(df_plot, aes(x = POP, y = PI, fill = POP)) +
  geom_violin(trim = TRUE, alpha = 0.8) +
  geom_point(data = summary_stats_plot, mapping = aes(x = POP, y = mean_pi), inherit.aes = FALSE, color = "black", size = 2) +
  geom_errorbar(data = summary_stats_plot,
                mapping = aes(x = POP, ymin = mean_pi - se_pi, ymax = mean_pi + se_pi),
                inherit.aes = FALSE, width = 0.2) +
  scale_fill_manual(values = pop_colors) +
  labs(title = "Distribution of windowed π (violin) with mean ± SE",
       x = "Population", y = expression(pi~"(per window)")) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(outdir, "pi_violin_meanSE.png"), p2, width = 10, height = 6, dpi = DPI)
message("Wrote: ", file.path(outdir, "pi_violin_meanSE.png"))

# PLOT 3: Density overlays
p3 <- ggplot(df_plot, aes(x = PI, color = POP, fill = POP)) +
  geom_density(alpha = 0.25) +
  scale_color_manual(values = pop_colors) +
  scale_fill_manual(values = pop_colors) +
  labs(title = "Density of windowed π by population", x = expression(pi~"(per window)")) +
  theme_minimal(base_size = 13)
ggsave(file.path(outdir, "pi_density.png"), p3, width = 10, height = 6, dpi = DPI)
message("Wrote: ", file.path(outdir, "pi_density.png"))

# PLOT 4: Per-contig tracks - choose top contigs in df_plot
top_contigs <- df_plot %>% count(CHROM, sort = TRUE) %>% slice_head(n = 6) %>% pull(CHROM)
if (length(top_contigs) > 0) {
  df_top <- df_plot %>% filter(CHROM %in% top_contigs) %>% arrange(CHROM, BIN_START) %>% mutate(bin_mid = (BIN_START + BIN_END)/2)
  p4 <- ggplot(df_top, aes(x = bin_mid/1e3, y = PI, color = POP)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ CHROM, scales = "free_x", ncol = 1) +
    scale_color_manual(values = pop_colors) +
    labs(x = "Position (kb)", y = expression(pi~"(window)"),
         title = "Windowed π along top contigs (per population)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  ggsave(file.path(outdir, "pi_top_contigs_tracks.png"), p4, width = 10, height = max(6, 2 * length(top_contigs)), dpi = DPI)
  message("Wrote: ", file.path(outdir, "pi_top_contigs_tracks.png"))
} else {
  message("[INFO] No contigs found to plot tracks.")
}

# PLOT 5: Mean π barplot with error bars
p5 <- ggplot(summary_stats_plot, aes(x = POP, y = mean_pi, fill = POP)) +
  geom_col(width = 0.6, alpha = 0.9) +
  geom_errorbar(data = summary_stats_plot, mapping = aes(x = POP, ymin = mean_pi - se_pi, ymax = mean_pi + se_pi), inherit.aes = FALSE, width = 0.15) +
  scale_fill_manual(values = pop_colors) +
  labs(x = "Population", y = "Mean π (per window)", title = "Mean windowed π by population (±SE)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(outdir, "pi_mean_barplot.png"), p5, width = 9, height = 6, dpi = DPI)
message("Wrote: ", file.path(outdir, "pi_mean_barplot.png"))

# Save final CSVs
write_csv(summary_stats_analysis, file.path(outdir, "pi_window_summary_stats_analysis_including_EG35_remapped.csv"))
write_csv(summary_stats_plot, file.path(outdir, "pi_window_summary_stats_for_plots_BKN_HYD_only.csv"))
write_csv(df_plot, file.path(outdir, "pi_merged_windowed_for_plotting_BKN_HYD_only.csv"))
message("[INFO] Wrote all summary and plotting CSVs to: ", outdir)

message("Done. EG_35 was remapped into BKN for analysis (tagged as remapped_sample = 'BKN_UV_o_01') but excluded from plotting panels (plots only show genuine BKN & HYD).")
