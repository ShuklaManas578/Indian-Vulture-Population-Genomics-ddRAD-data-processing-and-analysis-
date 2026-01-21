#!/usr/bin/env Rscript


# plot_tajima.R (Rscript)
# Make boxplots, violin plots and genome-wide track of Tajima's D per set (ALL + populations)

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(ggplot2); library(readr); library(stringr); library(scales)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", default = "tajima_sliding_results/tajima_sliding_allsets_combined.csv",
              help="combined CSV produced by sliding_tajima_windows.sh or compute_tajimas_vcftools.sh"),
  make_option(c("-o","--outdir"), type="character", default = "plots_tajima", help="output directory for plots"),
  make_option(c("--dpi"), type="integer", default = 1200, help="PNG DPI"),
  make_option(c("--min_variants"), type="integer", default = 5, help="min variants in window to include in plots")
)
opt <- parse_args(OptionParser(option_list=option_list))
infile <- opt$input; OUTDIR <- opt$outdir; DPI <- opt$dpi; MIN_VARIANTS <- opt$min_variants
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# Check if input file exists
if (!file.exists(infile)) {
  stop("Input file not found: ", infile)
}

message("Reading data from: ", infile)
dat <- read_csv(infile, col_types = cols(
  chrom = col_character(), start = col_double(), end = col_double(),
  n_variants = col_double(), tajimaD = col_double(), set = col_character()
), progress = FALSE)

message("Data dimensions: ", nrow(dat), " rows x ", ncol(dat), " columns")
message("Sets found: ", paste(unique(dat$set), collapse = ", "))

# drop windows with extremely low variant counts or NA tajima
dat2 <- dat %>% filter(n_variants >= MIN_VARIANTS & !is.na(tajimaD))
message("After filtering (min variants >= ", MIN_VARIANTS, "): ", nrow(dat2), " windows")

# summary stats per set
summ <- dat2 %>% group_by(set) %>% summarise(
  n_windows = n(), 
  taj_mean = mean(tajimaD, na.rm=TRUE), 
  taj_median = median(tajimaD, na.rm=TRUE), 
  taj_sd = sd(tajimaD, na.rm=TRUE)
) %>% arrange(set)

write_csv(summ, file.path(OUTDIR, "tajima_summary_per_set.csv"))
message("Summary statistics written to: ", file.path(OUTDIR, "tajima_summary_per_set.csv"))

# Boxplot (per set)
message("Creating boxplot/violin plot...")
p1 <- ggplot(dat2, aes(x = set, y = tajimaD, fill = set)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.8, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.6) +
  labs(x = "Set (population)", y = "Tajima's D", title = "Tajima's D distribution per set (sliding windows)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1))

ggsave(file.path(OUTDIR, "tajima_box_violin_per_set.png"), p1, width = 10, height = 6, dpi = DPI)
ggsave(file.path(OUTDIR, "tajima_box_violin_per_set.svg"), p1, width = 10, height = 6)
message("Boxplot saved: ", file.path(OUTDIR, "tajima_box_violin_per_set.png"))

# genome-wide track: perturb the data to show per-chrom coordinates
# compute mid-point for plotting
dat2 <- dat2 %>% mutate(mid = (start + end)/2)

# make one plot per set (faceted) or combined colored by set
message("Creating genome-wide track plot...")
p2 <- ggplot(dat2, aes(x = mid/1e6, y = tajimaD, color = set)) +
  geom_point(size = 0.6, alpha = 0.7) +
  facet_wrap(~ set, ncol = 1, scales = "free_x") +
  labs(x = "Position (Mb)", y = "Tajima's D", title = "Genome-wide Tajima's D (per-window)") +
  theme_minimal(base_size = 12)

# Calculate appropriate height based on number of sets
n_sets <- length(unique(dat2$set))
plot_height <- max(3 * n_sets, 6)  # Minimum height of 6 inches

ggsave(file.path(OUTDIR, "tajima_genomewide_faceted.png"), p2, width = 12, height = plot_height, dpi = DPI)
ggsave(file.path(OUTDIR, "tajima_genomewide_faceted.svg"), p2, width = 12, height = plot_height)
message("Genome-wide track plot saved: ", file.path(OUTDIR, "tajima_genomewide_faceted.png"))

# per-contig summary (boxplot)
message("Creating per-contig summary plot...")
p3_data <- dat2 %>% 
  group_by(set, chrom) %>% 
  summarise(medianT = median(tajimaD, na.rm=TRUE), nwin = n(), .groups = 'drop')

p3 <- ggplot(p3_data, aes(x = reorder(chrom, -medianT), y = medianT, fill = set)) +
  geom_col(position = position_dodge(width=0.8)) +
  labs(x = "Contig (ordered by median TajimaD)", y = "median Tajima's D", title = "Per-contig median Tajima's D (per set)") +
  theme_minimal(base_size = 12) + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))

ggsave(file.path(OUTDIR, "tajima_per_contig_median.png"), p3, width = 12, height = 6, dpi = DPI)
ggsave(file.path(OUTDIR, "tajima_per_contig_median.svg"), p3, width = 12, height = 6)
message("Per-contig summary plot saved: ", file.path(OUTDIR, "tajima_per_contig_median.png"))

message("[DONE] All plots written to: ", normalizePath(OUTDIR))
message("[INFO] Summary CSV: ", file.path(OUTDIR, "tajima_summary_per_set.csv"))
