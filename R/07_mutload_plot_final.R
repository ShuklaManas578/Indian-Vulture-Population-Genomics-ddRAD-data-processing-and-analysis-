#!/usr/bin/env Rscript

# 15 Mutation_load Plots (Rscript)
# mutload_plot_final.R
# Produces mutation-load plots with error bar representation


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(ggridges)
  library(scales)
  library(tibble)
})

# parse args
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i + 1 > length(args)) return(default)
  return(args[i + 1])
}
summary_f <- getarg("--summary")
boot_f    <- getarg("--bootstrap")
master_f  <- getarg("--master", NA)
outdir    <- getarg("--outdir", "plots_final")
dpi_arg   <- getarg("--dpi", "1500")

if (is.null(summary_f) || is.null(boot_f)) {
  stop("Usage: Rscript mutload_plot_final_corrected.R --summary <file> --bootstrap <file> [--master <file>] --outdir <dir> --dpi <n>")
}
dpi <- as.integer(dpi_arg); if (is.na(dpi) || dpi <= 0) dpi <- 1500
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(outdir, "plot_log.txt")

cat(sprintf("plot_log: %s\n", Sys.time()), file = logf)
cat(sprintf("summary: %s\nbootstrap: %s\nmaster: %s\noutdir: %s\ndpi: %s\n\n",
            summary_f, boot_f, master_f, outdir, dpi), file = logf, append = TRUE)

# read tables
cat("Reading summary and bootstrap tables...\n", file = logf, append = TRUE)
sumdf <- fread(summary_f, data.table = FALSE)
bootdf <- fread(boot_f, data.table = FALSE)

# Clean column names carefully
clean_names <- function(df) {
  colnames(df) <- gsub("^X\\.\\.", "", colnames(df))
  colnames(df) <- gsub("^\\.\\.", "", colnames(df))
  colnames(df) <- gsub("\\.\\.", "\\.", colnames(df))
  colnames(df) <- gsub("^X", "", colnames(df))
  colnames(df) <- gsub("^\\.", "", colnames(df))
  return(df)
}

sumdf <- clean_names(sumdf)
bootdf <- clean_names(bootdf)

# Check column names
cat("Columns in summary:", paste(colnames(sumdf), collapse=", "), "\n", file=logf, append=TRUE)
cat("Columns in bootstrap:", paste(colnames(bootdf), collapse=", "), "\n", file=logf, append=TRUE)

# Ensure required columns exist and are properly named
if (!"sample" %in% colnames(sumdf)) {
  stop("No 'sample' column found in summary file. Available columns: ", paste(colnames(sumdf), collapse=", "))
}

# Ensure numeric types for key columns - check different possible column names
numeric_cols <- c("total_nonref", "total_hom", "total_nonref_perMb", 
                  "callable_bp")
impact_cols <- c("hom_HIGH", "hom_MODERATE", "hom_LOW", "hom_MODIFIER")

# Check which impact columns actually exist
existing_impact_cols <- impact_cols[impact_cols %in% colnames(sumdf)]
if (length(existing_impact_cols) < 4) {
  cat("Warning: Not all impact columns found. Found:", paste(existing_impact_cols, collapse=", "), "\n", file=logf, append=TRUE)
}

# Convert all numeric columns
for (col in c(numeric_cols, existing_impact_cols)) {
  if (col %in% colnames(sumdf)) {
    sumdf[[col]] <- as.numeric(as.character(sumdf[[col]]))
  }
}

# Read master table if provided
master_present <- FALSE
master_tbl <- NULL
if (!is.na(master_f) && nzchar(master_f) && file.exists(master_f)) {
  cat("Reading master table...\n", file = logf, append = TRUE)
  master_tbl <- tryCatch(fread(master_f, data.table = FALSE), 
                         error = function(e) {
                           cat("ERROR reading master table:", conditionMessage(e), "\n", 
                               file = logf, append = TRUE)
                           NULL
                         })
  if (!is.null(master_tbl) && nrow(master_tbl) > 0) {
    master_tbl <- clean_names(master_tbl)
    master_present <- TRUE
  }
}

# Prepare data for plotting 
# Merge summary and bootstrap data
sumdf <- left_join(sumdf, bootdf, by = "sample")

# Clean up column names from merge
colnames(sumdf) <- gsub("\\.x$", "", colnames(sumdf))
colnames(sumdf) <- gsub("\\.y$", "_bootstrap", colnames(sumdf))

# Calculate heterozygous derived alleles (total_nonref - total_hom)
sumdf$total_het <- sumdf$total_nonref - sumdf$total_hom

# Calculate homozygous mutation load per Mb
sumdf$total_hom_perMb <- (sumdf$total_hom / sumdf$callable_bp) * 1e6

# Calculate total impact sum for checking
impact_sum_cols <- c("hom_HIGH", "hom_MODERATE", "hom_LOW", "hom_MODIFIER")
if (all(impact_sum_cols %in% colnames(sumdf))) {
  sumdf$impact_sum <- rowSums(sumdf[, impact_sum_cols], na.rm = TRUE)
  # Check if impact_sum equals total_hom
  mismatch <- any(abs(sumdf$impact_sum - sumdf$total_hom) > 1, na.rm = TRUE)
  if (mismatch) {
    cat("Warning: Sum of impact categories does not equal total_hom for some samples\n", file=logf, append=TRUE)
  }
}

# Order samples by location and total mutation load
sumdf <- sumdf %>%
  arrange(location, desc(total_nonref_perMb))
sumdf$sample <- factor(sumdf$sample, levels = unique(sumdf$sample))
sumdf$location <- factor(sumdf$location, levels = c("BKN", "HYD"))

# Color palettes
impact_palette <- c(
  HIGH = "#DE3163",        # Cerise
  MODERATE = "#6E0B3A",    # Bordeaux
  LOW = "#D4AF37",         # Gold
  MODIFIER = "#006B6B"     # Teal
)

# For stacked bars
impact_palette_stacked <- c(
  MODIFIER = "#006B6B",
  LOW = "#D4AF37",
  MODERATE = "#6E0B3A",
  HIGH = "#DE3163"
)

# Palette without MODIFIER
impact_palette_noMOD <- impact_palette[names(impact_palette) != "MODIFIER"]
impact_palette_stacked_noMOD <- impact_palette_stacked[names(impact_palette_stacked) != "MODIFIER"]

derived_palette <- c(
  het = "#D95D39",         
  hom = "#9B2345"          
)

location_palette <- c(
  BKN = "#6E0B3A",        
  HYD = "#006B6B"        
)

# Sunset gradient for heatmaps and consequence plots
sunset_gradient <- colorRampPalette(c("#2D0421", "#5B0F47", "#9B2345", 
                                      "#D95D39", "#FF8C42", "#FFC857","#FFDDAA","#FFF3D6","#F7B7A3","#E07A5F"))(100)

sunset_discrete <- colorRampPalette(c("#2D0421", "#5B0F47", "#9B2345", 
                                      "#D95D39", "#FF8C42", "#FFC857","#FFDDAA","#FFF3D6","#F7B7A3","#E07A5F"))(10)

# Helper function to save plots
save_plot <- function(p, name, width = 10, height = 6, dpi = 1500) {
  pngfile <- file.path(outdir, paste0(name, ".png"))
  pdffile <- file.path(outdir, paste0(name, ".pdf"))
  
  png(filename = pngfile, width = width * dpi, height = height * dpi, res = dpi)
  print(p); dev.off()
  
  ggsave(pdffile, plot = p, width = width, height = height, device = "pdf")
  
  cat(sprintf("Saved: %s and %s\n", pngfile, pdffile), file = logf, append = TRUE)
}

# Create impact data frames
# Full impact data
impact_long <- sumdf %>%
  select(sample, location, any_of(impact_sum_cols)) %>%
  pivot_longer(cols = any_of(impact_sum_cols),
               names_to = "impact", values_to = "count") %>%
  mutate(impact = gsub("hom_", "", impact),
         impact = factor(impact, levels = c("MODIFIER", "LOW", "MODERATE", "HIGH")),
         count = as.numeric(count))

# Impact data without MODIFIER
impact_long_noMOD <- impact_long %>%
  filter(impact != "MODIFIER") %>%
  mutate(impact = factor(impact, levels = c("LOW", "MODERATE", "HIGH")))

# Plot 1: Total mutation load (derived alleles per Mb)
p1_data <- sumdf
# Ensure CI columns are present
ci_cols <- c("total_nonref_perMb_ci2.5", "total_nonref_perMb_ci97.5")
if (all(ci_cols %in% colnames(p1_data))) {
  p1_data$ci_lower <- as.numeric(p1_data$total_nonref_perMb_ci2.5)
  p1_data$ci_upper <- as.numeric(p1_data$total_nonref_perMb_ci97.5)
} else {
  cat("Warning: Missing CI columns for mutation load plot\n", file=logf, append=TRUE)
  p1_data$ci_lower <- NA
  p1_data$ci_upper <- NA
}

p1 <- ggplot(p1_data, aes(x = sample, y = total_nonref_perMb, fill = location)) +
  geom_col(width = 0.7) +
  {if(!all(is.na(p1_data$ci_lower))) 
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  width = 0.3, color = "black", linewidth = 0.5, na.rm = TRUE)} +
  scale_fill_manual(values = location_palette, 
                    breaks = names(location_palette),
                    labels = names(location_palette)) +
  theme_minimal(base_size = 11) +
  labs(title = "Total Mutation Load: Derived Alleles per Megabase",
       subtitle = "Error bars show 95% bootstrap CI",
       x = "Sample", y = "Derived Alleles per Mb", fill = "Location") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 10)) +
  coord_cartesian(ylim = c(0, max(p1_data$total_nonref_perMb, na.rm = TRUE) * 1.2))

save_plot(p1, "01_total_mutation_load_perMb_with_CI", width = 12, height = 6)

# Plot 2: Homozygous-derived counts by impact category (ALL)
p2 <- ggplot(impact_long, aes(x = sample, y = count, fill = impact)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = impact_palette_stacked,
                    breaks = names(impact_palette_stacked),
                    labels = names(impact_palette_stacked)) +
  theme_minimal(base_size = 11) +
  labs(title = "Homozygous-Derived Alleles by Functional Impact",
       subtitle = "Stacked by functional consequence",
       x = "Sample", y = "Homozygous-Derived Count", fill = "Impact") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

save_plot(p2, "02_homozygous_derived_by_impact_ALL", width = 12, height = 6)

# Plot 2a: Homozygous-derived counts by impact category (without MODIFIER)
p2a <- ggplot(impact_long_noMOD, aes(x = sample, y = count, fill = impact)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = impact_palette_stacked_noMOD,
                    breaks = names(impact_palette_stacked_noMOD),
                    labels = names(impact_palette_stacked_noMOD)) +
  theme_minimal(base_size = 11) +
  labs(title = "Homozygous-Derived Alleles by Functional Impact (excluding MODIFIER)",
       subtitle = "Stacked by functional consequence",
       x = "Sample", y = "Homozygous-Derived Count", fill = "Impact") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

save_plot(p2a, "02a_homozygous_derived_by_impact_noMOD", width = 12, height = 6)

# Plot 3: Proportion of homozygous-derived alleles by impact (ALL)
p3 <- ggplot(impact_long, aes(x = sample, y = count, fill = impact)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = impact_palette_stacked,
                    breaks = names(impact_palette_stacked),
                    labels = names(impact_palette_stacked)) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal(base_size = 11) +
  labs(title = "Proportion of Homozygous-Derived Alleles by Impact",
       subtitle = "Relative contribution of each functional category",
       x = "Sample", y = "Proportion", fill = "Impact") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

save_plot(p3, "03_homozygous_impact_proportions_ALL", width = 12, height = 6)

# Plot 3a: Proportion of homozygous-derived alleles by impact (without MODIFIER)
p3a <- ggplot(impact_long_noMOD, aes(x = sample, y = count, fill = impact)) +
  geom_col(position = "fill", width = 0.7) +
  scale_fill_manual(values = impact_palette_stacked_noMOD,
                    breaks = names(impact_palette_stacked_noMOD),
                    labels = names(impact_palette_stacked_noMOD)) +
  scale_y_continuous(labels = percent_format()) +
  theme_minimal(base_size = 11) +
  labs(title = "Proportion of Homozygous-Derived Alleles by Impact (excluding MODIFIER)",
       subtitle = "Relative contribution of functional categories",
       x = "Sample", y = "Proportion", fill = "Impact") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

save_plot(p3a, "03a_homozygous_impact_proportions_noMOD", width = 12, height = 6)

# Plot 4: Homozygous vs heterozygous derived alleles
if ("total_het" %in% colnames(sumdf) && "total_hom" %in% colnames(sumdf)) {
  geno_long <- sumdf %>%
    select(sample, location, total_hom, total_het) %>%
    pivot_longer(cols = c(total_hom, total_het),
                 names_to = "genotype", values_to = "count") %>%
    mutate(genotype = gsub("total_", "", genotype),
           genotype = factor(genotype, levels = c("het", "hom")),
           count = as.numeric(count))
  
  p4 <- ggplot(geno_long, aes(x = sample, y = count, fill = genotype)) +
    geom_col(position = "stack", width = 0.7) +
    scale_fill_manual(values = derived_palette,
                      breaks = names(derived_palette),
                      labels = c("Heterozygous", "Homozygous")) +
    theme_minimal(base_size = 11) +
    labs(title = "Derived Allele Composition: Homozygous vs Heterozygous",
         subtitle = "Stacked by zygosity",
         x = "Sample", y = "Derived Allele Count", fill = "Zygosity") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 12))
  
  save_plot(p4, "04_derived_allele_zygosity", width = 12, height = 6)
} else {
  cat("Warning: Missing total_het or total_hom columns for zygosity plot\n", file=logf, append=TRUE)
}

# Plot 5: Comparison between locations (mutation load)
p5a <- ggplot(sumdf, aes(x = location, y = total_nonref_perMb, fill = location)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  scale_fill_manual(values = location_palette,
                    breaks = names(location_palette),
                    labels = names(location_palette)) +
  theme_minimal(base_size = 12) +
  labs(title = "Mutation Load Comparison: BKN vs HYD",
       subtitle = "",
       x = "Location", y = "Derived Alleles per Mb", fill = "Location") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

save_plot(p5a, "05a_location_comparison_mutation_load", width = 8, height = 6)

# Plot 5b: Impact category comparison between locations (ALL)
if (nrow(impact_long) > 0) {
  impact_avg <- impact_long %>%
    group_by(location, impact) %>%
    summarise(mean_count = mean(count, na.rm = TRUE),
              se = sd(count, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(upper = mean_count + 1.96 * se,
           lower = pmax(mean_count - 1.96 * se, 0))
  
  # Reorder impact for plotting
  impact_avg$impact <- factor(impact_avg$impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
  
  p5b <- ggplot(impact_avg, aes(x = impact, y = mean_count, fill = location)) +
    geom_col(position = position_dodge(0.9), width = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(0.9), width = 0.3) +
    scale_fill_manual(values = location_palette,
                      breaks = names(location_palette),
                      labels = names(location_palette)) +
    theme_minimal(base_size = 12) +
    labs(title = "Average Homozygous-Derived Counts by Impact and Location",
         subtitle = "Bars show mean ± 95% CI",
         x = "Functional Impact", y = "Mean Homozygous-Derived Count",
         fill = "Location") +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  save_plot(p5b, "05b_impact_comparison_by_location_ALL", width = 10, height = 6)
}

# Plot 5c: Impact category comparison between locations (without MODIFIER)
if (nrow(impact_long_noMOD) > 0) {
  impact_avg_noMOD <- impact_long_noMOD %>%
    group_by(location, impact) %>%
    summarise(mean_count = mean(count, na.rm = TRUE),
              se = sd(count, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(upper = mean_count + 1.96 * se,
           lower = pmax(mean_count - 1.96 * se, 0))
  
  # Reorder impact for plotting
  impact_avg_noMOD$impact <- factor(impact_avg_noMOD$impact, levels = c("HIGH", "MODERATE", "LOW"))
  
  p5c <- ggplot(impact_avg_noMOD, aes(x = impact, y = mean_count, fill = location)) +
    geom_col(position = position_dodge(0.9), width = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = position_dodge(0.9), width = 0.3) +
    scale_fill_manual(values = location_palette,
                      breaks = names(location_palette),
                      labels = names(location_palette)) +
    theme_minimal(base_size = 12) +
    labs(title = "Average Homozygous-Derived Counts by Impact and Location (excluding MODIFIER)",
         subtitle = "Bars show mean ± 95% CI",
         x = "Functional Impact", y = "Mean Homozygous-Derived Count",
         fill = "Location") +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  save_plot(p5c, "05c_impact_comparison_by_location_noMOD", width = 10, height = 6)
}

# Plot 6: Correlation between mutation load metrics
if ("total_nonref" %in% colnames(sumdf) && "total_hom" %in% colnames(sumdf)) {
  p6 <- ggplot(sumdf, aes(x = total_nonref, y = total_hom, color = location)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, alpha = 0.3) +
    geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20, 
                    box.padding = 0.5, point.padding = 0.2) +
    scale_color_manual(values = location_palette,
                       breaks = names(location_palette),
                       labels = names(location_palette)) +
    theme_minimal(base_size = 12) +
    labs(title = "Correlation: Total Derived vs Homozygous-Derived Alleles",
         subtitle = "Linear relationship between mutation load metrics",
         x = "Total Derived Alleles", y = "Homozygous-Derived Alleles",
         color = "Location") +
    theme(plot.title = element_text(face = "bold"))
  
  save_plot(p6, "06_correlation_total_vs_homozygous", width = 9, height = 7)
}

# Plot 7: Heatmap of samples by impact categories (ALL)
if (all(impact_sum_cols %in% colnames(sumdf))) {
  heatmap_data <- sumdf %>%
    select(sample, all_of(impact_sum_cols)) %>%
    as.data.frame()
  
  # Set row names and convert to matrix
  rownames(heatmap_data) <- heatmap_data$sample
  heatmap_data$sample <- NULL
  heatmap_mat <- as.matrix(heatmap_data)
  
  # Normalize by row for better visualization
  heatmap_norm <- t(apply(heatmap_mat, 1, function(x) x/sum(x)))
  
  # Create the heatmap
  png(file.path(outdir, "07_heatmap_impact_profiles_ALL.png"), 
      width = 10 * dpi, height = 8 * dpi, res = dpi)
  pheatmap(heatmap_norm,
           color = sunset_gradient,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           main = "Sample Impact Profiles: Normalized Homozygous-Derived Counts",
           fontsize_row = 8,
           fontsize_col = 10)
  dev.off()
  
  # Also save as PDF
  pdf(file.path(outdir, "07_heatmap_impact_profiles_ALL.pdf"), width = 10, height = 8)
  pheatmap(heatmap_norm,
           color = sunset_gradient,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           main = "Sample Impact Profiles: Normalized Homozygous-Derived Counts",
           fontsize_row = 8,
           fontsize_col = 10)
  dev.off()
  
  cat("Saved heatmap: 07_heatmap_impact_profiles_ALL.png and .pdf\n", file=logf, append=TRUE)
}

# Plot 7a: Heatmap of samples by impact categories (without MODIFIER)
if (all(c("hom_HIGH", "hom_MODERATE", "hom_LOW") %in% colnames(sumdf))) {
  heatmap_data_noMOD <- sumdf %>%
    select(sample, hom_HIGH, hom_MODERATE, hom_LOW) %>%
    as.data.frame()
  
  # Set row names and convert to matrix
  rownames(heatmap_data_noMOD) <- heatmap_data_noMOD$sample
  heatmap_data_noMOD$sample <- NULL
  heatmap_mat_noMOD <- as.matrix(heatmap_data_noMOD)
  
  # Normalize by row for better visualization
  heatmap_norm_noMOD <- t(apply(heatmap_mat_noMOD, 1, function(x) x/sum(x)))
  
  # Create the heatmap
  png(file.path(outdir, "07a_heatmap_impact_profiles_noMOD.png"), 
      width = 10 * dpi, height = 8 * dpi, res = dpi)
  pheatmap(heatmap_norm_noMOD,
           color = sunset_gradient,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           main = "Sample Impact Profiles (excluding MODIFIER): Normalized Homozygous-Derived Counts",
           fontsize_row = 8,
           fontsize_col = 10)
  dev.off()
  
  # Also save as PDF
  pdf(file.path(outdir, "07a_heatmap_impact_profiles_noMOD.pdf"), width = 10, height = 8)
  pheatmap(heatmap_norm_noMOD,
           color = sunset_gradient,
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           main = "Sample Impact Profiles (excluding MODIFIER): Normalized Homozygous-Derived Counts",
           fontsize_row = 8,
           fontsize_col = 10)
  dev.off()
  
  cat("Saved heatmap: 07a_heatmap_impact_profiles_noMOD.png and .pdf\n", file=logf, append=TRUE)
}

# Plot 8: Distribution of mutation loads
p8 <- ggplot(sumdf, aes(x = total_nonref_perMb, y = location, fill = location)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9) +
  scale_fill_manual(values = location_palette,
                    breaks = names(location_palette),
                    labels = names(location_palette)) +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Mutation Load per Location",
       subtitle = "",
       x = "Derived Alleles per Mb", y = "Location", fill = "Location") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

save_plot(p8, "08_mutation_load_distribution", width = 9, height = 6)

# Plot 9: Bootstrap confidence intervals for impact categories (ALL)
boot_impact_cols <- c("HIGH_hom_ci2.5", "HIGH_hom_ci97.5",
                      "MODERATE_hom_ci2.5", "MODERATE_hom_ci97.5",
                      "LOW_hom_ci2.5", "LOW_hom_ci97.5",
                      "MODIFIER_hom_ci2.5", "MODIFIER_hom_ci97.5")

if (all(boot_impact_cols %in% colnames(sumdf))) {
  boot_impact_data <- data.frame()
  for (impact in c("HIGH", "MODERATE", "LOW", "MODIFIER")) {
    lower_col <- paste0(impact, "_hom_ci2.5")
    upper_col <- paste0(impact, "_hom_ci97.5")
    mean_col <- paste0("hom_", impact)
    
    if (all(c(lower_col, upper_col, mean_col) %in% colnames(sumdf))) {
      temp <- sumdf %>%
        select(sample, location, all_of(c(lower_col, upper_col, mean_col))) %>%
        rename(lower = !!sym(lower_col), upper = !!sym(upper_col), mean = !!sym(mean_col)) %>%
        mutate(impact = impact)
      
      boot_impact_data <- rbind(boot_impact_data, temp)
    }
  }
  
  if (nrow(boot_impact_data) > 0) {
    boot_impact_data$impact <- factor(boot_impact_data$impact, 
                                      levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"))
    
    p9 <- ggplot(boot_impact_data, aes(x = sample, y = mean, color = impact)) +
      geom_point(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    width = 0.3, position = position_dodge(width = 0.5)) +
      facet_wrap(~ impact, scales = "free_y", ncol = 2) +
      scale_color_manual(values = impact_palette,
                         breaks = names(impact_palette),
                         labels = names(impact_palette)) +
      theme_minimal(base_size = 11) +
      labs(title = "Bootstrap Confidence Intervals for Homozygous-Derived Counts by Impact",
           subtitle = "Error bars show 95% bootstrap CI",
           x = "Sample", y = "Homozygous-Derived Count", color = "Impact") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(face = "bold", size = 12),
            strip.text = element_text(face = "bold"))
    
    save_plot(p9, "09_bootstrap_CI_by_impact_ALL", width = 14, height = 10)
  }
}

# Plot 9a: Bootstrap confidence intervals for impact categories (without MODIFIER)
boot_impact_cols_noMOD <- c("HIGH_hom_ci2.5", "HIGH_hom_ci97.5",
                            "MODERATE_hom_ci2.5", "MODERATE_hom_ci97.5",
                            "LOW_hom_ci2.5", "LOW_hom_ci97.5")

if (all(boot_impact_cols_noMOD %in% colnames(sumdf))) {
  boot_impact_data_noMOD <- data.frame()
  for (impact in c("HIGH", "MODERATE", "LOW")) {
    lower_col <- paste0(impact, "_hom_ci2.5")
    upper_col <- paste0(impact, "_hom_ci97.5")
    mean_col <- paste0("hom_", impact)
    
    if (all(c(lower_col, upper_col, mean_col) %in% colnames(sumdf))) {
      temp <- sumdf %>%
        select(sample, location, all_of(c(lower_col, upper_col, mean_col))) %>%
        rename(lower = !!sym(lower_col), upper = !!sym(upper_col), mean = !!sym(mean_col)) %>%
        mutate(impact = impact)
      
      boot_impact_data_noMOD <- rbind(boot_impact_data_noMOD, temp)
    }
  }
  
  if (nrow(boot_impact_data_noMOD) > 0) {
    boot_impact_data_noMOD$impact <- factor(boot_impact_data_noMOD$impact, 
                                            levels = c("HIGH", "MODERATE", "LOW"))
    
    p9a <- ggplot(boot_impact_data_noMOD, aes(x = sample, y = mean, color = impact)) +
      geom_point(size = 2, position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    width = 0.3, position = position_dodge(width = 0.5)) +
      facet_wrap(~ impact, scales = "free_y", ncol = 2) +
      scale_color_manual(values = impact_palette_noMOD,
                         breaks = names(impact_palette_noMOD),
                         labels = names(impact_palette_noMOD)) +
      theme_minimal(base_size = 11) +
      labs(title = "Bootstrap Confidence Intervals for Homozygous-Derived Counts by Impact (excluding MODIFIER)",
           subtitle = "Error bars show 95% bootstrap CI",
           x = "Sample", y = "Homozygous-Derived Count", color = "Impact") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(face = "bold", size = 12),
            strip.text = element_text(face = "bold"))
    
    save_plot(p9a, "09a_bootstrap_CI_by_impact_noMOD", width = 14, height = 10)
  }
}

# Plot 10: Master table analyses (if available)
if (master_present) {
  cat("Analyzing master table data...\n", file = logf, append = TRUE)
  
  # Check for required columns
  required_master_cols <- c("CONSEQUENCE", "GENE", "WORST_IMPACT")
  if (all(required_master_cols %in% colnames(master_tbl))) {
    
    # Extract derived allele columns
    derived_cols <- grep("^DERIVED_", colnames(master_tbl), value = TRUE)
    
    if (length(derived_cols) > 0) {
      cat(sprintf("Found %d derived allele columns\n", length(derived_cols)), file=logf, append=TRUE)
      
      # Create long format for consequences
      cons_data <- master_tbl %>%
        select(CONSEQUENCE, WORST_IMPACT, GENE, all_of(derived_cols)) %>%
        pivot_longer(cols = all_of(derived_cols),
                     names_to = "sample",
                     values_to = "derived_count") %>%
        mutate(sample = gsub("^DERIVED_", "", sample),
               derived_count = suppressWarnings(as.numeric(derived_count))) %>%
        # Remove rows where derived_count is NA or 0
        filter(!is.na(derived_count), derived_count > 0)
      
      if (nrow(cons_data) > 0) {
        # Summarize top consequences
        top_cons <- cons_data %>%
          group_by(CONSEQUENCE) %>%
          summarise(total_derived = sum(derived_count), .groups = "drop") %>%
          arrange(desc(total_derived)) %>%
          slice_head(n = 10)
        
        # Use sunset discrete palette for consequences
        top_cons$CONSEQUENCE <- factor(top_cons$CONSEQUENCE, 
                                       levels = top_cons$CONSEQUENCE[order(top_cons$total_derived, decreasing = TRUE)])
        
        p10a <- ggplot(top_cons, aes(x = reorder(CONSEQUENCE, total_derived), 
                                     y = total_derived, fill = CONSEQUENCE)) +
          geom_col() +
          coord_flip() +
          scale_fill_manual(values = sunset_discrete,
                            breaks = levels(top_cons$CONSEQUENCE),
                            labels = levels(top_cons$CONSEQUENCE)) +
          theme_minimal(base_size = 11) +
          labs(title = "Top 10 Consequence Categories by Derived Allele Count",
               subtitle = "",
               x = "Consequence", y = "Total Derived Alleles", fill = "Consequence") +
          theme(plot.title = element_text(face = "bold"),
                legend.position = "none") 
        
        save_plot(p10a, "10a_top_consequences_sunset", width = 10, height = 6)
        
        # Also create version with individual colors for each consequence (not gradient)
        p10a_individual <- ggplot(top_cons, aes(x = reorder(CONSEQUENCE, total_derived), 
                                                y = total_derived)) +
          geom_col(fill = "#9B2345") +
          coord_flip() +
          theme_minimal(base_size = 11) +
          labs(title = "Top 10 Consequence Categories by Derived Allele Count",
               subtitle = "",
               x = "Consequence", y = "Total Derived Alleles") +
          theme(plot.title = element_text(face = "bold"))
        
        save_plot(p10a_individual, "10a_top_consequences_single", width = 10, height = 6)
        
        # Consequence profiles by sample
        cons_sample <- cons_data %>%
          filter(CONSEQUENCE %in% top_cons$CONSEQUENCE) %>%
          group_by(sample, CONSEQUENCE) %>%
          summarise(total = sum(derived_count), .groups = "drop") %>%
          complete(sample, CONSEQUENCE, fill = list(total = 0))
        
        # Ensure samples are in the same order as before
        cons_sample$sample <- factor(cons_sample$sample, levels = levels(sumdf$sample))
        cons_sample$CONSEQUENCE <- factor(cons_sample$CONSEQUENCE, 
                                          levels = levels(top_cons$CONSEQUENCE))
        
        p10b <- ggplot(cons_sample, aes(x = sample, y = total, fill = CONSEQUENCE)) +
          geom_col(position = "fill") +
          scale_y_continuous(labels = percent_format()) +
          scale_fill_manual(values = sunset_discrete,
                            breaks = levels(cons_sample$CONSEQUENCE),
                            labels = levels(cons_sample$CONSEQUENCE)) +
          theme_minimal(base_size = 10) +
          labs(title = "Consequence Profiles by Sample",
               subtitle = "Relative proportions of top 10 consequence categories",
               x = "Sample", y = "Proportion", fill = "Consequence") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                plot.title = element_text(face = "bold"),
                legend.text = element_text(size = 8),
                legend.position = "right")
        
        save_plot(p10b, "10b_consequence_profiles_sunset", width = 12, height = 6)
        
        # Gene-level analysis
        gene_data <- cons_data %>%
          filter(!is.na(GENE), GENE != "", GENE != ".") %>%
          group_by(GENE, WORST_IMPACT) %>%
          summarise(total_derived = sum(derived_count),
                    num_samples = n_distinct(sample),
                    .groups = "drop") %>%
          arrange(desc(total_derived))
        
        if (nrow(gene_data) > 0) {
          top_genes <- gene_data %>%
            slice_head(n = 20)
          
          # Create color palette for impacts
          gene_impact_palette <- impact_palette[names(impact_palette) %in% unique(top_genes$WORST_IMPACT)]
          
          p10c <- ggplot(top_genes, aes(x = reorder(GENE, total_derived), 
                                        y = total_derived, fill = WORST_IMPACT)) +
            geom_col() +
            coord_flip() +
            scale_fill_manual(values = gene_impact_palette,
                              breaks = names(gene_impact_palette),
                              labels = names(gene_impact_palette)) +
            theme_minimal(base_size = 11) +
            labs(title = "Top 20 Genes by Derived Allele Burden",
                 subtitle = "Colored by functional impact",
                 x = "Gene", y = "Total Derived Alleles", fill = "Impact") +
            theme(plot.title = element_text(face = "bold"))
          
          save_plot(p10c, "10c_top_genes_by_impact", width = 10, height = 8)
          
          # Also create gene plot with sunset gradient based on total derived count
          top_genes_simple <- gene_data %>%
            group_by(GENE) %>%
            summarise(total_derived = sum(total_derived),
                      .groups = "drop") %>%
            arrange(desc(total_derived)) %>%
            slice_head(n = 20)
          
          # Create a gradient color based on total_derived values
          p10d <- ggplot(top_genes_simple, aes(x = reorder(GENE, total_derived), 
                                               y = total_derived, fill = total_derived)) +
            geom_col() +
            coord_flip() +
            scale_fill_gradientn(colors = sunset_gradient,
                                 name = "Derived Allele Count") +
            theme_minimal(base_size = 11) +
            labs(title = "Top 20 Genes by Derived Allele Burden",
                 subtitle = "Colored by derived allele count",
                 x = "Gene", y = "Total Derived Alleles") +
            theme(plot.title = element_text(face = "bold"))
          
          save_plot(p10d, "10d_top_genes_gradient", width = 10, height = 8)
        }
      } else {
        cat("No valid consequence data found in master table\n", file=logf, append=TRUE)
      }
    } else {
      cat("No derived allele columns found in master table\n", file=logf, append=TRUE)
    }
  } else {
    missing_cols <- setdiff(required_master_cols, colnames(master_tbl))
    cat(sprintf("Missing required columns in master table: %s\n", 
                paste(missing_cols, collapse=", ")), file=logf, append=TRUE)
  }
} else {
  cat("No master table provided or master failed to read -> skipping consequence/gene plots\n", 
      file = logf, append = TRUE)
}

# Plot 11: Key metrics summary
# Calculate high impact ratio only if we have all necessary columns
required_for_ratio <- c("hom_HIGH", "hom_MODERATE", "hom_LOW", "hom_MODIFIER")
if (all(required_for_ratio %in% colnames(sumdf))) {
  key_metrics <- sumdf %>%
    select(sample, location, 
           total_nonref_perMb, total_hom_perMb,
           hom_HIGH, hom_MODERATE, hom_LOW, hom_MODIFIER) %>%
    mutate(high_impact_ratio = hom_HIGH / (hom_HIGH + hom_MODERATE + hom_LOW + hom_MODIFIER))
  
  # Create a combined metrics plot
  metrics_melt <- key_metrics %>%
    select(sample, location, total_nonref_perMb, total_hom_perMb, high_impact_ratio) %>%
    pivot_longer(cols = c(total_nonref_perMb, total_hom_perMb, high_impact_ratio),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, 
                           levels = c("total_nonref_perMb", "total_hom_perMb", "high_impact_ratio"),
                           labels = c("Total Derived/Mb", "Homozygous Derived/Mb", "High Impact Ratio")))
  
  p11 <- ggplot(metrics_melt, aes(x = location, y = value, fill = location)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = location_palette,
                      breaks = names(location_palette),
                      labels = names(location_palette)) +
    theme_minimal(base_size = 11) +
    labs(title = "Key Mutation Load Metrics by Location",
         x = "Location", y = "Value", fill = "Location") +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 12),
          strip.text = element_text(face = "bold"),
          axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  save_plot(p11, "11_key_metrics_summary", width = 14, height = 5)
}

# Create summary statistics table
summary_stats <- sumdf %>%
  group_by(location) %>%
  summarise(
    n_samples = n(),
    mean_mutation_load = mean(total_nonref_perMb, na.rm = TRUE),
    sd_mutation_load = sd(total_nonref_perMb, na.rm = TRUE),
    median_mutation_load = median(total_nonref_perMb, na.rm = TRUE),
    mean_homozygous = mean(total_hom, na.rm = TRUE),
    sd_homozygous = sd(total_hom, na.rm = TRUE),
    median_homozygous = median(total_hom, na.rm = TRUE),
    mean_heterozygous = ifelse("total_het" %in% colnames(sumdf), mean(total_het, na.rm = TRUE), NA),
    sd_heterozygous = ifelse("total_het" %in% colnames(sumdf), sd(total_het, na.rm = TRUE), NA),
    median_heterozygous = ifelse("total_het" %in% colnames(sumdf), median(total_het, na.rm = TRUE), NA),
    .groups = "drop"
  )

write.csv(summary_stats, file.path(outdir, "location_summary_statistics.csv"), 
          row.names = FALSE, quote = FALSE)

# Save full data with CI
full_data_cols <- c("sample", "location", "callable_bp", "total_nonref", "total_hom", 
                    "total_nonref_perMb", "total_hom_perMb")
if ("total_het" %in% colnames(sumdf)) full_data_cols <- c(full_data_cols, "total_het")
full_data_cols <- c(full_data_cols, existing_impact_cols)

# Add bootstrap CI columns if they exist
ci_cols_to_add <- c("total_nonref_ci2.5", "total_nonref_ci97.5",
                    "total_nonref_perMb_ci2.5", "total_nonref_perMb_ci97.5",
                    "HIGH_hom_ci2.5", "HIGH_hom_ci97.5",
                    "MODERATE_hom_ci2.5", "MODERATE_hom_ci97.5",
                    "LOW_hom_ci2.5", "LOW_hom_ci97.5",
                    "MODIFIER_hom_ci2.5", "MODIFIER_hom_ci97.5")

full_data_cols <- c(full_data_cols, ci_cols_to_add[ci_cols_to_add %in% colnames(sumdf)])

full_data <- sumdf %>%
  select(any_of(full_data_cols))

write.csv(full_data, file.path(outdir, "full_sample_data_with_CI.csv"), 
          row.names = FALSE, quote = FALSE)

# Finalize
cat("\n=== PLOTTING SUMMARY ===\n", file = logf, append = TRUE)
cat(sprintf("Total samples plotted: %d\n", nrow(sumdf)), file = logf, append = TRUE)
cat(sprintf("Samples per location: BKN=%d, HYD=%d\n",
            sum(sumdf$location == "BKN"), sum(sumdf$location == "HYD")),
    file = logf, append = TRUE)

if (sum(sumdf$location == "BKN") > 0) {
  bkn_mean <- mean(sumdf$total_nonref_perMb[sumdf$location == "BKN"], na.rm = TRUE)
  cat(sprintf("BKN mean mutation load: %.1f derived alleles/Mb\n", bkn_mean), file=logf, append=TRUE)
}

if (sum(sumdf$location == "HYD") > 0) {
  hyd_mean <- mean(sumdf$total_nonref_perMb[sumdf$location == "HYD"], na.rm = TRUE)
  cat(sprintf("HYD mean mutation load: %.1f derived alleles/Mb\n", hyd_mean), file=logf, append=TRUE)
}

cat(sprintf("Total plots generated: %d\n", length(list.files(outdir, pattern="\\.(png|pdf)$"))), 
    file=logf, append=TRUE)
cat(sprintf("Plots saved to: %s\n", normalizePath(outdir)), file = logf, append = TRUE)
cat(sprintf("Log file: %s\n", normalizePath(logf)), file = logf, append = TRUE)

# Print summary to console
cat("\n=== MUTATION LOAD ANALYSIS COMPLETE ===\n")
cat(sprintf("Samples analyzed: %d (BKN: %d, HYD: %d)\n", 
            nrow(sumdf), sum(sumdf$location == "BKN"), sum(sumdf$location == "HYD")))
if (exists("bkn_mean") && exists("hyd_mean")) {
  cat(sprintf("Mean mutation load: BKN=%.1f, HYD=%.1f derived alleles/Mb\n", bkn_mean, hyd_mean))
}
cat(sprintf("Summary statistics saved to: %s\n", file.path(outdir, "location_summary_statistics.csv")))
cat(sprintf("Full data with CI saved to: %s\n", file.path(outdir, "full_sample_data_with_CI.csv")))
cat(sprintf("All plots saved to: %s\n", normalizePath(outdir)))
cat(sprintf("Log file: %s\n", normalizePath(logf)))
cat("\nPlot files created:\n")
plot_files <- list.files(outdir, pattern="\\.(png|pdf)$")
for (file in plot_files) {
  cat(sprintf("  - %s\n", file))
}
