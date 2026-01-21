#!/usr/bin/env Rscript


#8 RoH-based plots (Rscript)
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(viridis)   # cividis, magma, etc.
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
workdir <- if (length(args) >= 1) args[[1]] else "."
workdir <- normalizePath(workdir)
setwd(workdir)
cat("Working directory:", workdir, "\n")

# output folder for plots
outdir <- file.path(workdir, "roh_plots212")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat("Outputs will be written to:", outdir, "\n\n")

# find PLINK ROH .hom files
hom_files <- list.files(pattern = "_roh\\.hom$", ignore.case = TRUE)
if (length(hom_files) == 0) stop("No '*_roh.hom' files found in: ", workdir)

# robust reader (tries header; fallback to no-header heuristics)
read_roh_data <- function(filepath) {
  popname <- sub("_roh\\.hom$", "", basename(filepath), ignore.case = TRUE)
  df <- tryCatch(
    read.table(filepath, header = TRUE, stringsAsFactors = FALSE, comment.char = ""),
    error = function(e) {
      warning("Failed to read ", filepath, " with header=TRUE. Trying without header...")
      read.table(filepath, header = FALSE, stringsAsFactors = FALSE)
    }
  )
  # If headerless, attempt to assign sensible names
  if (!"IID" %in% names(df) && ncol(df) >= 2) {
    possible_names <- c("FID","IID","CHR","POS1","POS2","KB","NSNP","SNP_COUNT")
    names(df)[1:min(length(possible_names), ncol(df))] <- possible_names[1:min(length(possible_names), ncol(df))]
  }
  if (!"IID" %in% names(df)) stop("File ", filepath, " doesn't contain an IID column (and auto-recovery failed).")
  # compute KB if missing
  if (!"KB" %in% names(df)) {
    if (all(c("POS1","POS2") %in% names(df))) {
      df$KB <- (as.numeric(df$POS2) - as.numeric(df$POS1)) / 1000
      warning("Computed KB from POS1/POS2 for ", filepath)
    } else {
      stop("File ", filepath, " lacks KB and POS1/POS2 columns; cannot compute ROH length.")
    }
  } else {
    df$KB <- as.numeric(df$KB)
  }
  # remove non-positive
  df <- df %>% filter(!is.na(KB) & KB > 0)
  df$Population <- popname
  df$Individual <- df$IID
  return(df)
}

# read all files
cat("Reading files:\n")
roh_list <- lapply(hom_files, function(f) { cat("  ", f, "\n"); read_roh_data(f) })
all_roh_full <- bind_rows(roh_list)
if (nrow(all_roh_full) == 0) stop("No ROH segments found after reading files.")

# Colors: explicit values so remapped sample color is stable
col_HYD   <- "#DE3163"  # cerise
col_BKN   <- "#008B8B"  # dark cyan
col_EG35  <- "#B5A642"  # brass (the color to reuse for remapped BKN_UV_o_01)

# Define category levels and colors
category_levels <- c("0-500 kb\n(Ancient)",
                     "500-1000 kb\n(Old)",
                     "1000-2000 kb\n(Intermediate)",
                     "2000-5000 kb\n(Recent)",
                     ">5000 kb\n(Very Recent)")

cat_colors <- c(
  "0-500 kb\n(Ancient)" = "#ECECEC",       
  "500-1000 kb\n(Old)" = "#9AD0EA",        
  "1000-2000 kb\n(Intermediate)" = "#66C2A5", 
  "2000-5000 kb\n(Recent)" = "#008B8B",    
  ">5000 kb\n(Very Recent)" = "#D3AF37"    
)

# Use findInterval for more robust categorization
all_roh_full$ROH_Size_Category <- factor(
  findInterval(all_roh_full$KB, 
               c(0, 500, 1000, 2000, 5000, Inf),
               rightmost.closed = FALSE),
  levels = 1:5,
  labels = category_levels
)

# Remove any NA categories
all_roh_full <- all_roh_full %>% filter(!is.na(ROH_Size_Category))

# Ensure factor levels
all_roh_full$ROH_Size_Category <- factor(all_roh_full$ROH_Size_Category, levels = category_levels)

# Prepare filtered datasets

# 1) For "all other plots": use only genuine BKN and HYD (exclude any EG and exclude any BKN_o)
all_roh_filtered <- all_roh_full %>%
  filter(Population %in% c("BKN", "HYD"))

# 2) For the two special plots:

all_roh_remapped <- all_roh_full %>%
  mutate(
    remapped_flag = ifelse(Population == "EG_35", TRUE, FALSE),
    Population = ifelse(remapped_flag, "BKN_o", Population),
    Individual = ifelse(remapped_flag, "BKN_UV_o_01", Individual)
  )

# Ensure factor levels for size category in remapped frames
all_roh_remapped$ROH_Size_Category <- factor(all_roh_remapped$ROH_Size_Category, levels = category_levels)

# Prepare summary tables for "other plots" (BKN & HYD only)
individual_totals_filtered <- all_roh_filtered %>%
  group_by(Population, Individual) %>%
  summarise(Total_ROH_KB = sum(KB, na.rm = TRUE),
            Number_Segments = n(),
            Mean_Segment_KB = mean(KB, na.rm = TRUE),
            Max_Segment_KB = max(KB, na.rm = TRUE),
            .groups = "drop")

category_totals_filtered <- all_roh_filtered %>%
  group_by(Population, Individual, ROH_Size_Category) %>%
  summarise(Category_KB = sum(KB, na.rm = TRUE),
            Segments_in_Category = n(),
            .groups = "drop")

population_category_summary_filtered <- category_totals_filtered %>%
  group_by(Population, ROH_Size_Category) %>%
  summarise(Total_KB = sum(Category_KB, na.rm = TRUE),
            Mean_KB_per_Ind = mean(Category_KB, na.rm = TRUE),
            .groups = "drop")

# Prepare summary tables for remapped/included case (for special plots)
individual_totals_remapped <- all_roh_remapped %>%
  group_by(Population, Individual) %>%
  summarise(Total_ROH_KB = sum(KB, na.rm = TRUE),
            Number_Segments = n(),
            Mean_Segment_KB = mean(KB, na.rm = TRUE),
            Max_Segment_KB = max(KB, na.rm = TRUE),
            .groups = "drop")

category_totals_remapped <- all_roh_remapped %>%
  group_by(Population, Individual, ROH_Size_Category) %>%
  summarise(Category_KB = sum(KB, na.rm = TRUE),
            Segments_in_Category = n(),
            .groups = "drop")

population_category_summary_remapped <- category_totals_remapped %>%
  group_by(Population, ROH_Size_Category) %>%
  summarise(Total_KB = sum(Category_KB, na.rm = TRUE),
            Mean_KB_per_Ind = mean(Category_KB, na.rm = TRUE),
            .groups = "drop")

# For all "other" plots (only BKN & HYD)
pop_color_map_filtered <- c("BKN" = col_BKN, "HYD" = col_HYD)

# For remapped case: include BKN_o colored as original EG color
unique_pops_remapped <- unique(all_roh_remapped$Population)

remapped_color_map <- list()
for (p in unique_pops_remapped) {
  if (p == "BKN") remapped_color_map[[p]] <- col_BKN
  else if (p == "HYD") remapped_color_map[[p]] <- col_HYD
  else if (p == "BKN_o") remapped_color_map[[p]] <- col_EG35
  else {
    # assign remaining a viridis palette
    remapped_color_map[[p]] <- viridis::cividis(10)[(which(unique_pops_remapped == p) %% 10) + 1]
  }
}
remapped_color_map <- unlist(remapped_color_map)

# plotting params
dpi_val <- 1200
plot_bg <- "white"
base_size <- 14


# PLOT 1: Total ROH by Population (filtered: BKN & HYD only)

p1 <- ggplot(individual_totals_filtered, aes(x = Population, y = Total_ROH_KB, fill = Population)) +
  geom_boxplot(alpha = 0.85, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.85, size = 3, aes(color = Population)) +
  scale_fill_manual(values = pop_color_map_filtered) +
  scale_color_manual(values = pop_color_map_filtered) +
  labs(title = "Total ROH burden by population (BKN & HYD only)",
       subtitle = paste0("ROH segments (", sum(individual_totals_filtered$Number_Segments), " segments total for BKN & HYD)"),
       x = "Population", y = "Total ROH (KB)") +
  theme_minimal(base_size = base_size) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(filename = file.path(outdir, "custom_total_roh_by_population.png"),
       plot = p1, width = 10, height = 7, dpi = dpi_val, bg = plot_bg)


# PLOT 2: Proportional stacked bar by category (filtered)

p2 <- ggplot(population_category_summary_filtered,
             aes(x = Population, y = Total_KB, fill = ROH_Size_Category)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_fill_manual(values = cat_colors, drop = FALSE, na.value = NA) +
  labs(title = "ROH size category distribution by population",
       x = "Population", y = "Proportion of total ROH") +
  theme_minimal(base_size = base_size) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold")) +
  scale_y_continuous(labels = percent_format())

ggsave(filename = file.path(outdir, "custom_roh_categories_proportional.png"),
       plot = p2, width = 12, height = 7, dpi = dpi_val, bg = plot_bg)


# PLOT 3: Absolute stacked bars (filtered)

p3 <- ggplot(population_category_summary_filtered,
             aes(x = Population, y = Total_KB, fill = ROH_Size_Category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = cat_colors, drop = FALSE, na.value = NA) +
  labs(title = "Absolute ROH length by size category and population",
       x = "Population", y = "Total ROH (KB)") +
  theme_minimal(base_size = base_size) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"))

ggsave(filename = file.path(outdir, "custom_roh_categories_absolute.png"),
       plot = p3, width = 12, height = 7, dpi = dpi_val, bg = plot_bg)


# PLOT 4: Individual stacked profiles (two versions)

# Prepare ordering for excludeEG35
individual_order_excl <- individual_totals_filtered %>%
  arrange(Population, desc(Total_ROH_KB)) %>%
  group_by(Population) %>%
  mutate(Order = row_number()) %>%
  ungroup() %>%
  select(Population, Individual, Total_ROH_KB, Order)

category_totals_excl2 <- category_totals_filtered %>%
  left_join(individual_order_excl, by = c("Population", "Individual")) %>%
  arrange(Population, desc(Total_ROH_KB))

max_order_excl <- ifelse(nrow(individual_order_excl) > 0, max(individual_order_excl$Order, na.rm = TRUE), 0)
category_totals_excl2 <- category_totals_excl2 %>%
  mutate(Order = ifelse(is.na(Order), max_order_excl + 1, Order)) %>%
  filter(!is.na(ROH_Size_Category))  # Remove any NA categories

p4_excl <- ggplot(category_totals_excl2, aes(x = reorder(Individual, -Order), y = Category_KB, fill = ROH_Size_Category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Population, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cat_colors, drop = FALSE, na.value = NA) +
  labs(title = "Individual ROH profiles (stacked by size category)",
       x = "Individual", y = "ROH length (KB)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

out_p4_excl <- file.path(outdir, "custom_individual_roh_profiles_excludeEG35.png")
ggsave(filename = out_p4_excl,
       plot = p4_excl, width = 18, height = 8, dpi = dpi_val, bg = plot_bg)
cat("Wrote:", out_p4_excl, "\n")

# Prepare ordering for remapped include version (use all_roh_remapped but keep all populations)
individual_order_remap <- individual_totals_remapped %>%
  arrange(Population, desc(Total_ROH_KB)) %>%
  group_by(Population) %>%
  mutate(Order = row_number()) %>%
  ungroup() %>%
  select(Population, Individual, Total_ROH_KB, Order)

category_totals_remap2 <- category_totals_remapped %>%
  left_join(individual_order_remap, by = c("Population", "Individual")) %>%
  arrange(Population, desc(Total_ROH_KB))

max_order_remap <- ifelse(nrow(individual_order_remap) > 0, max(individual_order_remap$Order, na.rm = TRUE), 0)
category_totals_remap2 <- category_totals_remap2 %>%
  mutate(Order = ifelse(is.na(Order), max_order_remap + 1, Order)) %>%
  filter(!is.na(ROH_Size_Category))  # Remove any NA categories

p4_remap <- ggplot(category_totals_remap2, aes(x = reorder(Individual, -Order), y = Category_KB, fill = ROH_Size_Category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~Population, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cat_colors, drop = FALSE, na.value = NA) +
  labs(title = "Individual ROH profiles (stacked)",
       x = "Individual", y = "ROH length (KB)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

out_p4_remap <- file.path(outdir, "custom_individual_roh_profiles_withEG35_remapped_BKN_UV_o_01.png")
ggsave(filename = out_p4_remap,
       plot = p4_remap, width = 18, height = 8, dpi = dpi_val, bg = plot_bg)
cat("Wrote:", out_p4_remap, "\n")

# Also save SVG versions for both
ggsave(filename = sub("\\.png$", ".svg", out_p4_excl), plot = p4_excl, width = 18, height = 8)
ggsave(filename = sub("\\.png$", ".svg", out_p4_remap), plot = p4_remap, width = 18, height = 8)


# PLOT 5: Ancient ROH (>=5000 KB) - two versions
-
ancient_excl <- all_roh_filtered %>% filter(KB >= 5000)
if (nrow(ancient_excl) > 0) {
  p5_excl <- ggplot(ancient_excl, aes(x = Individual, y = KB, fill = Population)) +
    geom_col() +
    facet_grid(~Population, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pop_color_map_filtered) +
    labs(title = "Very Recent ROH segments (>= 5000 kb)", x = "Individual", y = "ROH length (KB)") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold"))
  out_p5_excl <- file.path(outdir, "custom_ancient_roh_segments_excludeEG35.png")
  ggsave(filename = out_p5_excl, plot = p5_excl, width = 14, height = 6, dpi = dpi_val, bg = plot_bg)
  cat("Wrote:", out_p5_excl, "\n")
} else {
  message("No ancient ROH segments (>= 5000 kb) found in BKN & HYD (excludeEG35).")
}

# remapped include (EG_35 remapped to BKN_o)
ancient_remap <- all_roh_remapped %>% filter(KB >= 5000)
if (nrow(ancient_remap) > 0) {
  p5_remap <- ggplot(ancient_remap, aes(x = Individual, y = KB, fill = Population)) +
    geom_col() +
    facet_grid(~Population, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = remapped_color_map) +
    labs(title = "Very Recent ROH segments (>= 5000 kb)",
         x = "Individual", y = "ROH length (KB)") +
    theme_minimal(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          strip.text = element_text(face = "bold"))
  out_p5_remap <- file.path(outdir, "custom_ancient_roh_segments_withEG35_remapped_BKN_UV_o_01.png")
  ggsave(filename = out_p5_remap, plot = p5_remap, width = 14, height = 6, dpi = dpi_val, bg = plot_bg)
  cat("Wrote:", out_p5_remap, "\n")
} else {
  message("No Very Recent ROH segments (>= 5000 kb) found in remapped dataset.")
}


# PLOT 6: ROH length density (filtered only BKN & HYD)

p6 <- ggplot(all_roh_filtered, aes(x = KB, fill = Population)) +
  geom_density(alpha = 0.55) +
  scale_fill_manual(values = pop_color_map_filtered) +
  labs(title = "ROH segment length distribution (density)", x = "ROH length (KB)", y = "Density") +
  theme_minimal(base_size = base_size) +
  theme(plot.title = element_text(face = "bold")) +
  xlim(0, quantile(all_roh_filtered$KB, 0.99, na.rm = TRUE))

ggsave(filename = file.path(outdir, "custom_roh_length_distribution.png"),
       plot = p6, width = 12, height = 7, dpi = dpi_val, bg = plot_bg)


# Save CSV summaries
write.csv(individual_totals_filtered, file.path(outdir, "custom_roh_individual_totals_BKN_HYD.csv"), row.names = FALSE)
write.csv(population_category_summary_filtered, file.path(outdir, "custom_roh_category_summary_BKN_HYD.csv"), row.names = FALSE)
write.csv(category_totals_excl2, file.path(outdir, "custom_roh_category_totals_per_ind_BKN_HYD.csv"), row.names = FALSE)

# CSVs for remapped outputs (for special plots)
write.csv(individual_totals_remapped, file.path(outdir, "custom_roh_individual_totals_withEG35_remapped.csv"), row.names = FALSE)
write.csv(population_category_summary_remapped, file.path(outdir, "custom_roh_category_summary_withEG35_remapped.csv"), row.names = FALSE)
write.csv(category_totals_remap2, file.path(outdir, "custom_roh_category_totals_per_ind_withEG35_remapped.csv"), row.names = FALSE)

# Save mapping file that shows the remapping of EG_35 -> BKN_UV_o_01
remap_table <- all_roh_full %>%
  filter(Population == "EG_35") %>%
  distinct(Population, IID) %>%
  mutate(new_population = "BKN_o", new_individual = "BKN_UV_o_01")
if (nrow(remap_table) > 0) write.csv(remap_table, file.path(outdir, "EG35_to_BKN_UV_o_01_remap_table.csv"), row.names = FALSE)

# Text summary
sink(file.path(outdir, "custom_roh_analysis_report.txt"))
cat("CUSTOM ROH ANALYSIS REPORT\n")
cat("==========================\n\n")
cat("Input files:\n"); cat(paste(hom_files, collapse = "\n"), "\n\n")
cat("Total ROH segments (all files):", nrow(all_roh_full), "\n")
cat("Individuals analysed (all files):", length(unique(all_roh_full$IID)), "\n")
cat("Populations in all files:", paste(sort(unique(all_roh_full$Population)), collapse = ", "), "\n\n")
cat("Note: For most plots we excluded EG_35 and only used BKN & HYD.\n")
cat("Two special plots (individual profiles and ancient ROH) were written twice:\n")
cat("  - *_excludeEG35.png  : excludes EG_35 (only BKN & HYD present)\n")
cat("  - *_withEG35_remapped_BKN_UV_o_01.png : includes EG_35 rows remapped to Individual = BKN_UV_o_01 and Population = BKN_o\n\n")
cat("Color assigned for remapped sample (BKN_UV_o_01) is: ", col_EG35, " (same as previous EG_35 color)\n\n")
cat("Generated PNGs (1200 dpi) in directory:", outdir, "\n")
sink()

cat("\nDONE — plots and CSVs written to:", outdir, "\n")
