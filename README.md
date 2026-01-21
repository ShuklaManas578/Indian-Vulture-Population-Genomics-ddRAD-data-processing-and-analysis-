# Indian-Vulture-Population-Genomics-ddRAD-data-processing-and-analysis-
Analysis code for ddRAD sequencing of Indian vultures (15 samples). This repository contains ordered bash pipeline for preprocessing, joint genotyping and population genetic analyses, R scripts for population-genetic summaries and plotting, Python helpers for allele polarization/dadi, and an environment.yml for reproducibility


---

## Overview of scripts

> NOTE: filenames are shown as in this repository (numeric prefixes to indicate recommended execution order).

### `scripts/bash/` — pipeline & QC (batch / cluster scripts)
These scripts implement the main preprocessing, mapping, variant calling, population-level VCF preparation and population genetic meteric estimation steps.

- `01_merge_and_trim_fastp.sh`
  - Purpose: merge raw reads (biological duplicates), run `fastp` for adapter trimming and quality filtering.
  - Inputs: raw FASTQ pairs. Outputs: cleaned FASTQ, per-sample QC reports.
  - Popgen note: good trimming and adapter removal reduce false SNPs caused by low-quality bases.

- `02_map_and_stats_bwa.sh`
  - Purpose: map reads to reference genome assmebly of Golden eagle with `bwa mem`, sort/index BAM with `samtools`, collect mapping stats.
  - Outputs: `.bam` files, `.bai`, mapping summary (mapping rate, coverage).
  - Popgen note: consistent mapping parameters are essential for comparable genotype calls.

- `03_regroup_and_QC_vcf.sh`
  - Purpose: regroup VCFs if the population clustering is to be tweaked, run basic VCF QC (vcftools/bcftools) to remove low-depth/low-quality sites.
  - Outputs: filtered VCFs ready for joint genotyping or further filtering.

- `04_run_GATK_parallel_genotyping.sh`
  - Purpose: run GATK joint-calling / GenotypeGVCFs or similar in parallel across the genome and apply genotype filters as needed.
  - Outputs: high-quality VCF for population analyses.
  - Popgen note: joint genotyping prevents genotyping differences caused by calling samples separately.

- `05_Gatk_filter_varification.sh`
  - Purpose: apply GATK hard/variant filtration and validation checks (e.g., VQSR, manual thresholds).

- `06_ld_prune_and_lddecay.sh`
  - Purpose: perform LD pruning and produce LD decay metrics/plots (e.g., using `plink` or `vcftools`).
  - Outputs: pruned sets for PCA,ADMIXTURE etc. and LD decay summaries.

- `07_fst_from_pruned_bfiles.sh`
  - Purpose: convert VCF → PLINK, LD-prune, compute pairwise FST with PLINK/VCFtools.
  - Popgen note: pruning before FST reduces linkage bias.

- `08_sliding_window_tajimasD.sh`
  - Purpose: compute Tajima’s D in sliding windows across the genome (vcftools or custom scripts).
  - Popgen note: windows with extreme Tajima's D may indicate selection or demographic change.

- `09_compute_overlapping_tajimas_vcftools.sh`
  - Purpose: compute per-site / per-window diversity statistics (π, segregating sites) in overlapping windows.

- `10_Nucleotide_diversity_pi.sh`
  - Purpose: compute nucleotide diversity (π) across populations or windows.
  - Per-site pi_windows_summary_pop.csv with standard deviation also reported.

- `11_run_roh_all_pops.sh`
  - Purpose: call runs of homozygosity (RoH) per sample (plink/roh tools).

- `12_run_admixture.sh`
  - Purpose: run ADMIXTURE across a range of K values and collect CV error for model selection.
  - Outputs: P and Q files per-replicate for each K run and plink files with admixture.stdout and admixture.sterr also generated.

- `13_Constructing_folded_SFS.sh`
  - Purpose: Prepares a 2-D joint folded spectra BKN and HYD populations from best projection using `angsd` and `realSFS`.


---

### `scripts/R/` — plotting and population-stat visualizations (R scripts)
These scripts generate the figures used in the manuscript (PCA, FST heatmaps, heterozygosity plots, RoH plots, diversity plots, mutation load plots, network/structure plots).

- `01_pca_from_pruned_sets.R` — PCA plot generation from pruned PLINK datasets.
- `02_make_all_fst_plots.R` — FST heatmaps, pairwise barplots or matrices.
- `03_make_admixture_plots_with_cvSD.R` — Admixture barplots with cross-validation error and standard deviation.
- `04_heterozygosity_plots.R` — per-sample heterozygosity calculation and violin/boxplots.
- `05_plot_pi_v2.R` — nucleotide diversity (π) plotting across sliding windows tuned to optimize SNP density.
- `06_plot_tajima.R` — Overlapping Tajima’s D summaries & genomic distributions.
- `07_mutload_plot_final.R` — mutation load and derived allele frequency spectrum plots.
- `08_netstruct_genetic_similarity.R` — network/graph-based visualization of genetic similarity using threshold clustering.
- `09_RoH-based_plots.R` — RoH length distributions and stacked plots and per-sample RoH summaries.

**R packages commonly used here**: `ggplot2`, `dplyr`, `tidyr`, `cowplot/patchwork`, `adegenet`, `pegas`, `hierfstat`, `ggtree`/`igraph` (where networks used), `cowplot`. (See `environment.yml` below for starter list.)

---

### `scripts/python/` — helpers and demographic inference
- `01_make_outgroup_vcf_fixed2.py` — helper to construct outgroup VCF or filter variants for outgroup-based polarizations.
- `02_polarize_and_mutload_fixed.py` — polarize alleles to ancestral state and compute derived allele counts / mutation load metrics.
- `03_dadi_pipeline.py` — dadi inference script to build SFS and run demographic model fitting.

**Python packages commonly used here**: `numpy`, `pandas`, `dadi`, `scikit-allel` or `cyvcf2`, `matplotlib`, `seaborn`.

---

## How to use this repo (recommended quick steps)
1. **Inspect** every script and change absolute paths or dataset-specific variables to your local layout before running anything.
2. Create the conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate vulture_analysis
