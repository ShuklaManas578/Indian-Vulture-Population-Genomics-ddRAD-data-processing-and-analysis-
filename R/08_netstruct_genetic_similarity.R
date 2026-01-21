#!/usr/bin/env Rscript

# 21 Genetic Similarity Networks Visualization (Rscript)
# netstruct_genetic_similarity.R

suppressPackageStartupMessages({
  library(optparse); library(parallel); library(igraph); library(ape)
})

option_list <- list(
  make_option(c("--pruned_dir"), type="character",
              default="Path/to/ld_prune_results",
              help="Directory containing pruned_relaxed and pruned_stringent files"),
  make_option(c("--dataset"), type="character", default="both",
              help="Which dataset to process: 'relaxed', 'stringent', or 'both'"),
  make_option(c("--thresholds"), type="character",
              default="0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.70,0.73,0.75,0.77,0.8,0.83,0.87,0.9,0.92,0.95,0.98",
              help="Comma-separated thresholds or 'quantiles:k"),
  make_option(c("--nperm"), type="integer", default=1000, help="Number of permutations"),
  make_option(c("--perm.method"), type="character", default="label", help="'label' or 'geno'"),
  make_option(c("--ncores"), type="integer", default = max(1, parallel::detectCores()-1),
              help="Parallel cores"),
  make_option(c("--dpi"), type="integer", default=1200, help="PNG DPI"),
  make_option(c("--svg"), type="logical", default=TRUE, help="Write SVG"),
  make_option(c("--pdf"), type="logical", default=FALSE, help="Write PDF"),
  make_option(c("--verbose"), type="logical", default=TRUE)
)
opt <- parse_args(OptionParser(option_list=option_list))

PRUNED_DIR <- normalizePath(opt$pruned_dir)
DATASET <- opt$dataset
NCORES <- opt$ncores
VERBOSE <- opt$verbose
dpi <- opt$dpi

logmsg <- function(...) if (VERBOSE) message("[", format(Sys.time(), "%F %T"), "] ", ...)


ORIG_EG35_NAME <- "EG_35"
RENAMED_EG35 <- "BKN_UV_o_01"
RENAMED_EG35_POP <- "BKN_o"  

COL_EG35 <- "#B5A642"

# Custom Color Scheme
get_community_colors <- function(n_communities) {
  
  base_colors <- c(
    "#DE3163",  
    "#008B8B",  
    "#B5A642",  
    "#C9A0DC", 
    "#FDBCB4"   
  )
  
  if (n_communities <= length(base_colors)) {
    return(base_colors[1:n_communities])
  } else {
    # Extend with additional distinct colors if needed
    additional_colors <- c("#FF7F00", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628")
    return(c(base_colors, additional_colors)[1:n_communities])
  }
}

# Enhanced Plotting Function with Larger Networks
create_professional_network_plot <- function(g, comm_df, metrics, label, threshold, mod, pval,
                                             filepath, width = 16, height = 12, dpi = 1200,
                                             format = "png") {
  
  if (igraph::vcount(g) == 0) {
    if (format == "png") png(filepath, width = width, height = height, units = "in", res = dpi)
    if (format == "svg") svg(filepath, width = width, height = height)
    if (format == "pdf") pdf(filepath, width = width, height = height)
    
    plot.new()
    title(main = paste0(label, " - Threshold: ", round(threshold, 4)),
          cex.main = 1.5, font.main = 2)
    text(0.5, 0.6, "No edges at this threshold", cex = 1.2, col = "darkred")
    text(0.5, 0.5, "Network is completely disconnected", cex = 1.0, col = "darkgray")
    text(0.5, 0.4, paste("Similarity threshold too high:", round(threshold, 4)),
         cex = 0.9, col = "darkgray")
    
    dev.off()
    return()
  }
  

  if (format == "png") png(filepath, width = width, height = height, units = "in", res = dpi)
  if (format == "svg") svg(filepath, width = width, height = height)
  if (format == "pdf") pdf(filepath, width = width, height = height)
  

  layout_matrix <- matrix(c(1, 1, 1, 1, 2, 3), nrow = 3, ncol = 2, byrow = TRUE)
  layout(layout_matrix, heights = c(0.7, 0.15, 0.15))
  
  # MAIN NETWORK PLOT
  par(mar = c(1, 1, 3, 1), bg = "white")
  
  # Calculate layout - use Fruchterman-Reingold for better spacing
  L <- igraph::layout_with_fr(g, niter = 1000)
  
  # Prepare community colors
  communities <- na.omit(unique(comm_df$community))
  n_comms <- length(communities)
  community_colors <- get_community_colors(n_comms)
  
  # Assign colors to nodes
  node_colors <- rep("#B3B3B3", igraph::vcount(g))  
  node_border <- rep("black", igraph::vcount(g))
  
  if (n_comms > 0) {
    for (i in seq_along(communities)) {
      comm_nodes <- which(comm_df$community == communities[i])
            if (length(comm_nodes) > 0) {
        comm_samples <- comm_df$sample[comm_nodes]
        v_idx <- match(comm_samples, igraph::V(g)$name)
        v_idx <- v_idx[!is.na(v_idx)]
        if (length(v_idx) > 0) node_colors[v_idx] <- community_colors[i]
      }
    }
  }
  

  if (RENAMED_EG35 %in% igraph::V(g)$name) {
    vpos <- which(igraph::V(g)$name == RENAMED_EG35)
    if (length(vpos) == 1) {
      node_colors[vpos] <- COL_EG35
    } else if (length(vpos) > 1) {
      node_colors[vpos] <- rep(COL_EG35, length(vpos))
    }
  }
  
  # Calculate node degrees for sizing - INCREASED SIZE RANGE (approx 20% larger)
  degrees <- igraph::degree(g)
  if (max(degrees) > min(degrees)) {
    node_size <- 6 + (degrees - min(degrees)) / (max(degrees) - min(degrees)) * 8
  } else {
    node_size <- rep(10, igraph::vcount(g))
  }
  

  edge_weights <- igraph::E(g)$weight
  
  if (!is.null(edge_weights) && length(edge_weights) > 0) {
    min_weight <- min(edge_weights)
    max_weight <- max(edge_weights)
    
    if (max_weight > min_weight) {
      edge_widths <- 1 + (edge_weights - min_weight) / (max_weight - min_weight) * 9
    } else {
      edge_widths <- rep(4, length(edge_weights))
    }
    
    edge_alpha <- 0.3 + (edge_weights - min_weight) / (max_weight - min_weight) * 0.5
    edge_alpha[is.na(edge_alpha)] <- 0.4
    edge_colors <- sapply(edge_alpha, function(alpha) {
      adjustcolor("#2F4F4F", alpha.f = alpha)
    })
    
  } else {
    edge_widths <- rep(3, igraph::ecount(g))
    edge_colors <- adjustcolor("#2F4F4F", alpha.f = 0.4)
  }

  plot(g, layout = L,
       vertex.color = node_colors,
       vertex.frame.color = node_border,
       vertex.size = node_size,
       vertex.label = igraph::V(g)$name,
       vertex.label.cex = 0.9,
       vertex.label.color = "black",
       vertex.label.family = "sans",
       vertex.label.font = 2,
       edge.color = edge_colors,
       edge.width = edge_widths,
       edge.curved = 0.15,
       main = paste0(label, " Network\nThreshold: ", round(threshold, 4)),
       cex.main = 1.3,
       frame = FALSE)
  
  # Add network statistics to plot
  legend("topleft",
         legend = c(
           paste("Samples:", igraph::vcount(g)),
           paste("Edges:", igraph::ecount(g)),
           paste("Communities:", n_comms),
           paste("Modularity:", round(mod, 4)),
           paste("Mean Degree:", round(mean(degrees), 2))
         ),
         bty = "o",
         box.col = "white",
         bg = adjustcolor("white", alpha.f = 0.9),
         cex = 1.0,
         text.col = "darkblue")
  
  # COMMUNITY LEGEND (Top Right)
  par(mar = c(2, 2, 2, 2))
  plot.new()
  if (n_comms > 0) {
    legend("center",
           legend = paste("Community", communities),
           pch = 19,
           col = community_colors[1:n_comms],
           pt.cex = 1.8,
           bty = "n",
           ncol = 1,
           cex = 1.1,
           title = "Communities",
           title.adj = 0.2,
           title.col = "darkred",
           title.cex = 1.1)
  } else {
    text(0.5, 0.5, "No communities", cex = 1.2, col = "darkgray")
  }
  
  # DEGREE LEGEND (Bottom Right) 
  par(mar = c(2, 2, 2, 2))
  plot.new()
  
  if (length(unique(degrees)) > 1) {
    degree_vals <- c(min(degrees), median(degrees), max(degrees))
    degree_sizes <- 6 + (degree_vals - min(degrees)) / (max(degrees) - min(degrees)) * 8
    
    legend("center",
           legend = paste("Degree:", degree_vals),
           pch = 19,
           col = "black",
           pt.cex = degree_sizes / 10,
           bty = "n",
           ncol = 1,
           cex = 1.1,
           title = "Node Degree",
           title.adj = 0.2,
           title.col = "darkred",
           title.cex = 1.1)
  } else {
    text(0.5, 0.5, paste("Degree:", degrees[1]), cex = 1.2, col = "darkgray")
  }
  
  dev.off()
}

# Core Functions 
is_plink_prefix <- function(p) {
  file.exists(paste0(p,".bed")) && file.exists(paste0(p,".bim")) && file.exists(paste0(p,".fam"))
}

clean_sample <- function(s) gsub("\\s+","", s)

get_geno_matrix <- function(plink_prefix, tmpdir = NULL) {
  if (is.null(tmpdir)) {
    tmpdir <- tempdir()
  } else {
    if (!dir.exists(tmpdir)) {
      dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  logmsg("Using temporary directory: ", tmpdir)
  
  if (!requireNamespace("SNPRelate", quietly=TRUE)) {
    stop("Install 'SNPRelate' for PLINK reading")
  }
  
  library(SNPRelate); library(gdsfmt)
  
  bed <- paste0(plink_prefix,".bed")
  bim <- paste0(plink_prefix,".bim")
  fam <- paste0(plink_prefix,".fam")
  
  gdsfile <- file.path(tmpdir, paste0("tmp_", basename(plink_prefix), ".gds"))
  if (!dir.exists(dirname(gdsfile))) dir.create(dirname(gdsfile), recursive = TRUE, showWarnings = FALSE)
  
  logmsg("SNPRelate BED2GDS -> ", gdsfile)
  
  SNPRelate::snpgdsBED2GDS(bed.fn=bed, bim.fn=bim, fam.fn=fam,
                           out.gdsfn=gdsfile, cvt.chr="char", verbose=FALSE)
  
  genofile <- SNPRelate::snpgdsOpen(gdsfile)
  on.exit({
    try(SNPRelate::snpgdsClose(genofile), silent=TRUE)
    if (file.exists(gdsfile)) unlink(gdsfile)
  }, add = TRUE)
  
  G <- SNPRelate::snpgdsGetGeno(genofile, with.id = FALSE)
  snp_ids <- read.gdsn(index.gdsn(genofile,"snp.id"))
  sample_ids <- read.gdsn(index.gdsn(genofile,"sample.id"))
  
  geno <- if (nrow(G) == length(snp_ids) && ncol(G) == length(sample_ids)) t(G) else G
  rownames(geno) <- sample_ids
  colnames(geno) <- snp_ids
  geno[geno < 0] <- NA_real_
  
  geno <- matrix(as.numeric(geno),
                 nrow=nrow(geno), ncol=ncol(geno),
                 dimnames=list(rownames(geno), colnames(geno)))
  
  logmsg("Genotype matrix: ", nrow(geno), " samples x ", ncol(geno), " SNPs")
  return(list(samples = rownames(geno), geno = geno))
}

pairwise_distance_matrix <- function(geno, ncores=1) {
  geno <- as.matrix(geno)
  n <- nrow(geno)
  samp <- rownames(geno)
  
  calc_row <- function(i) {
    gi <- geno[i,]
    diffs <- abs(sweep(geno, 2, gi, FUN="-"))/2
    rs <- rowMeans(diffs, na.rm=TRUE)
    rs[is.nan(rs)] <- NA_real_
    rs
  }
  
  if (ncores > 1 && .Platform$OS.type != "windows") {
    rows_list <- parallel::mclapply(seq_len(n), calc_row, mc.cores = ncores)
  } else {
    rows_list <- lapply(seq_len(n), calc_row)
  }
  
  D <- do.call(rbind, rows_list)
  rownames(D) <- samp
  colnames(D) <- samp
  diag(D) <- 0
  D
}

build_graph_at_threshold <- function(Smat, t) {
  A <- Smat
  A[A < t] <- 0
  diag(A) <- 0
  g <- igraph::graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE, diag=FALSE)
  igraph::simplify(g)
}

analyze_graph <- function(g) {
  if (igraph::ecount(g) == 0) {
    return(list(nc = 0, mod = NA_real_, membership = rep(NA_integer_, igraph::vcount(g))))
  }
  cm <- igraph::cluster_fast_greedy(g, weights = igraph::E(g)$weight)
  list(nc = length(unique(igraph::membership(cm))),
       mod = igraph::modularity(cm),
       membership = igraph::membership(cm))
}

# Main Processing Function
process_dataset <- function(plink_prefix, label, outdir_base, thresholds, nperm, ncores, dpi, write_svg, write_pdf) {
  logmsg("Processing dataset: ", label)
  
  OUTDIR <- file.path(outdir_base, label)
  
  if (!dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)) {
    stop("Failed to create output directory: ", OUTDIR)
  }
  logmsg("Output directory: ", OUTDIR)
  
  tmpdir <- file.path(tempdir(), paste0("netstruct_", label))
  if (!dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Failed to create temporary directory: ", tmpdir)
  }
  
  logmsg("Reading genotypes from: ", plink_prefix)
  Gobj <- get_geno_matrix(plink_prefix = plink_prefix, tmpdir = tmpdir)
  geno <- Gobj$geno
  
 
  rownames(geno) <- vapply(rownames(geno), clean_sample, FUN.VALUE = character(1))
  samples <- rownames(geno)
  
  match_eg35 <- which(samples == ORIG_EG35_NAME | grepl(paste0("\\b", ORIG_EG35_NAME, "\\b"), samples) |
                        grepl(paste0("^", ORIG_EG35_NAME, "_", ORIG_EG35_NAME, "$"), samples))
  if (length(match_eg35) > 0) {
    logmsg("Found sample(s) matching EG_35 pattern:", paste(samples[match_eg35], collapse = ", "))
    new_names <- samples
    new_names[match_eg35] <- RENAMED_EG35
    # ensure uniqueness after renaming
    if (any(duplicated(new_names))) new_names <- make.unique(new_names, sep = "_")
    rownames(geno) <- new_names
    logmsg("Renamed to:", paste(new_names[match_eg35], collapse = ", "))
  } else {
    rownames(geno) <- samples
  }
  
  logmsg("Computing pairwise distance matrix...")
  D <- pairwise_distance_matrix(geno, ncores = ncores)
  
  dist_file <- file.path(OUTDIR, paste0(label, "_pairwise_distances.csv"))
  write.csv(D, dist_file, row.names=TRUE)
  logmsg("Wrote distance matrix: ", dist_file)
  
  mean_dist <- round(rowMeans(D, na.rm=TRUE), 4)
  nn_val <- rep(NA_real_, length(mean_dist))
  nn_partner <- rep(NA_character_, length(mean_dist))
  
  for (i in seq_along(nn_val)) {
    di <- D[i,]
    di[i] <- NA
    if (!all(is.na(di))) {
      j <- which.min(di)
      nn_val[i] <- di[j]
      nn_partner[i] <- names(di)[j]
    }
  }
  
  metrics <- data.frame(
    sample = names(mean_dist),
    mean_dist = mean_dist,
    nn = round(nn_val,4),
    nn_partner = nn_partner,
    stringsAsFactors = FALSE
  )
  
  metrics_file <- file.path(OUTDIR, paste0(label, "_sample_metrics.csv"))
  write.csv(metrics, metrics_file, row.names=FALSE)
  logmsg("Wrote sample metrics: ", metrics_file)
  
  # Similarity matrix and thresholds
  S <- 1 - D
  S[S < 0] <- 0
  
  # Parse thresholds
  if (grepl("^quantiles:", thresholds)) {
    k <- as.numeric(sub("^quantiles:", "", thresholds))
    vals <- S[upper.tri(S)]
    probs <- seq(0.99, 0.99 - 0.01*(k-1), length.out = k)
    thresholds_vec <- as.numeric(quantile(vals, probs = probs, na.rm=TRUE))
  } else {
    thresholds_vec <- as.numeric(unlist(strsplit(thresholds, ",")))
  }
  
  logmsg("Processing ", length(thresholds_vec), " thresholds")
  logmsg("Threshold range: ", round(min(thresholds_vec),4), " to ", round(max(thresholds_vec),4))
  
  # Process each threshold
  for (t in thresholds_vec) {
    logmsg("Processing threshold ", round(t,4))
    g <- build_graph_at_threshold(S, t)
    an <- analyze_graph(g)
    obs_mod <- an$mod
    obs_nc <- an$nc
    
    logmsg("  Nodes: ", igraph::vcount(g),
           " Edges: ", igraph::ecount(g),
           " Communities: ", obs_nc,
           " Modularity: ", ifelse(is.na(obs_mod),"NA", round(obs_mod,4)))
    
    # Community assignment
    if (igraph::vcount(g) == 0) {
      comm_df <- data.frame(sample = rownames(D), community = NA_integer_, stringsAsFactors = FALSE)
    } else {
      memb <- an$membership
      comm_df <- data.frame(sample = igraph::V(g)$name, community = NA_integer_, stringsAsFactors = FALSE)
      if (!is.null(memb) && length(memb) > 0) {
        comm_df$community[match(names(memb), comm_df$sample)] <- as.integer(memb)
      }
    }
    
    # Create network plots
    base_name <- paste0(label, "_network_t", round(t,3))
    png_file <- file.path(OUTDIR, paste0(base_name, ".png"))
    svg_file <- file.path(OUTDIR, paste0(base_name, ".svg"))
    pdf_file <- file.path(OUTDIR, paste0(base_name, ".pdf"))
    
    # Create plot in multiple formats
    create_professional_network_plot(g, comm_df, metrics, label, t, obs_mod, NA,
                                     png_file, dpi = dpi, format = "png")
    logmsg("  Created PNG: ", basename(png_file))
    
    if (write_svg) {
      create_professional_network_plot(g, comm_df, metrics, label, t, obs_mod, NA,
                                       svg_file, format = "svg")
      logmsg("  Created SVG: ", basename(svg_file))
    }
    
    if (write_pdf) {
      create_professional_network_plot(g, comm_df, metrics, label, t, obs_mod, NA,
                                       pdf_file, format = "pdf")
      logmsg("  Created PDF: ", basename(pdf_file))
    }
    

    comm_file <- file.path(OUTDIR, paste0(label, "_communities_t", format(t, scientific=FALSE), ".csv"))
    write.csv(comm_df, comm_file, row.names=FALSE)
    
    summary_file <- file.path(OUTDIR, paste0(label, "_summary_t", format(t, scientific=FALSE), ".csv"))
    write.csv(data.frame(
      threshold = t,
      nodes = igraph::vcount(g),
      edges = igraph::ecount(g),
      ncomms = an$nc,
      modularity = an$mod,
      p_value = NA
    ), summary_file, row.names=FALSE)
  }
  
  # Clean up temporary directory
  unlink(tmpdir, recursive = TRUE)
  
  logmsg("Completed processing for: ", label)
  return(OUTDIR)
}


logmsg("Starting professional network analysis for Vultures datasets")
logmsg("Input directory: ", PRUNED_DIR)

# Define datasets
datasets <- list()
if (DATASET %in% c("relaxed", "both")) {
  datasets$relaxed <- list(
    prefix = file.path(PRUNED_DIR, "pruned_relaxed"),
    label = "pruned_relaxed"
  )
}
if (DATASET %in% c("stringent", "both")) {
  datasets$stringent <- list(
    prefix = file.path(PRUNED_DIR, "pruned_stringent"),
    label = "pruned_stringent"
  )
}

# Create main output directory
main_outdir <- file.path(PRUNED_DIR, "network_analysis_results")
if (!dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)) {
  stop("Failed to create main output directory: ", main_outdir)
}

logmsg("Main output directory: ", main_outdir)

# Process each dataset
output_dirs <- list()
for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  if (is_plink_prefix(ds$prefix)) {
    logmsg("=== Processing ", ds$label, " ===")
    tryCatch({
      output_dir <- process_dataset(
        plink_prefix = ds$prefix,
        label = ds$label,
        outdir_base = main_outdir,
        thresholds = opt$thresholds,
        nperm = opt$nperm,
        ncores = NCORES,
        dpi = opt$dpi,
        write_svg = opt$svg,
        write_pdf = opt$pdf
      )
      output_dirs[[ds_name]] <- output_dir
    }, error = function(e) {
      logmsg("ERROR processing ", ds$label, ": ", e$message)
    })
  } else {
    logmsg("Skipping ", ds$label, " - PLINK files not found at: ", ds$prefix)
  }
}

logmsg('\n ANALYSIS COMPLETE!')
logmsg('Output directory: ', normalizePath(main_outdir))
for (ds_name in names(output_dirs)) {
  logmsg(' ', ds_name, ': ', output_dirs[[ds_name]])
}

if (length(output_dirs) == 0) {
  logmsg(' No datasets were successfully processed. Check errors above.')
} else {
  logmsg('\n network plots created with:')
  logmsg('   • LARGER FONT SIZES for statistics and legends')
}
