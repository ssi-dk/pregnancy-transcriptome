#!/usr/bin/env Rscript

#************************************
# WEIGHTED GRAPH CORRELATION NETWORK ANALYSIS
#
# DESCRIPTION: This scripts performs WGCNA on the genes found to be
# significantly differently expressed between the pregnancy weeks and pp,
# with the threshold selection included.
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# OBS: Needed: module load openssl/1.1.1s

# Imports and initialization
#-------------------------------------------------------------------------------
# Librariesx
library(tidyverse)
library(WGCNA)
library(DESeq2)
library(rtracklayer)
library(gridExtra)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))
source(file.path(path_project, "scripts", "functions_plotting.R"))

# Load list of significant genes
sig_genes_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "replicated_significant_genes_and_p_values.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Load model coeffecient for module plots
model_coef_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Load (wrangled/cleaned) count files
counts_b2 <- read.table(file = file.path(path_data, "03_gene_filtration", "b2_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Design (keys with meta data)
design_data <- read.table(file = file.path(path_data, "04_normalization", "design_finished.tsv"), as.is = TRUE, header = TRUE)

# Gencode v28 primary assembly GTF for gene symbol
gencode_gtf <- rtracklayer::import("PROJECT_PATH/reference/gencode.v28.primary_assembly.annotation.gtf")

# Model pvalues
pvals_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "significant_genes_table.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Data wrangling
#-------------------------------------------------------------------------------
# Ensure correct order of design/meta data
design_b2 <- design_data %>%
  filter(batch == "2") %>%
  arrange(match(sample, colnames(counts_b2)))

rownames(design_b2) <- pull(design_b2, sample)

## Extract list of significant genes
significant_genes <- pull(sig_genes_df, gene)

# Obtain ENSG and gene symbol/id from gtf file
ENSG_gene_symbol <- gencode_gtf %>% 
  as.data.frame() %>%
  select(gene_id, gene_name) %>% 
  unique() %>%
  mutate(gene = gene_id) %>%
  select(-gene_id) 
  
# Normalization: Variance stabilization transformation
#-------------------------------------------------------------------------------
# NB: We perform variance stabilization on all genes before extracting the genes of interest
# Create dds
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_b2),
                              colData = design_b2,
                              design = ~ 1) # We are not specifying a model

# Perform variance stabilization
dds_norm <- vst(dds)

# Extract normalized counts 
counts_vst_b2 <- assay(dds_norm) 

# Filter to only contain significant genes
counts_vst_sig_b2 <- as.matrix(counts_vst_b2[significant_genes,])


# Make PCA after VST normalization colored by PP/pregnant to check bio variability 
pca_plot_after_VST_norm <- make_PCA_plot(rownames_to_column(as.data.frame(counts_vst_b2[significant_genes,]), var = 'gene'), design_b2, 'postdelivery')$plot


# Normalization: Linear model residual extraction
#-------------------------------------------------------------------------------
# Run model for each gene and extract residual as normalized count
counts_norm_sig_b2 <- sapply(X = 1:nrow(counts_vst_sig_b2), FUN = function(gene_number) extract_gene_count_residuals(gene_number, counts_vst_sig_b2, design_b2))

# Rename
rownames(counts_norm_sig_b2) <- colnames(counts_vst_sig_b2)
colnames(counts_norm_sig_b2) <- rownames(counts_vst_sig_b2)

# Make PCA after final normalization
pca_plot_after_norm <- make_PCA_plot(rownames_to_column(as.data.frame(t(counts_norm_sig_b2)), var = "gene"), design_b2, 'postdelivery')$plot


# Network construction: Choosing soft power threshold
#-------------------------------------------------------------------------------
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))
network_type = 'signed'

# Call the network topology analysis function
soft <- pickSoftThreshold(counts_norm_sig_b2,
                          powerVector = power,
                          networkType = network_type,
                          verbose = 5)

soft_data <- soft$fitIndices

# Visualisation
a1 <- ggplot(soft_data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(soft_data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


soft_plot <- grid.arrange(a1, a2, nrow = 2)

# Prepare source data
soft_source_data <- soft_data %>% 
  select(Power, SFT.R.sq, mean.k.)



# Network construction: Module identification (block wise)
#-------------------------------------------------------------------------------

# Convert matrix to numerix
counts_norm_sig_b2[] <- sapply(counts_norm_sig_b2, as.numeric)

# Define selected soft power threshold
soft_power <- 16

# Temporarily redefine cor function for WGCNA (will cause error otherwise)
tmp_cor <- cor
cor <- WGCNA::cor

# Memory estimate w.r.t blocksize
#block_wise_network <- blockwiseModules(counts_norm_sig_b2,
#                          maxBlockSize = 5000, # Maximum number of genes in each block
#                          TOMType = 'signed',
#                          power = soft_power,
#                          mergeCutHeight = 0.25, # Threshold to merge similar modules
#                          numericLabels = FALSE, # For color labels instead of numbers
#                          randomSeed = 1234,
#                          verbose = 3)
# Reset cor
#cor <- tmp_cor



# Try different thresholds

for (cut_threshold in c(0.14, 0.16, 0.20)){
  print(cut_threshold)

  block_wise_network <- blockwiseModules(counts_norm_sig_b2,
                                         maxBlockSize = 5000, # Maximum number of genes in each block
                                         TOMType = 'signed',
                                         networkType = network_type,
                                         power = soft_power,
                                         mergeCutHeight = cut_threshold, # Threshold to merge similar modules
                                         numericLabels = FALSE, # For color labels instead of numbers
                                         randomSeed = 1234,
                                         verbose = 3)
  
  print(table(block_wise_network$colors))
  
  png(paste0(path_data, '/06_WGCNA/threshold_plots/WGCNA_dendrogram_merged_and_unmerged_threshold_', as.character(cut_threshold), '.png'), width = 900, height = 1200, units = 'px', res=120)
  print(plotDendroAndColors(block_wise_network$dendrograms[[1]], cbind(block_wise_network$unmergedColors, block_wise_network$colors),
                            c('unmerged', 'merged'),
                            dendroLabels = FALSE,
                            addGuide = TRUE,
                            hang= 0.03,
                            guideHang = 0.05))
  dev.off()
  
  cluster_sizes <- table(block_wise_network$colors) %>%
    as.data.frame() 

  col <- as.character(cluster_sizes$Var1)
  names(col) <- as.character(cluster_sizes$Var1)

  cluster_sizes_plot <- cluster_sizes %>%
    ggplot(aes(x = Var1, y = Freq, fill = Var1)) + 
    geom_bar(stat='identity') +
    scale_fill_identity() +
    theme(axis.text.y=element_text(size=14)) + 
    xlab('Module') + ylab('Number of genes') +
    labs(title = paste('mergeCutHeight =', as.character(cut_threshold)))

  png(paste0(path_data, '/06_WGCNA/threshold_plots/WGCNA_cluster_sizes_threshold_', as.character(cut_threshold), '.png'), width = 900, height = 500, units = 'px', res=120)
  print(cluster_sizes_plot)
  dev.off()

}


# With selected threshold

cut_threshold <- 0.14

print(paste("Selected merge cut threshold:", cut_threshold))

block_wise_network <- blockwiseModules(counts_norm_sig_b2,
                                       maxBlockSize = 5000, # Maximum number of genes in each block
                                       TOMType = 'signed',
                                       networkType = network_type,
                                       power = soft_power,
                                       mergeCutHeight = cut_threshold, # Threshold to merge similar modules
                                       numericLabels = FALSE, # For color labels instead of numbers
                                       randomSeed = 1234,
                                       verbose = 3)

png(paste0(path_data, '/06_WGCNA/WGCNA_dendrogram_merged_and_unmerged_threshold_', as.character(cut_threshold), '.png'), width = 900, height = 1200, units = 'px', res=120)
print(plotDendroAndColors(block_wise_network$dendrograms[[1]], cbind(block_wise_network$unmergedColors, block_wise_network$colors),
                          c('unmerged', 'merged'),
                          dendroLabels = FALSE,
                          addGuide = TRUE,
                          hang= 0.03,
                          guideHang = 0.05))
dev.off()

##### Extract source data #####
# Dendrogram gene position and bottom plot bars
dendro_position_bar <- cbind(block_wise_network$unmergedColors, block_wise_network$colors)

# Dendrogram structure
hc <- as.hclust(block_wise_network$dendrograms[[1]])

# Create table describing each merge
dendrogram_merge_data <- tibble(
  merge_step = seq_len(nrow(hc$merge)),
  left_branch = hc$merge[, 1],
  right_branch = hc$merge[, 2],
  height = hc$height
)

cor <- tmp_cor
print("finished WGCNA")

# Plot module sizes
#-------------------------------------------------------------------------------
cluster_sizes <- table(block_wise_network$colors) %>% # Contains which color/cluster each gene is assigned to
  as.data.frame() 

col <- as.character(cluster_sizes$Var1)
names(col) <- as.character(cluster_sizes$Var1)

cluster_sizes_plot <- cluster_sizes %>%
  ggplot(aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat='identity') +
  scale_fill_identity() +
  theme(axis.text.y=element_text(size=14)) + 
  xlab('Module') + ylab('Number of genes')



# Extract genes and their assigned module color
#-------------------------------------------------------------------------------
module_genes <- as.data.frame(block_wise_network$colors) 
colnames(module_genes) <- c("color") # Overwrite $ format in colname

module_genes <- module_genes %>%
  rownames_to_column(var = "gene")

# Extract hub genes for modules
#-------------------------------------------------------------------------------
# Extract hub genes
hub_genes <- WGCNA::chooseTopHubInEachModule(datExpr = counts_norm_sig_b2, colorh = block_wise_network$colors, type = network_type, power = soft_power) # Note old documentation says to use "labels2color", but we already have color assignments in network

print("hub gene formation")
hub_genes <- hub_genes %>%
  as.data.frame() %>%
  rownames_to_column(var = "hub_color") 

colnames(hub_genes) <- c("hub_color", "gene") 

hub_genes <- hub_genes %>%
  left_join(ENSG_gene_symbol, by = "gene")


# Extract most significant gene in each module
#-------------------------------------------------------------------------------
print("sig gene formation")

sig_module_genes <- module_genes %>%
  full_join(pvals_df, by = "gene") %>%
  select(gene, color, pvalue_discovery, gene_name) %>%
  group_by(color) %>%
  arrange(pvalue_discovery) %>%
  filter(row_number() == 1) %>% # slice doesn't work 
  ungroup()

print(sig_module_genes)

# Plot fold change segmental plots for all genes in each module and highlight hub and most significant gene
#-------------------------------------------------------------------------------
print("plotting")
# Wrangle data for segmental plotting
model_estimates_modules <- module_genes %>%
  #rownames_to_column(var = "gene") %>%
  left_join(ENSG_gene_symbol, by = "gene") %>%
  left_join(model_coef_b2, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect)))
  )


for (module_color in pull(distinct(model_estimates_modules, color))){

  # Extract hub, significant gene and pvalue for this module
  hub_gene_color <- filter(hub_genes, hub_color == module_color) 
  hub_gene <- pull(hub_gene_color, gene)

  sig_module_genes_color <- filter(sig_module_genes, color == module_color)
  sig_gene <- pull(sig_module_genes_color, gene)

  pval <- pull(sig_module_genes_color, pvalue_discovery) %>% format(digits = 3) %>% as.character()

  # Define colors for plot
  hub_gene_str <- paste0("Hub gene\n(",  as.character(pull(hub_gene_color, gene_name)), ")")
  sig_gene_str <- paste0("Top gene\n(", as.character(pull(sig_module_genes_color, gene_name)), ",\np-value\n", pval, ")")
  #sig_gene_str <- paste0("Significant gene (", as.character(pull(sig_module_genes_color, gene_name)), ")")


  #hub_colors <- c(paste("Hub", "gene") = "grey88", 
  #                paste("Significant", "gene") = "black", 
  #                "Module gene" = as.character(module_color))
  
  # Define coloring in data
  if (module_color != "grey"){ # No hub gene for grey modules
    module_data <- model_estimates_modules %>%
    filter(color == module_color) %>%
    mutate(hub_status = case_when(gene == hub_gene ~ hub_gene_str,
                                  gene == sig_gene ~ sig_gene_str,
                                  TRUE ~ "Module gene")) %>%
    mutate(hub_status = factor(hub_status, levels = c("Module gene", hub_gene_str, sig_gene_str))) %>%
    arrange(hub_status) 
  } else {
    module_data <- model_estimates_modules %>%
    filter(color == module_color) %>%
    mutate(hub_status = case_when(gene == sig_gene ~ sig_gene_str,
                                  TRUE ~ "Module gene")) %>%
    mutate(hub_status = factor(hub_status, levels = c("Module gene", sig_gene_str))) %>%
    arrange(hub_status) 
  }

  print(module_color)
  print(hub_gene_str)
  print(sig_gene_str)


  # Plot log fold change compared to post partum for each module
  module_plot <- module_data %>%
    ggplot() +
    geom_segment(aes(x = week_number, xend = week_end, y = Estimate, yend = Estimate, color = hub_status, linewidth = hub_status)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    #scale_color_manual(values=hub_colors) +
    scale_color_manual(breaks = c("Module gene", sig_gene_str, hub_gene_str),
                       values = c(module_color, "ivory4", "ivory3")) +
    scale_linewidth_manual(breaks = c("Module gene", sig_gene_str, hub_gene_str),
                       values = c(1, 3, 3)) +
    theme_bw() +
    labs(title = paste0("Module ", module_color), x = "Pregnancy week", y = "Log(fold change) to PP level") +
    theme(plot.title = element_text(size = 20),
          axis.text = element_text(size = 17),
          axis.title = element_text(size = 18), 
          legend.position = "right", #"bottom",
          legend.title = element_blank(), 
          legend.text = element_text(size = 18),
          legend.spacing.y = unit(0.5, "cm")) +
          #legend.box.background = element_rect(fill = "transparent", color = "transparent"))
    guides(color = guide_legend(byrow = TRUE))
  
  # Print module plot
  print(module_plot)
  #png(paste0(path_, "06_WGCNA/LogFoldChange_module_", as.character(module_color), ".png"), width = 800, height = 600)
  #print(module_plot)
  #dev.off()  


}



# Output
#-------------------------------------------------------------------------------
# PCA plot after VST
png(file.path(path_results, '06_WGCNA", "PCA_after_VST_norm_postdelivery.png'))
print(pca_plot_after_VST_norm)
dev.off()

# PCA plot after VST + residual normalization
png(file.path(path_results, "06_WGCNA", "PCA_after_norm_postdelivery.png"))
print(pca_plot_after_norm)
dev.off()

# Soft threshold plot
ggsave(filename=file.path(path_data, "06_WGCNA", "threshold_plots", "soft_threshold_plot.png"), plot=soft_plot)
# Source data 
write.table(x=soft_source_data, file = file.path(path_results, "06_WGCNA", "plot_tables", "soft_threshold_plot_source_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


# WGCNA threshold plot
png(paste0(path_data, '/06_WGCNA/WGCNA_cluster_sizes_threshold_', as.character(cut_threshold), '.png'), width = 900, height = 500, units = 'px', res=120)
print(cluster_sizes_plot)
dev.off()
# Source data
write.table(x=cluster_sizes, file = file.path(path_results, "06_WGCNA", "plot_tables", "cluster_sizes_barplot_source_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

# Source data for dendrogram
write.table(x=dendro_position_bar, file = file.path(path_results, "06_WGCNA", "plot_tables", "dendro_gene_position.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)
write.table(x=dendrogram_merge_data, file = file.path(path_results, "06_WGCNA", "plot_tables", "dendro_merge_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "06_WGCNA", "session_info", "06a_session_info.txt"))
