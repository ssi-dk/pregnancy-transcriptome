#!/usr/bin/env Rscript

#************************************
# WEIGHTED GRAPH CORRELATION NETWORK ANALYSIS
#
# DESCRIPTION: This scripts performs WGCNA on the genes found to be
# significantly differently expressed between the pregnancy weeks and pp
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************


# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(WGCNA)
library(DESeq2)
library(rtracklayer)
library(ggpubr)
library(scales)
library(grid)


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
gencode_gtf <- rtracklayer::import("/home/projects/ssi_gen1/data/pregnancy/reference/gencode.v28.primary_assembly.annotation.gtf")

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
# Make PCA plot before normalization 
pca_before_any_norm <- make_PCA_plot(rownames_to_column(as.data.frame(counts_b2[significant_genes,]), var = 'gene'), design_b2, 'postdelivery')
pca_plot_before_any_norm <- pca_before_any_norm$plot
pca_data_before_any_norm <- pca_before_any_norm$data

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
pca_after_VST_norm <- make_PCA_plot(rownames_to_column(as.data.frame(counts_vst_b2[significant_genes,]), var = 'gene'), design_b2, 'postdelivery')
pca_plot_after_VST_norm <- pca_after_VST_norm$plot
pca_data_after_VST_norm <- pca_after_VST_norm$data

# Normalization: Linear model residual extraction
#-------------------------------------------------------------------------------
# Run model for each gene and extract residual as normalized count
counts_norm_sig_b2 <- sapply(X = 1:nrow(counts_vst_sig_b2), FUN = function(gene_number) extract_gene_count_residuals(gene_number, counts_vst_sig_b2, design_b2))

# Rename
rownames(counts_norm_sig_b2) <- colnames(counts_vst_sig_b2)
colnames(counts_norm_sig_b2) <- rownames(counts_vst_sig_b2)

# Make PCA after final normalization
pca_after_norm <- make_PCA_plot(rownames_to_column(as.data.frame(t(counts_norm_sig_b2)), var = "gene"), design_b2, 'postdelivery')
pca_plot_after_norm <- pca_after_norm$plot
pca_data_after_norm <- pca_after_norm$data


# Network construction: Module identification (block wise)
#-------------------------------------------------------------------------------
# Temporarily redefine cor function for WGCNA (will cause error otherwise)
tmp_cor <- cor
cor <- WGCNA::cor

# Convert matrix to numerix
counts_norm_sig_b2[] <- sapply(counts_norm_sig_b2, as.numeric)

# With selected threshold
soft_power <- 16
cut_threshold <- 0.14
network_type = 'signed'

print(paste("Selected merge cut threshold:", cut_threshold))

# Create network
block_wise_network <- blockwiseModules(counts_norm_sig_b2,
                                       maxBlockSize = 5000, # Maximum number of genes in each block
                                       TOMType = 'signed',
                                       networkType = network_type,
                                       power = soft_power,
                                       mergeCutHeight = cut_threshold, # Threshold to merge similar modules
                                       numericLabels = FALSE, # For color labels instead of numbers
                                       randomSeed = 1234,
                                       verbose = 3)

# Redefine cor 
cor <- tmp_cor
print("finished WGCNA")


# Extract genes and their assigned module color
#-------------------------------------------------------------------------------
module_genes_b1 <- as.data.frame(block_wise_network$colors) 
colnames(module_genes_b1) <- c("color") # Overwrite $ format in colname

module_genes_b1 <- module_genes_b1 %>%
  rownames_to_column(var = "gene")


# Plot fold change segmental plots for all genes in each module and highlight hub and most significant gene
#-------------------------------------------------------------------------------
# Wrangle data for segmental plotting
model_estimates_modules <- module_genes_b1 %>%
  left_join(ENSG_gene_symbol, by = "gene") %>%
  left_join(model_coef_b2, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect)))
  ) 


# Data wrangling
#-------------------------------------------------------------------------------
# Extract the pregnancy-associated genes
sig_module_genes <- module_genes_b1 %>%
  full_join(pvals_df, by = "gene") %>%
  select(gene, color, pvalue_discovery, gene_name) %>%
  group_by(color) %>%
  arrange(pvalue_discovery) %>%
  filter(row_number() == 1) %>% # slice doesn't work 
  ungroup()

# Prepare data for log fold plot
model_estimates_modules <- module_genes_b1 %>% 
  left_join(ENSG_gene_symbol, by = "gene") %>%
  left_join(model_coef_b2, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))))


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

# Prepare plot data
#-------------------------------------------------------------------------------
plot_data <- model_estimates_modules

# Log fold change plots
#-------------------------------------------------------------------------------
# List for plot objects
module_plot_list <- list()
module_plot_data_list <- list()

# Color lists
# module_colors <- pull(distinct(model_estimates_modules, color))
module_colors <- c("green", "yellow", "black", "pink", "red", "brown", "turquoise", "blue", "grey")
nice_module_color <- setNames(c("green3", "gold1", "grey5", "pink2", "indianred2", "tan4" , "mediumturquoise", "dodgerblue2", "grey"), module_colors)
highlight_color <- setNames(c("green4", "tan3", "blue3", "pink4", "darkred", "#3B1F0B", "turquoise4", "midnightblue", "azure4"), module_colors)



# Loop through modules to create log-fold plots for all
for (module_color in module_colors){
  # Extract hub, significant gene and pvalue for this module
  hub_gene_color <- filter(hub_genes, hub_color == module_color) 
  hub_gene <- pull(hub_gene_color, gene_name)
  
  sig_module_genes_color <- filter(sig_module_genes, color == module_color)
  sig_gene <- pull(sig_module_genes_color, gene_name)
  
  pval <- pull(sig_module_genes_color, pvalue_discovery) %>% format(digits = 3) %>% as.character()
  
  # Define colors for plot
  hub_gene_str <- paste0("Hub gene\n(",  as.character(pull(hub_gene_color, gene_name)), ")")
  sig_gene_str <- paste0("Top gene (", as.character(pull(sig_module_genes_color, gene_name)), ",\np = ", pval, ")")

  
  # Define coloring in data
  if (module_color != "grey"){ # No hub gene for grey modules
    module_data <- model_estimates_modules %>%
      filter(color == module_color) %>%
      mutate(hub_status = case_when(gene_name == hub_gene ~ hub_gene_str,
                                    gene_name == sig_gene ~ sig_gene_str,
                                    TRUE ~ "Module\n. gene")) %>%
      mutate(hub_status = factor(hub_status, levels = c("Module\n. gene", hub_gene_str, sig_gene_str))) %>%
      arrange(hub_status) 
  } else {
    module_data <- model_estimates_modules %>%
      filter(color == module_color) %>%
      mutate(hub_status = case_when(gene_name == sig_gene ~ sig_gene_str,
                                    TRUE ~ "Module\n. gene")) %>%
      mutate(hub_status = factor(hub_status, levels = c("Module\n. gene", sig_gene_str))) %>%
      arrange(hub_status) 
  }
  
  # Number of genes in module
  n_module <- module_data %>% 
    dplyr::select(gene) %>% 
    unique() %>% 
    nrow()
  
  # Plot log fold change compared to post partum for each module
  module_plot_list[[module_color]] <- module_data %>%
    ggplot() +
    geom_segment(aes(x = week_number, xend = week_end, y = Estimate, yend = Estimate, color = hub_status, linewidth = hub_status, alpha = hub_status)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    #scale_color_manual(values=hub_colors) +
    scale_color_manual(breaks = c("Module\n. gene", hub_gene_str, sig_gene_str),
                       values = c(nice_module_color[[module_color]], highlight_color[[module_color]], "ivory4")) +
    #values = c(module_color, "#D7BFBF","#A9B88D")) +
    scale_linewidth_manual(breaks = c("Module\n. gene", hub_gene_str, sig_gene_str),
                           values = c(1, 3, 2)) +
    scale_alpha_manual(breaks = c("Module\n. gene", hub_gene_str, sig_gene_str),
                       values = c(0.7, 1, 1)) +
    theme_bw() +
    labs(title = paste0("Module ", module_color, " (n = ", n_module, ")")) + #, x = "Pregnancy week", y = "Log(fold change) to PP level") +
    theme(plot.title = element_text(size = 25),
          axis.text = element_text(size = 18),
          axis.title = element_blank(),
          #axis.title = element_text(size = 19), 
          legend.position = "bottom",
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          legend.spacing.y = unit(0.5, "cm"),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.box.just = "center",                # centers legend contents
          legend.box.margin = margin(t = 5),
    ) +
    guides(color = guide_legend(byrow = TRUE))
  
  # Save source data
  module_plot_data_list[[module_color]] <- module_data
  
}

module_plot_all <- do.call(ggarrange, c(module_plot_list, ncol = 3, nrow = 3, align="hv"))#, common.legend = TRUE, legend="right"))

module_plot_all_label <-annotate_figure(module_plot_all , left = textGrob("Log(fold change) to PP level", rot = 90, vjust = 1, gp = gpar(cex = 2.3)),
                                        bottom = textGrob("Gestational age (weeks)", gp = gpar(cex = 2.3)))

# Gather all source data
module_plot_source_data <- bind_rows(module_plot_data_list)


# Output
#-------------------------------------------------------------------------------
# PCA before any normalization
# PCA plot before any normalization
png(file.path(path_results, "06_WGCNA", "PCA_before_any_normalization.png"))
print(pca_plot_before_any_norm)
dev.off()
# PCA data before any normalization
write.table(x = pca_data_before_any_norm, file = file.path(path_results, "06_WGCNA", "plot_tables", "PCA_before_any_normalization_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


# PCA plot after VST
png(file.path(path_results, "06_WGCNA", "PCA_after_VST_norm_postdelivery.png"))
print(pca_plot_after_VST_norm)
dev.off()
# PCA data before any normalization
write.table(x = pca_data_after_VST_norm, file = file.path(path_results, "06_WGCNA", "plot_tables", "PCA_after_VST_norm_postdelivery_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


# PCA plot after final normalization
png(file.path(path_results, "06_WGCNA", "PCA_after_norm_postdelivery.png"))
print(pca_plot_after_norm)
dev.off()
# PCA data final normalization
write.table(x = pca_data_after_norm, file = file.path(path_results, "06_WGCNA", "plot_tables", "PCA_after_norm_postdelivery_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


# Module plots
png(file.path(path_results, "06_WGCNA", "LogFoldChange_module_all.png"), width = 1600, height = 1600)
print(module_plot_all_label)
dev.off() 
# Module plot source data
write.table(x = module_plot_source_data, file = file.path(path_results, "06_WGCNA", "plot_tables", "module_plot_source_data.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)


# Write df with genes, their assigned modules and module estimates
write.table(x = model_estimates_modules, file = file.path(path_data, "06_WGCNA", "gene_module_model_estimates.tsv"), row.names = FALSE, sep = "\t", quote = FALSE)

# Plot data

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "06_WGCNA", "session_info", "06b_session_info.txt"))

