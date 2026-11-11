#!/usr/bin/env Rscript

#************************************
# MAKE COMPARISON WITH GENE LISTS
#
# This script loads tables from other papers to compare the genes they have
# found to our own gene list
# 
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization 
#-------------------------------------------------------------------------------
# Libraries 
library("officer")
library("tidyverse")
library("readxl")
library(docxtractr)
library(plyr)
library(ggVennDiagram)


# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))


# Read in files (own Wright, and Gomez-Lopez)
wright_doc <- read_docx(file.path(path_data, "100_data_exploration", "Table_1.DOCX"))
gomez_table <- read_excel(file.path(path_data, "100_data_exploration", "Data_Sheet.XLS"), sheet = 'SupplementaryFile1', skip = 1, col_names = TRUE)
sig_genes <- read.table(file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "significant_genes_table.tsv"), header = TRUE)

# Load pvalues (for ordering)
pvals_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "all_pvalues_with_gene_name.tsv"), as.is = TRUE, header = TRUE)

# Wrangling
#-------------------------------------------------------------------------------
# Prepare office/word data
wright_table <- docx_extract_tbl(wright_doc) # Extract first (and only) table in file
wright_list <- assign_colnames(wright_table, 1) %>% 
  dplyr::select("Gene name") %>% 
  plyr::rename(replace = c("Gene name" = "gene_name")) %>% 
  pull(gene_name)

gomez_list <- gomez_table %>% 
  plyr::rename(replace = c("SYMBOL" = "gene_name")) %>% 
  pull(gene_name)

sig_genes_list <- pull(sig_genes, gene_name)

# Create common table
#-------------------------------------------------------------------------------
all_genes <- union(wright_list, union(gomez_list, sig_genes_list))

# Build the table using tibble and mutate
full_gene_table <- tibble(Gene_name = all_genes) %>%
  mutate(
    Significant = if_else(Gene_name %in% sig_genes_list, "Yes", "No"),
    Wright = if_else(Gene_name %in% wright_list, "Yes", "No"),
    Gomez_Lopez = if_else(Gene_name %in% gomez_list, "Yes", "No"),
    Significant = factor(Significant, levels = c("Yes", "No"))
) %>%
  left_join(pvals_df, by = c("Gene_name" = "gene_name")) %>% 
  arrange(pvalue_combined, Significant) %>% 
  mutate(pvalue_discovery = signif(pvalue_discovery, 3), 
         pvalue_replication = signif(pvalue_replication, 3), 
         pvalue_combined = signif(pvalue_combined, 3))

# Find overlap 
#-------------------------------------------------------------------------------
# Intersection with Wright
cat(intersect(wright_list, sig_genes_list), sep = "\t")

# Intersection with Gomez-Lopes
cat(intersect(gomez_list, sig_genes_list), sep = "\t")

# Intersection with both
cat(intersect(intersect(wright_list, sig_genes_list), gomez_list), sep = "\t")

# Venn diagram
#-------------------------------------------------------------------------------
venn_data <- list("Wright et al.\n(2023)" = wright_list, "Gomez-Lopez\net al.\n(2019)" = gomez_list, "Robustly Replicated Transcripts" = sig_genes_list)
venn_plot <- ggVennDiagram(venn_data, label_size = 6, set_color = "grey40", edge_size = 0.2, set_size = 4.5, label = "count") + 
  scale_fill_gradient(low = "cornsilk", high = "cornsilk3") 

# Output
#-------------------------------------------------------------------------------
# Write 
write.table(x = full_gene_table, file = file.path(path_results, "100_data_exploration", "Paper_overlap_table.tsv"), quote = FALSE, sep = "\t")


# Save venn plot
png(file.path(path_results, "100_data_exploration", "Paper_overlap_venn.png"))
print(venn_plot)
dev.off()

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "100_data_exploration", "session_info", "100d_session_info.txt"))

