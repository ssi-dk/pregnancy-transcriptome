#!/usr/bin/env Rscript

#************************************
# DATA FILTRATION
#
# DESCRIPTION: This scripts filters our count data to remove lowly
# expressed genes and ribosomal genes according to database and
# genes correlates to ribosomal genes.
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(edgeR)
library(ggfortify)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))
source(file.path(path_project, "scripts", "functions_plotting.R"))

# Reference with ribosomal RNAs
ribo <- read.table(file = file.path(path_input, "gencode.v28.primary_assembly.annotation.rRNAs.txt"), as.is = TRUE, header = FALSE)[,1]

# (Wrangled/cleaned) count files
counts_b1 <- read.table(file = file.path(path_data, "01_data_wrangling", "b1_clean_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)
counts_b2 <- read.table(file = file.path(path_data, "01_data_wrangling", "b2_clean_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Design (keys with meta data)
design_data <- read.table(file = file.path(path_data, "02_mapping_summaries", "design_data_with_mapping_summaries.tsv"), as.is = TRUE, header = TRUE)

# Create output file with filtration stats
outputFileFiltration <- file(file.path(path_results, "03_gene_filtration", "n_genes_removed_during_filtration.txt"))


# Data wrangling
#-------------------------------------------------------------------------------
# Wrangle count data and remove summary lines
counts_both <- counts_b1 %>%
  mutate(gene = rownames(counts_b1)) %>%
  full_join(mutate(counts_b2, gene = rownames(counts_b2)), by = "gene") %>%
  filter(!grepl("N_", gene)) %>%
  select(gene, everything())

ribo_genes <- data.frame(ribo_genes = ribo[ribo%in%pull(counts_both, gene)])


# Data filtration
#-------------------------------------------------------------------------------
status = c()

for (batch_number in c(1,2)){
  
  # Extract counts
  design_batch <-  filter(design_data, batch == as.character(batch_number))
  samples_batch <- design_batch %>% pull(sample)
  counts_batch <- counts_both %>% select(gene, all_of(samples_batch))
  
  # Update status for final printing in file
  status <- c(status, paste("BATCH", as.character(batch_number)))
  
  status <- c(status, paste("Number of genes prior to filtration:",
                            nrow(counts_both)))
  
  # Make and save PCA plot before any filtration
  pca_plot_no_filtration <- make_PCA_plot(counts_batch, design_batch, 'library')
  png(paste0(path_results, "/03_gene_filtration/PCA_before_any_filtration_b", as.character(batch_number), ".png"))
  print(pca_plot_no_filtration)
  dev.off()
  
  # Make and save PCA plot before any filtration colored according to ribo content
  design_ribo <- design_batch %>% 
    mutate("Ribo_class" = cut_number(ribo_fraction, n = 3, labels = c("low", "medium", "high")))
  
  pca_plot_no_filtration_ribo <- make_PCA_plot(counts_batch, design_ribo, 'Ribo_class')
  png(paste0(path_results, "/03_gene_filtration/PCA_before_any_filtration_ribo_color_b", as.character(batch_number), ".png"))
  print(pca_plot_no_filtration_ribo)
  dev.off()
  
  
  #### REMOVE RIBO TRANSCRIPTS ####
  print(paste0("Removing ribo, batch ", as.character(batch_number)))
  counts_no_ribo <- remove_ribo_transcripts(counts_batch, ribo_genes)
  
  status <- c(status, paste0("Number of genes after removing ribo transcripts: ",
                             nrow(counts_no_ribo),
                             " (",
                             nrow(counts_batch)-nrow(counts_no_ribo),
                             " genes removed)"))

  # Write removed genes into file
  genes_removed_ribo <- counts_batch %>% 
    anti_join(counts_no_ribo, by = "gene") %>%
    select(gene) 
  write.table(x = genes_removed_ribo, file = paste0(path_data, "/03_gene_filtration/genes_removed_during_ribo_transcript_filtration_b", batch_number, ".tsv"), quote = FALSE, row.names = FALSE)


  # Make and save PCA plot
  pca_plot_ribo_transcripts <- make_PCA_plot(counts_no_ribo, design_batch, 'library')
  png(paste0(path_results, "/03_gene_filtration/PCA_after_removal_of_ribo_transcripts_b", as.character(batch_number), ".png"))
  print(pca_plot_ribo_transcripts)
  dev.off()
  
  # Make and save PCA plot colored according to ribo
  pca_plot_ribo_transcripts_ribo <- make_PCA_plot(counts_no_ribo, design_ribo, 'Ribo_class')
  png(paste0(path_results, "/03_gene_filtration/PCA_after_removal_of_ribo_transcripts_ribo_color_b", as.character(batch_number), ".png"))
  print(pca_plot_ribo_transcripts_ribo)
  dev.off()
  
  #### REMOVE LOWLY EXPRESSED GENES ####
  print(paste0("Removing low exprs, batch ", as.character(batch_number)))
  low_exprs_thresholds <- c(floor(0.75*10), floor(0.75*21))
  #counts_no_ribo_low_expr <- remove_lowly_expressed_genes_edgeR(counts_no_ribo, design_batch)
  # alternative
  counts_no_ribo_low_expr <- remove_lowly_expressed_genes_manual(counts_no_ribo, design_data, low_exprs_thresholds[batch_number])
  
  status <- c(status, paste0("Number of genes after removing lowly expressed genes: ",
                             nrow(counts_no_ribo_low_expr),
                             " (",
                             nrow(counts_no_ribo)-nrow(counts_no_ribo_low_expr),
                             " genes removed)"))
  
  # Write removed genes into file
  genes_removed_low_expr <- counts_no_ribo %>% 
    anti_join(counts_no_ribo_low_expr, by = "gene") %>%
    select(gene) 
  write.table(x = genes_removed_low_expr, file = paste0(path_data, "/03_gene_filtration/genes_removed_during_low_expression_filtration_b", batch_number, ".tsv"), quote = FALSE, row.names = FALSE)

  # Make and save PCA plots
  pca_plot_after_low_genes <- make_PCA_plot(counts_no_ribo_low_expr, design_batch, 'library')
  png(paste0(path_results, "/03_gene_filtration/PCA_after_low_genes_removal_b", as.character(batch_number), ".png"))
  print(pca_plot_after_low_genes)
  dev.off()
  
  # Extract final genes
  assign(paste0("batch_", as.character(batch_number), "_genes"), pull(counts_no_ribo_low_expr, gene))
  
}

#### UPDATE ORIGINAL DATA ####
counts_filtered_b1 <- counts_b1 %>%
  mutate(gene = rownames(counts_b1)) %>%
  filter(gene %in% batch_1_genes) %>%
  select(-gene)

counts_filtered_b2 <- counts_b2 %>%
  mutate(gene = rownames(counts_b2)) %>%
  filter(gene %in% batch_2_genes) %>%
  select(-gene)


# Output
#-------------------------------------------------------------------------------
# Write into and close output file with filtration stats
writeLines(status, con = outputFileFiltration)
close(outputFileFiltration)

# Writing cleaned counts files
write.table(x = counts_filtered_b1, file = file.path(path_data, "03_gene_filtration", "b1_filtered_counts.tsv"), row.names = TRUE, quote = FALSE, sep = "\t")
write.table(x = counts_filtered_b2, file = file.path(path_data, "03_gene_filtration", "b2_filtered_counts.tsv"), row.names = TRUE, quote = FALSE, sep = "\t")


# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "03_gene_filtration", "session_info", "03_session_info.txt"))
