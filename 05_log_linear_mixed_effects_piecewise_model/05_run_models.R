#!/usr/bin/env Rscript

#************************************
# LINEAR MIXED EFFECTS MODEL
#
# DESCRIPTION: This scripts runs linear mixed effects models for RNAseq
# counts to detect genes differentially expressed given the post partum
# samples as reference against the different pregnancy weeks. It is run
# for batches separately and then combined (both using batch as random
# effect and not) for confirmation of p-values of significant genes
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(lmerTest)

# Functions
source('/PATH/scripts/functions_models.R')

# Path
path <- '/PATH/data/'

# (Wrangled/cleaned) count files
counts_b1 <- read.table(file = paste0(path, "03_gene_filtration/b1_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)
counts_b2 <- read.table(file = paste0(path, "03_gene_filtration/b2_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Design (keys with meta data)
design_data_full <- read.table(file = paste0(path, "04_normalization/design_finished.tsv"), as.is = TRUE, header = TRUE)


# Data wrangling
#-------------------------------------------------------------------------------
# Extract relevant data from design including batch
design_needed_combined <- design_data_full %>%
  select(sample, id, ga_week_adj, ribo_fraction, lib_size, norm_factor, batch) %>%
  mutate(id = as.factor(id),
         id = forcats::fct_inorder(id),
         ga_week_adj = factor(ga_week_adj, levels = c("pp", "wfirst", sprintf("w%02d", seq(8,40)))),
         batch = as.factor(batch))

# Extact relevant data from design excluding batch
design_needed <- design_needed_combined %>%
  select(-batch)

# Model
#-------------------------------------------------------------------------------

#### BATCH 1 (replication batch)
# Convert to matrix
counts_b1_matrix <- as.matrix(counts_b1)

design_data_b1 <- data.frame(sample = colnames(counts_b1)) %>%
  left_join(design_needed, by = "sample")


model_dfs_b1 <- lapply(1:nrow(counts_b1), function(gene_number) run_linear_mixed_effects_model(gene_number, counts_b1_matrix, design_data_b1))

final_model_results_b1 <- do.call(rbind, model_dfs_b1)


#### BATCH 2 (discovery batch)
# Convert to matrix
counts_b2_matrix <- as.matrix(counts_b2)

design_data_b2 <- data.frame(sample = colnames(counts_b2)) %>%
  left_join(design_needed, by = "sample")


model_dfs_b2 <- lapply(1:nrow(counts_b2), function(gene_number) run_linear_mixed_effects_model(gene_number, counts_b2_matrix, design_data_b2))

final_model_results_b2 <- do.call(rbind, model_dfs_b2)


#### COMBINED (common genes between the discovery and replication batch)
counts_b2_to_merge <- counts_b2 %>% 
  rownames_to_column(var = "gene") 

counts_both_matrix <- counts_b1 %>% 
  rownames_to_column(var = "gene") %>% 
  full_join(counts_b2_to_merge, by = "gene") %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix() %>% 
  na.omit() # To only keep genes common between both data sets

design_data_both <- data.frame(sample = colnames(counts_both_matrix)) %>%
  left_join(design_needed_combined, by = "sample")

model_dfs_both <- lapply(1:nrow(counts_both_matrix), function(gene_number) run_linear_mixed_effects_model_combined(gene_number, counts_both_matrix, design_data_both))

final_model_results_both <- do.call(rbind, model_dfs_both)

#### COMBINED --> But not with batch as random effect in model
design_data_both_no_batch <- data.frame(sample = colnames(counts_both_matrix)) %>%
  left_join(design_needed, by = "sample")

model_dfs_both_no_batch <- lapply(1:nrow(counts_both_matrix), function(gene_number) run_linear_mixed_effects_model(gene_number, counts_both_matrix, design_data_both_no_batch))

final_model_results_both_no_batch <- do.call(rbind, model_dfs_both_no_batch)



# Output
#-------------------------------------------------------------------------------
# Write model p values to file
final_model_results_b1 %>%
  select(gene, pvalue, singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_pvalues_b1.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

final_model_results_b2 %>%
  select(gene, pvalue, singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_pvalues_b2.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

final_model_results_both %>%
  select(gene, pvalue, singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_pvalues_combined.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

final_model_results_both_no_batch %>%
  select(gene, pvalue, singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_pvalues_combined_no_batch.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Write model coeffecients values to file
final_model_results_b1 %>%
  select(-pvalue, -singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_coef_b1.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

final_model_results_b2 %>%
  select(-pvalue, -singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_coef_b2.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

final_model_results_both %>%
  select(-pvalue, -singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_coef_combined.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

final_model_results_both_no_batch %>%
  select(-pvalue, -singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "05_log_linear_mixed_effects_piecewise_model/model_coef_combined_no_batch.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

