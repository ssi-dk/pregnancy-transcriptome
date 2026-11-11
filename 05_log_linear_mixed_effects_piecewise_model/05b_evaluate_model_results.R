#!/usr/bin/env Rscript

#************************************
# EVALUTE SIGNIFICANT GENES
#
# DESCRIPTION: This script evaluates which genes have been found
# to be significant from the lme models
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("rtracklayer")
library(rtracklayer)
suppressMessages(library(tidyverse))


# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))

# Load model p values for each gene
pvalues_b1 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_pvalues_b1.tsv"), as.is = TRUE, header = TRUE)
pvalues_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_pvalues_b2.tsv"), as.is = TRUE, header = TRUE)
pvalues_combined <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_pvalues_combined_no_batch.tsv"), as.is = TRUE, header = TRUE)

# Load model coeffecients for each gene
model_coef_b1 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b1.tsv"), as.is = TRUE, header = TRUE, sep = "\t")
model_coef_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")
model_coef_combined <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_combined_no_batch.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Gencode v28 primary assembly GTF for gene symbol
gencode_gtf <- rtracklayer::import("PROJECT_PATH/reference/gencode.v28.primary_assembly.annotation.gtf")

# Wrangling gtf 
#-------------------------------------------------------------------------------
# Obtain ENSG and gene symbol/id from gtf file
ENSG_gene_symbol <- gencode_gtf %>% 
  as.data.frame() %>%
  select(gene_id, gene_name) %>% 
  dplyr::rename(gene = gene_id) %>%
  unique()

# Statistics
#-------------------------------------------------------------------------------
# Make data with bonferroni corrected p values for both batches and combined
pvalues_b1_bonferroni_sig <- pvalues_b1 %>%
  filter(pvalue < 0.05/nrow(pvalues_b1))

pvalues_b2_bonferroni_sig <- pvalues_b2 %>%
  filter(pvalue< 0.05/nrow(pvalues_b2))

pvalues_combined_bonferroni_sig <- pvalues_combined %>%
  filter(pvalue< 0.05/nrow(pvalues_combined)) 

# The replicated confirmed genes of significant discovery
confirmed_sig_genes <- pvalues_b2_bonferroni_sig %>%
  inner_join(pvalues_b1_bonferroni_sig, by = "gene") %>% 
  dplyr::rename(pvalue_discovery = pvalue.x,
         pvalue_replication = pvalue.y)

# Make final table of significant genes
most_sig_genes <- confirmed_sig_genes %>%
  #inner_join(pvalues_combined_bonferroni_sig, by = "gene") %>% 
  inner_join(pvalues_combined, by = "gene") %>% #  We are essentially not checking if they are significant in combined, but as it has to be lower than the sig. p-values later it will be the case nontheless
  dplyr::rename(pvalue_combined = pvalue) %>% 
  select(-contains("singularity")) %>%
  arrange(pvalue_combined) %>% 
  left_join(ENSG_gene_symbol, by = "gene") 

# Make final table all p values and gene names 
all_pvalues_ENSG <- pvalues_b2 %>%
  dplyr::rename(pvalue_discovery = pvalue) %>%
  full_join(pvalues_b1, by = "gene") %>%
  dplyr::rename(pvalue_replication = pvalue) %>%
  full_join(pvalues_combined, by = "gene") %>%
  dplyr::rename(pvalue_combined = pvalue) %>%
  select(-contains("singularity")) %>%
  arrange(pvalue_combined) %>% 
  left_join(ENSG_gene_symbol, by = "gene")


# Output
#-------------------------------------------------------------------------------
# Print statements
print("Replication, soft threshold")
pvalues_b1 %>%
  filter(pvalue < 0.05) %>%
  nrow()

print(paste0("Replication, Bonferroni (", as.character(0.05/nrow(pvalues_b1)), " for ", as.character(nrow(pvalues_b1)), " tested transcripts"))
pvalues_b1_bonferroni_sig %>%
  nrow()

print("Discovery, soft threshold")
pvalues_b2 %>%
  filter(pvalue < 0.05) %>%
  nrow()

print(paste0("Discovery, Bonferroni (", as.character(0.05/nrow(pvalues_b2)), " for ", as.character(nrow(pvalues_b2)), " tested transcripts"))
pvalues_b2_bonferroni_sig %>%
  nrow()

print("Combined, soft threshold")
pvalues_combined %>%
  filter(pvalue < 0.05) %>%
  nrow()

print(paste0("Combined, Bonferroni  (", as.character(0.05/nrow(pvalues_combined)), " for ", as.character(nrow(pvalues_combined)), " tested transcripts"))
pvalues_combined_bonferroni_sig %>%
  nrow()

print("Replication confirmed significant genes")
confirmed_sig_genes %>%
  nrow()

print("Median of p value in entire discovery batch")
print(median(pull(pvalues_b2, pvalue)))
# Evaluate significance criteria - do any genes have higher combined p-values than discovery and rep?
#-------------------------------------------------------------------------------
# Count how many genes this is the case for
n_genes_soft_threshold <- most_sig_genes %>% 
  filter(!((pvalue_discovery >= pvalue_combined) & (pvalue_replication >= pvalue_combined))) 

print("Number of genes removed with hard threshold (i.e. require that combined model p value is lower than for separate batches), i.e. the number of genes that will be filtered out with the more conservative threshold")
print(nrow(n_genes_soft_threshold))
print(n_genes_soft_threshold)


# Look into plot of example genes
#example_soft_thredshold_genes <- c("ENSG00000135218.17", "ENSG00000173559.12", "ENSG00000038427.15")
example_soft_thredshold_genes <- c("ENSG00000104447.12", "ENSG00000171791.12", "ENSG00000038427.15")

# Wrangle data for plots
replication_coef <- model_coef_b1 %>%
    mutate(model = "Replication")

discovery_coef <- model_coef_b2 %>%
    mutate(model = "Discovery")

combined_coef <- model_coef_combined %>%
    mutate(model = "Combined")

model_estimates_soft_thresholds <- replication_coef %>%
    bind_rows(discovery_coef) %>%
    bind_rows(combined_coef) %>%
    filter(gene %in% example_soft_thredshold_genes) %>% 
    filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
    mutate(
      week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                              coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
      week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                           coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect)))
    ) %>%
    left_join(most_sig_genes, by = "gene") %>%
    mutate(gene_pvalue = paste0(gene, "\nr: ", as.character(pvalue_replication), "\nd: ", as.character(pvalue_discovery), "\nc: ", as.character(pvalue_combined)))


fold_change_plots <- model_estimates_soft_thresholds %>%
ggplot() +
  geom_segment(aes(x = week_number, xend = week_end, y = Estimate, yend = Estimate, color = model), linewidth = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Replication" = "blue", "Discovery" = "red", "Combined" = "green")) +
  theme_bw() +
  labs(x = "Pregnancy week", y = "Log(fold change) to PP level") +
  theme(plot.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 14), 
        legend.position = "bottom",
        legend.title = element_blank(), 
        legend.text = element_text(size = 14)) +
  facet_wrap(~gene_pvalue, ncol = length(example_soft_thredshold_genes))



# Make final data
#-------------------------------------------------------------------------------
most_sig_genes <- most_sig_genes %>% 
  filter((pvalue_discovery > pvalue_combined) & (pvalue_replication > pvalue_combined)) 

confirmed_sig_genes <- confirmed_sig_genes %>%
  right_join(select(most_sig_genes, gene), by = "gene")

print("Number of replication confirmed and lower combined p-value confirmed significant genes (hard threshold)")
confirmed_sig_genes %>%
  nrow()

# Output data
#-------------------------------------------------------------------------------
# Save plots
png(file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "LogFoldChange_soft_thredshold_plots.png"))
print(fold_change_plots)
dev.off() 

# Result files
write.table(x = confirmed_sig_genes, file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "replicated_significant_genes_and_p_values.tsv"), quote = FALSE, sep = "\t")
write.table(x = most_sig_genes, file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "significant_genes_table.tsv"), quote = FALSE, sep = "\t")
write.table(x = all_pvalues_ENSG, file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "all_pvalues_with_gene_name.tsv"), quote = FALSE, sep = "\t")

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "session_info", "05b_session_info.txt"))
