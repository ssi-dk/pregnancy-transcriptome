#!/usr/bin/env Rscript

#************************************
# GATHER FINAL DATA
#
# DESCRIPTION: Creates final data overview, which is also the data used for 
# the shiny app
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
#library(biomaRt)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))
source(file.path(path_project, "scripts", "functions_plotting.R"))

# Load model coeffecient for module plots
model_coef_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")
model_coef_b1 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b1.tsv"), as.is = TRUE, header = TRUE, sep = "\t")
model_coef_combined <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_combined_no_batch.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Model pvalues
pvals_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "all_pvalues_with_gene_name.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Module data 
gene_module_WGCNA_results <- read.table(file = file.path(path_data, "06_WGCNA", "gene_module_model_estimates.tsv"), as.is = TRUE, header = TRUE)

# NB: As biomart is not working currently, we will use previously downloaded data
mart_info <- read.table(file = file.path(path_data, "100_data_exploration", "biomart_temp_data.tsv"), as.is = TRUE, header = FALSE, sep = "\t", col.names = c("gene_name", "full_info")) 

# Wrangling 
#-------------------------------------------------------------------------------
# Limit module data
gene_module_data_slim <- gene_module_WGCNA_results %>% 
  dplyr::select(gene, color) %>% 
  unique()

# Duplicated genes
dup_genes <- pvals_df %>% 
  dplyr::select(gene, gene_name) %>% 
  unique() %>% 
  group_by(gene_name) %>% 
  dplyr::summarize(count = n()) %>% 
  filter(count > 1) %>% 
  dplyr::select(gene_name)

# Create data sets for all batches
BSG_data_disc <- pvals_df %>%
  arrange(pvalue_discovery) %>% 
  left_join(model_coef_b2, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))),
    CI_lower = Estimate - 1.96 * Std..Error,
    CI_upper = Estimate + 1.96 * Std..Error) %>% 
  left_join(gene_module_data_slim, by = "gene") %>% 
  dplyr::select(-Std..Error, -coef_effect, -gene) %>% 
  mutate(pvalue_discovery = signif(pvalue_discovery, digits = 3),
         pvalue_replication = signif(pvalue_replication, digits = 3),
         pvalue_combined = signif(pvalue_combined, digits = 3)) %>% 
  anti_join(dup_genes, by = "gene_name") %>% 
  filter(!is.na(pvalue_replication) & !is.na(pvalue_combined)) 

BSG_data_rep <- pvals_df %>%
  arrange(pvalue_replication) %>% 
  left_join(model_coef_b1, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))),
    CI_lower = Estimate - 1.96 * Std..Error,
    CI_upper = Estimate + 1.96 * Std..Error)

BSG_data_comb <- pvals_df %>%
  left_join(model_coef_combined, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))),
    CI_lower = Estimate - 1.96 * Std..Error,
    CI_upper = Estimate + 1.96 * Std..Error) %>% 
  arrange(pvalue_combined) 

# Collect data
BSG_data_comb_tmp <- BSG_data_comb %>% 
  dplyr::rename(Estimate_combined = Estimate,
         CI_lower_combined = CI_lower,
         CI_upper_combined = CI_upper) %>% 
  dplyr::select("gene_name", "week_number", "week_end", "Estimate_combined", "CI_lower_combined", "CI_upper_combined")

BSG_data_rep_tmp <- BSG_data_rep %>% 
  mutate(Estimate_replication = Estimate,
         CI_lower_replication = CI_lower,
         CI_upper_replication = CI_upper) %>% 
  dplyr::select("gene_name", "week_number", "week_end", "Estimate_replication", "CI_lower_replication", "CI_upper_replication")


BSG_data_all <- BSG_data_disc %>% 
  left_join(BSG_data_comb_tmp, by = c("gene_name", "week_number", "week_end")) %>% 
  left_join(BSG_data_rep_tmp, by = c("gene_name", "week_number", "week_end"))  


# Gene data extraction
#-------------------------------------------------------------------------------
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# # 
# gene_info <- getBM(attributes = c("hgnc_symbol", "description", "chromosome_name"),
#                    filters = "hgnc_symbol",
#                    values = unique(pull(BSG_data_all, gene_name)),
#                    mart = mart) %>% 
#   as_tibble() %>% 
#   mutate(description_clean = sub(" \\[.*", "", description),
#          full_info = paste0(description_clean, ", chromosome ", chromosome_name),
#          gene_name = hgnc_symbol) %>% 
#   dplyr::select(full_info, gene_name)



# Add to data
BSG_data_full <- BSG_data_all %>% 
  left_join(mart_info, by = "gene_name", relationship = "many-to-many") %>% 
  mutate(full_info = case_when(is.na(full_info) ~ "-",
                               TRUE ~ sub("[\t\r\n].*$", "", full_info))) |> # Some have values added to description, so we remove these
  unique()

# Only for discovery
BSG_data_discovery <- BSG_data_full %>% 
  dplyr::select(pvalue_discovery, gene_name, Estimate, week_number, week_end, CI_lower, CI_upper, color, full_info)

# Extract source datasets
BSG_data_clean_WGCNA_model_results <- BSG_data_full %>% 
  dplyr::select(contains("pvalue"), gene_name, contains("Estimate"), contains("week"), color)

BSG_data_clean_WGCNA_model_results_with_CI <- BSG_data_full %>% 
  dplyr::select(contains("pvalue"), gene_name, contains("Estimate"), contains("week"), color, contains("CI"))

# Extract only pvalues and reformat pvalues
BSG_data_clean_WGCNA_pvals <- BSG_data_clean_WGCNA_model_results_with_CI %>% 
  dplyr::select(gene_name, color, contains("pvalue")) %>% 
  unique() %>% 
  left_join(dplyr::select(pvals_df, gene, gene_name), by = "gene_name") %>% 
  arrange(pvalue_combined) %>% 
  mutate(pvalue_discovery = sprintf("%.2e", pvalue_discovery),
         pvalue_replication = sprintf("%.2e", pvalue_replication),
         pvalue_combined = sprintf("%.2e", pvalue_combined)) 

# Just significant genes (that have been assigned to module)
BSG_data_clean_WGCNA_pvals_sig <- BSG_data_clean_WGCNA_pvals %>% 
  filter(!is.na(color))

# Output
#-------------------------------------------------------------------------------
# Full data for shiny
write.table(x = BSG_data_full, file = file.path(path_results, "100_data_exploration", "tables", "BSG_full.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = BSG_data_discovery, file = file.path(path_results, "100_data_exploration", "tables", "BSG_discovery.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Cleaned data for supplementary with and without CI
write.table(x = BSG_data_clean_WGCNA_model_results, file = file.path(path_results, "100_data_exploration", "tables", "BSG_clean_data_all_models_and_module.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = BSG_data_clean_WGCNA_model_results_with_CI, file = file.path(path_results, "100_data_exploration", "tables", "BSG_clean_data_all_models_and_module_and_CI.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Cleaned data with only pvalues
write.table(x = BSG_data_clean_WGCNA_pvals, file = file.path(path_results, "100_data_exploration", "tables", "BSG_clean_pvals_and_module.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(x = BSG_data_clean_WGCNA_pvals_sig, file = file.path(path_results, "100_data_exploration", "tables", "BSG_pregnancy_associated_transcripts_clean_pvals_and_module.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)


# Create source data for selected genes
gene_list <- c("H1F0", "STOM", "MMP8", "MSR1", "TAGAP", "TRIB2", "ANXA3", "CD177", "ABCG1", "OLFM4", "ABCG1")

for (gene in gene_list){
  gene_data <- BSG_data_clean_WGCNA_model_results_with_CI %>% 
    dplyr::select(gene_name, pvalue_discovery, Estimate, contains("week"), CI_lower, CI_upper) %>% 
    filter(gene_name == gene) 
  
  write.table(x = gene_data, file = paste0(path_results, "/100_data_exploration", "/tables/", gene, "_model_logfold_and_pvalue_with_CI.tsv"), quote = TRUE, sep = "\t", row.names = FALSE, col.names = TRUE)
              
}
# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "100_data_exploration", "session_info", "100e_session_info.txt"))





