#!/usr/bin/env Rscript

#************************************
# CREATE SIGNIFICANT GENE LOG FOLD HEAT MAP
#
# DESCRIPTION: This script makes heatmap of log fold change for the pregnancy-associated
# genes 
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(patchwork)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))

# Load model results
model_coef_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Load model pvalues for significant genes
pvals_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "significant_genes_table.tsv"), as.is = TRUE, header = TRUE, row.names = 1)


# Wrangling
#-------------------------------------------------------------------------------
# How many genes to include
n_genes <- 40

# Clean time, reorder data and select top genes
model_results_wrangled <- model_coef_b2 %>% 
  filter(grepl("ga_week_", coef_effect)) %>% 
  mutate(week_clean = case_when(
    coef_effect == "ga_week_adjwfirst" ~ "4-7",
    TRUE ~ sub("^ga_week_adjw0?", "", coef_effect))) %>% 
  mutate(week_num = case_when(
    week_clean == "4-7" ~ 0,
    week_clean == "PP" ~ 100,  # Arbitrary large number to put 'pp' last
    TRUE ~ as.numeric(week_clean))) %>% 
  mutate(week_clean = fct_reorder(week_clean, week_num)) %>% 
  right_join(pvals_df, by = "gene") %>% 
  mutate(gene_name = fct_reorder(gene_name, pvalue_combined, .desc = TRUE)) %>% 
  arrange(pvalue_combined) %>% 
  as_tibble() %>% 
  head(33*n_genes)


# Make heatmap
#-------------------------------------------------------------------------------
# Create heatmap plot 
gene_fold_change_heatmap <- model_results_wrangled %>% 
  ggplot(aes(x=week_clean, y=gene_name, fill=Estimate, width=0.95, height=0.95)) +
  geom_tile() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 13, angle = 90, hjust = 0.8, vjust = 0.5),
        axis.title = element_text(size = 15),
        legend.title.align = 0.5,
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12)) +
  scale_fill_gradient2(low="royalblue", mid = "grey90", high="firebrick", midpoint = 0, name = "Log(fold change) to PP level") +
  xlab("Gestational age (weeks)") + ylab("Gene") 
  
# Create -log p value dotplot for top genes
batch_colors <- c("Combined" = "#698B69", "Discovery" = "#8968CD", "Replication" = "#EE7600")

log_p_value_plot <- model_results_wrangled %>% 
  pivot_longer(contains("pvalue"), names_to = "batch", values_to = "pvalue") %>% 
  mutate(log_p_value = -log10(pvalue)) %>% 
  mutate(Batch = case_when(batch == "pvalue_combined" ~ "Combined",
                           batch == "pvalue_discovery" ~ "Discovery",
                           batch == "pvalue_replication" ~ "Replication")) %>% 
  ggplot(aes(x = log_p_value, y = gene_name, shape = Batch, fill = Batch, color = Batch)) +
  geom_point(size = 3) +
  scale_color_manual(values = batch_colors) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title.align = 0.5,
        legend.position = "bottom",
        legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = 12), 
        legend.title = element_blank(),
        legend.spacing = unit(0.02, 'cm')) +
  xlab("-log10(P-value)") + ylab("") 

# Create common aligned plot
combined_plot <- gene_fold_change_heatmap + log_p_value_plot +
  plot_layout(ncol = 2, widths = c(2.3, 1))

# Output
#-------------------------------------------------------------------------------
# Combined plot
png(file.path(path_results, "100_data_exploration", "figures", "sig_genes_heatmap_pvalues.png"), height = 1000, width = 1800, units = "px", res=120)
print(combined_plot)
dev.off()

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "100_data_exploration", "session_info", "100b_session_info.txt"))



