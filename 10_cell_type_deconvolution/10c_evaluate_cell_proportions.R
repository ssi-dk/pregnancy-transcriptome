#!/usr/bin/env Rscript

#************************************
# DECONVOLUTION ANALYSIS
#
# DESCRIPTION: This script visualizes the cell type proportions obtained via
# CIBERSORTx, frations, docker version
# 
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization 
#-------------------------------------------------------------------------------
# Libraries 
library(tidyverse)
library(factoextra) # For PCA visualization

# Functions
source('/PATH/scripts/functions_models.R')

# Path
path <- '/PATH/data/'

# Load cell type proportions
cell_prop <- read.table(file = paste0(path, "10_cell_type_deconvolution/CIBERSORTx_Adjusted.txt"), as.is = TRUE, header = TRUE, sep = "\t")

# Load design data 
design_data <- design_data_full <- read.table(file = paste0(path, "04_normalization/design_finished.tsv"), as.is = TRUE, header = TRUE)

# Count data
counts_b2 <- read.table(file = paste0(path, "03_gene_filtration/b2_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Data wrangling
#-------------------------------------------------------------------------------
# Needed design columns
design_needed <- design_data %>% 
  filter(batch == "2") %>% 
  select(sample, ga_week_adj, id, lib_size, norm_factor, ribo_fraction) %>% 
  mutate(id = as.factor(id),
         id = forcats::fct_inorder(id),
         ga_week_adj = factor(ga_week_adj, levels = c("pp", "wfirst", sprintf("w%02d", seq(8,40))))) 

# Add weeks to data and pivot to long tidy format
cell_prop_with_design <- design_needed %>% 
  right_join(cell_prop, by = c("sample" = "Mixture")) 

# Wrangle data to format and reorder 
cell_prop_per_week <- cell_prop_with_design %>% 
  select(-id) %>% 
  pivot_longer(cols = -c(sample, ga_week_adj, P.value, Correlation, RMSE), names_to = "cell_type", values_to = "rel_prop") %>% 
  mutate(week_clean = case_when(
    ga_week_adj == "wfirst" ~ "4-7",
    ga_week_adj == "pp" ~ "PP",
    TRUE ~ sub("^w0?", "", ga_week_adj))) %>% 
  mutate(week_num = case_when(
    week_clean == "4-7" ~ 0,
    week_clean == "PP" ~ 100,  # Arbitrary large number to put 'pp' last
    TRUE ~ as.numeric(week_clean))) %>% 
  mutate(week_clean = fct_reorder(week_clean, week_num)) %>%  # To get correct order in plot
  filter(cell_type != "norm_factor" & cell_type != "ribo_fraction")

# Prepare data for modeling
cell_prop_week_for_model <- cell_prop_with_design %>% 
  select(-P.value, -Correlation, -RMSE)

# Extract cell types 
cell_types_list <- unique(pull(cell_prop_per_week, cell_type))

# Extract null cells
null_cell_list <- colnames(cell_prop[, colSums(column_to_rownames(cell_prop, var = "Mixture")) == 0])

# # Sanity checks! --> Look through cell type proportions
#-------------------------------------------------------------------------------
# No negative values
cell_prop_per_week %>% group_by(cell_type) %>% summarise(min_prop = min(rel_prop), max_prop = max(rel_prop), sum_cell = sum(rel_prop)) %>% print(n = 30L)

# Props are relative --> Sum to 1 across samples
cell_prop_per_week %>% group_by(sample) %>% summarise(sum = sum(rel_prop)) %>% mutate(sum = round(sum, 4)) %>% select(sum) %>% unique() # mutate due to floating-point number storing in R 

# Make plots (for all cells)
#-------------------------------------------------------------------------------
for (cell in cell_types_list){
  prop_per_cell <- cell_prop_per_week %>% 
    filter(cell_type == cell) %>% 
    ggplot(aes(x = week_clean, y = rel_prop)) +
    geom_boxplot(alpha = 0.3, fill = "lightskyblue2", outlier.shape = NA) + # Outliers appear via jitter
    geom_jitter(width = 0.2, height = 0, size = 2, color = "lightskyblue4", alpha = 0.8) +
    theme_bw() +
    labs(title = gsub("\\.", " ", cell)) + 
    xlab("Week") + ylab("Relative cell proportion (%)") +
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 20),
          title = element_text(size=24)) +
    scale_y_continuous(expand = expansion(add = c(0.001, 0), mult = c(0, 0.005))) # For proportions very close to 0
  
  print(paste(cell, "plotted"))
  
  png(paste0(path, "10_cell_type_deconvolution/cell_plots/cell_props_", cell, ".png"), height = 600, width = 2000, units = "px", res=120)
  print(prop_per_cell)
  dev.off()
  
}

# Which cell types differ significantly compared to post partum
#-------------------------------------------------------------------------------
# Run models for all cell types
model_dfs_cell_types <- lapply(cell_types_list, function(cell_type) run_linear_mixed_effects_model_cell_props(cell_type, cell_prop_week_for_model))

# Gather results
final_model_results_cell_props <- do.call(rbind, model_dfs_cell_types) %>% 
  mutate(p_adj = p.adjust(pvalue, method = "BH")) # FDR-corrected via Benjamini–Hochberg

# Print significant cell types
significant_cell_types <- final_model_results_cell_props %>%
  filter(p_adj < 0.05) %>%
  select(model_cell_type) %>% 
  filter(model_cell_type != "norm_factor" & model_cell_type != "ribo_fraction") %>% 
  unique()

print(significant_cell_types)

# Run PCA on cell type proportions
#-------------------------------------------------------------------------------
# Remove all zero columns
cell_props_clean <- column_to_rownames(cell_prop, var = "Mixture")[, colSums(column_to_rownames(cell_prop, var = "Mixture")) > 0]  

# Run PCA
pca <- prcomp(cell_props_clean, center = TRUE, scale. = TRUE) # Centering and standardization to focus on variance (not raw values) and ensuring equal contribution from all cell types

# Show % variance explained by each PC
summary(pca)  

# Make group plots (for significant cells)
#-------------------------------------------------------------------------------
# Create plot
prop_sig_cells_bar <- cell_prop_per_week %>% 
  group_by(week_clean, cell_type) %>% 
  mutate(mean_rel_prop = mean(rel_prop)) %>% 
  select(cell_type, week_clean, mean_rel_prop, week_num) %>% 
  unique() %>% 
  right_join(significant_cell_types, by = c("cell_type" = "model_cell_type")) %>% 
  mutate(cell_type = fct_reorder(cell_type, mean_rel_prop)) %>% 
  ggplot(aes(x = week_clean, y = mean_rel_prop, fill = cell_type)) +
  geom_bar(position="stack", stat="identity") + 
  theme_bw() +
  guides(fill=guide_legend(title="Cell type")) +
  xlab("Week") + ylab("Mean relative cell proportion") +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 19),
        legend.text = element_text(size=16),
        legend.title = element_text(size = 18)) 


# Output
#-------------------------------------------------------------------------------
# Model results
# Pvalues
final_model_results_cell_props %>%
  select(gene, pvalue, singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "10_cell_type_deconvolution/model_cell_prop_pvalues_b2.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Coeffecients
final_model_results_cell_props %>%
  select(-pvalue, -singularity) %>%
  unique() %>%
  write.table(file = paste0(path, "10_cell_type_deconvolution/model_cell_prop_coefs_b2.tsv"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

# Group plot for significant cells
png(paste0(path, "10_cell_type_deconvolution/cell_plots/prop_sig_cells.png"), height = 1900, width = 2100, units = "px", res=120)
print(prop_sig_cells)
dev.off()

# Group bar plot for significant cells
png(paste0(path, "10_cell_type_deconvolution/cell_plots/prop_sig_cells_rel_bar.png"), height = 1200, width = 2100, units = "px", res=120)
print(prop_sig_cells_bar)
dev.off()

