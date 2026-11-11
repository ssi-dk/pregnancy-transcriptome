#!/usr/bin/env Rscript

#************************************
# LOG FOLD CHANGE FOR MOST SIGNIFICANT GENES
#
# DESCRIPTION: This scripts plots the log fold change estimates obtained from 
# lme models that are included in different figures in the paper
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(ggpubr)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))
source(file.path(path_project, "scripts", "functions_plotting.R"))

# Load model coeffecient for module plots
model_coef_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Model pvalues
pvals_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "significant_genes_table.tsv"), as.is = TRUE, header = TRUE, row.names = 1)


# Data wrangling
#-------------------------------------------------------------------------------
# Wrangle data 
log_fold_pval_data <- pvals_df %>%
  left_join(model_coef_b2, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))),
    CI_lower = Estimate - 1.96 * Std..Error,
    CI_upper = Estimate + 1.96 * Std..Error,
    gene_name_p_value = paste0(gene_name, ", p-value: ", as.character(format(pvalue_discovery, digits = 3))))


# Create plot for highlighted genes with different trajectories in heatmap
#-------------------------------------------------------------------------------
OLFM4_plot <- plot_gene_prediction(log_fold_pval_data, "OLFM4", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ABCG1_plot <- plot_gene_prediction(log_fold_pval_data, "ABCG1", "#8968CD") + ylab("") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ANXA3_plot <- plot_gene_prediction(log_fold_pval_data, "ANXA3", "#8968CD") + ylab("") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

heatmap_highlight_plots <- do.call(ggarrange, c(OLFM4_plot, ABCG1_plot, ANXA3_plot, ncol = 3, nrow = 1, align="hv"))#, common.legend = TRUE, legend="right"))

# Create plot for STRING network plots
#-------------------------------------------------------------------------------
MMP8_plot <- plot_gene_prediction(log_fold_pval_data, "MMP8", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
TRIB2_plot <- plot_gene_prediction(log_fold_pval_data, "TRIB2", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
STOM_plot <- plot_gene_prediction(log_fold_pval_data, "STOM", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
H1F0_plot <- plot_gene_prediction(log_fold_pval_data, "H1F0", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
TAGAP_plot <- plot_gene_prediction(log_fold_pval_data, "TAGAP", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
MSR1_plot <- plot_gene_prediction(log_fold_pval_data, "MSR1", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ANXA3_plot <- plot_gene_prediction(log_fold_pval_data, "ANXA3", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
CD177_plot <- plot_gene_prediction(log_fold_pval_data, "CD177", "#8968CD") + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))


# Output
#-------------------------------------------------------------------------------
# Combined plot for heatmap
png(file.path(path_results, "100_data_exploration", "figures", "sig_genes_heatmap_highligt_logfold_plots.png"), height = 500, width = 2000, units = "px", res=120)
print(heatmap_highlight_plots )
dev.off()

# STRING network plots
png(file.path(path_results, "100_data_exploration", "figures", "MMP8_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(MMP8_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "TRIB2_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(TRIB2_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "STOM_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(STOM_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "H1F0_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(H1F0_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "TAGAP_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(TAGAP_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "MSR1_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(MSR1_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "ANXA3_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(ANXA3_plot)
dev.off()

png(file.path(path_results, "100_data_exploration", "figures", "CD177_STRING_network_logfold.png"), height = 800, width = 900, units = "px", res=120)
print(CD177_plot)
dev.off()

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "100_data_exploration", "session_info", "100c_session_info.txt"))



