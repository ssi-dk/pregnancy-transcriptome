#!/usr/bin/env Rscript

#************************************
# LOG FOLD CHANGE FOR MOST SIGNIFICANT GENES
#
# DESCRIPTION: This scripts plots the log fold change estimates obtained from 
# lme models for the genes found to be most significantly differently expressed
# during pregnancy compared to the post partum samples.
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# OBS: Needed: module load openssl/1.1.1s

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))

# Load model coeffecient for module plots
model_coef_b2 <- read.table(file = file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Model pvalues
pvals_df <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "significant_genes_table.tsv"), as.is = TRUE, header = TRUE, row.names = 1)


# Data wrangling
#-------------------------------------------------------------------------------
n_genes <- 6
# Filter to obtain 5 most significant genes
data_most_sig_genes <- pvals_df %>%
  arrange(pvalue_discovery) %>% 
  filter(row_number() %in% seq(1,n_genes)) %>%
  left_join(model_coef_b2, by = "gene") %>%
  filter(coef_effect != "(Intercept)", coef_effect != "log(lib_size)", coef_effect != "norm_factor", coef_effect != "log(ribo_fraction)") %>%
  mutate(
    week_number = case_when(coef_effect == "ga_week_adjwfirst" ~ 4,
                            coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))) - 1,
    week_end = case_when(coef_effect == "ga_week_adjwfirst" ~ 7,
                         coef_effect != "ga_week_adjwfirst" ~ as.numeric(gsub("[^0-9]", "", coef_effect))),
    CI_lower = Estimate - 1.96 * Std..Error,
    CI_upper = Estimate + 1.96 * Std..Error,
    gene_name_p_value = paste(gene_name, "\np-value:", as.character(format(pvalue_discovery, digits = 3))))

print(head(data_most_sig_genes))

# Log fold change plots
#-------------------------------------------------------------------------------
# Plot log fold change compared to post partum for the 5 most significant genes
linewidth_size = c("Estimate" = 2, "Confidence Interval" = 0.5)

fold_change_plots <- data_most_sig_genes %>%
  mutate(gene_name_p_value = fct_reorder(gene_name_p_value, pvalue_discovery)) %>% # For correct ordering of genes 
  ggplot() +
  geom_segment(aes(x = week_number, xend = week_end, y = Estimate, yend = Estimate, color = gene_name_p_value, linewidth = "Estimate")) + 
  geom_segment(aes(x = week_number, xend = week_end, y = CI_lower, yend = CI_lower, color =  gene_name_p_value, linewidth = "Confidence Interval"), alpha = 0.5) + 
  geom_segment(aes(x = week_number, xend = week_end, y = CI_upper, yend = CI_upper, color = gene_name_p_value, linewidth = "Confidence Interval"), alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Pregnancy week", y = "Log(fold change) to PP level") +
  scale_linewidth_manual(values = linewidth_size, name = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.position = "bottom",
        strip.background = element_rect(fill = "aliceblue", colour = "transparent", size = 1),
        strip.text = element_text(size = 18)) +
        #legend.position = "right", #"bottom",
        #legend.title = element_blank(), 
        #legend.text = element_text(size = 18),
        #legend.spacing.y = unit(0.5, "cm")) +
        #legend.box.background = element_rect(fill = "transparent", color = "transparent"))
  #guides(color = guide_legend(byrow = TRUE)) +
  facet_wrap(~gene_name_p_value, ncol = n_genes/2) +
  guides(color = "none")



# Plot selected gene (MMP8)
#-------------------------------------------------------------------------------
# Plot log fold change compared to post partum for the 5 most significant genes
linewidth_size_one_gene = c("Estimate" = 3, "Confidence Interval" = 1)

MMP8_title <- data_most_sig_genes %>% filter(gene_name == "MMP8") %>% pull(gene_name_p_value) %>% as.character() 
MMP8_title <- gsub(" \n", ", ", MMP8_title)

MMP8_plot <- data_most_sig_genes %>%
  filter(gene_name == "MMP8") %>%
  mutate(gene_name_p_value = fct_reorder(gene_name_p_value, pvalue_discovery)) %>% # For correct ordering of genes 
  ggplot() +
  geom_segment(aes(x = week_number, xend = week_end, y = Estimate, yend = Estimate, linewidth = "Estimate"), color = "coral4") + 
  geom_segment(aes(x = week_number, xend = week_end, y = CI_lower, yend = CI_lower, linewidth = "Confidence Interval"), color =  "coral3", alpha = 0.5) + 
  geom_segment(aes(x = week_number, xend = week_end, y = CI_upper, yend = CI_upper, linewidth = "Confidence Interval"), color = "coral3", alpha = 0.5) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(x = "Pregnancy week", y = "Log(fold change) to PP level", title = MMP8_title) +
  scale_linewidth_manual(values = linewidth_size_one_gene, name = "") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 23),
        legend.position = "bottom",
        strip.background = element_rect(fill = "aliceblue", colour = "transparent", size = 1),
        strip.text = element_text(size = 18)) +
  guides(color = "none")



# Output
#-------------------------------------------------------------------------------
# Save top gene plot
png(file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "LogFoldChange_plot_most_sig_genes.png"), width = 1400, height = 800)
print(fold_change_plots)
dev.off()  

# Save MMP8 plot
png(file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "LogFoldChange_plot_MMP8.png"), width = 900, height = 600)
print(MMP8_plot)
dev.off()

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "05_log_linear_mixed_effects_piecewise_model", "session_info", "05c_session_info.txt"))
