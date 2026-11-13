#!/usr/bin/env Rscript

#************************************
# MAPPING SUMMARIES
#
# DESCRIPTION: This script creates summary statistics on the mapping of the 
# RNAseq data
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

# Reference with ribosomal RNAs
ribo <- read.table(file = file.path(path_input, "gencode.v28.primary_assembly.annotation.rRNAs.txt"), as.is = TRUE, header = FALSE)[,1]

# (Wrangled/cleaned) count files
counts_b1 <- read.table(file = file.path(path_data, "01_data_wrangling", "b1_clean_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)
counts_b2 <- read.table(file = file.path(path_data, "01_data_wrangling", "b2_clean_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Design (keys with meta data)
design <- read.table(file = file.path(path_data, "01_data_wrangling", "design_data.tsv"), as.is = TRUE, header = TRUE)


# Data wrangling
#-------------------------------------------------------------------------------
# Extract only the ribosomes we have in our data
ribo <- ribo[ribo%in%rownames(counts_b1)]

# Convert data format and join data to work with just one data set
mapping_summary_b1 <- t(counts_b1[1:4,]) %>% as.data.frame() 
mapping_summary_b2 <- t(counts_b2[1:4,]) %>% as.data.frame()

# Add additional summaries for both batches (add N_total, N_feature, N_ribo, N_noRiboFeatures and N_mapped for all samples)
mapping_summary_full_b1 <- calculate_sample_summaries(counts_b1, mapping_summary_b1, ribo)
mapping_summary_full_b2 <- calculate_sample_summaries(counts_b2, mapping_summary_b2, ribo)

# Collect data for plotting
mapping_summary_full_b1_batch <- mapping_summary_full_b1 %>% 
  mutate(batch = "Batch 1")
mapping_summary_full_both <- mapping_summary_full_b2 %>% mutate(batch = "Batch 2") %>% 
  bind_rows(mapping_summary_full_b1_batch)

# Make summary for paper 
mapping_stats <- mapping_summary_full_both %>% 
  bind_rows(., mapping_summary_full_both %>% mutate(batch = "Both")) %>% 
  group_by(batch) %>% 
  summarise(mean_reads = mean(N_total),
            perc_unique = mean(N_mapped)/mean_reads,
            perc_multimap = mean(N_multimapping)/mean_reads,
            perc_unmap = mean(N_unmapped)/mean_reads, 
            total = sum(N_total)) 

# Adding mapping summary to meta data
colnames(mapping_summary_full_both)
design_with_mapping_summary <- mapping_summary_full_both %>% 
  select(-batch) %>% 
  mutate(sample = rownames(mapping_summary_full_both)) %>% 
  left_join(mutate(design, sample = paste(library, barcode, sep = "_")), by = "sample")
  

# Plotting  
#-------------------------------------------------------------------------------
# Add library info to data
mapping_summary_df <- mapping_summary_full_both %>% 
  mutate(Library = gsub("^(.*)_[^_]*_[^_]*$", "\\1", rownames(mapping_summary_full_both))) 

# Distribution of total reads
total_reads_hist <- 
 mapping_summary_df %>% 
  ggplot() +
  geom_histogram(aes(x = N_total/1e6), fill = "lightblue", color = "grey", bins = 30) +
  xlab("Total number of reads (Millions)") +
  facet_wrap(~batch) 

# Boxplot of total reads per sample_library
boxplot_reads_per_sample <- 
mapping_summary_df %>%
  ggplot() +
  geom_boxplot(aes(x = N_total/1e6, y = Library), fill = "lightblue") +
  xlab("Total number of reads (Millions)") + ylab("Library") +
  facet_wrap(~batch)

# Reads in features
reads_in_features <- 
mapping_summary_df %>% 
  ggplot() +
  geom_histogram(aes(x = N_feature/1e6), fill = "forestgreen", color = "grey", bins = 30) +
  xlab("Number of feature reads (Millions)") + 
  facet_wrap(~batch) 

# Vs plots
total_vs_feature <- mapping_summary_df %>% 
  ggplot(aes(x = N_total/1e6, y = N_feature/1e6, color = Library)) +
  geom_point(alpha = 0.5, size = 0.8) +
  xlab("Total number of reads (Millions)") + ylab("Number of feature reads (Millions)") +
  facet_wrap(~batch) 
ribo_vs_no_ribo <- mapping_summary_df %>% 
  ggplot(aes(x = N_noRiboFeature/1e6, y = N_ribo/1e6, color = Library)) +
  geom_point(alpha = 0.5, size = 0.8) +
  xlab("Number of non-ribo feature reads (Millions)") + ylab("Number of ribo reads (Millions)") +
  facet_wrap(~batch) 
total_vs_ribo <- mapping_summary_df %>% 
  ggplot(aes(x = N_total/1e6, y = N_ribo/1e6, color = Library)) +
  geom_point(alpha = 0.5, size = 0.8) +
  xlab("Total number of reads (Millions)") + ylab("Number of ribo reads (Millions)") +
  facet_wrap(~batch)

# Distribution of ratio of mapped vs multimapping
mapping_ratio <- 
mapping_summary_df %>% 
  ggplot() +
  geom_histogram(aes(x = N_mapped/N_multimapping), fill = "gold", color = "grey", bins = 40) +
  xlab("Number of mapped reads (Millions) / Number of multimapped reads (Millions)") + 
  facet_wrap(~batch) 



# Output 
#-------------------------------------------------------------------------------
# Save mapping statistics plots
png(file.path(path_results, "02_mapping_summaries", "distribution_total_reads.png"))
print(total_reads_hist)
dev.off()

png(file.path(path_results, "02_mapping_summaries", "boxplot_read_per_sample.png"))
print(boxplot_reads_per_sample)
dev.off()

png(file.path(path_results, "02_mapping_summaries", "distribution_reads_features.png"))
print(reads_in_features)
dev.off()

# Save plots collected in one
png(file.path(path_results, "02_mapping_summaries", "count_plotted_against_each_other.png"), width = 900, height = 1200, units = "px", res=120)
print(ggarrange(total_vs_feature, ribo_vs_no_ribo, total_vs_ribo, ncol = 1, nrow = 3, common.legend = TRUE, legend = "right"))
dev.off()

# Save plot with mapping ratio
png(file.path(path_results, "02_mapping_summaries", "mapping_ratio.png"))
print(mapping_ratio)
dev.off()

# Mapping summary statistics
write.table(x = mapping_stats, file = file.path(path_results, "02_mapping_summaries", "mapping_summary_statistics.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")

# Design with added mapping summaries
write.table(x = design_with_mapping_summary, file = file.path(path_data, "02_mapping_summaries", "design_data_with_mapping_summaries.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "02_mapping_summaries", "session_info", "02_session_info.txt"))
