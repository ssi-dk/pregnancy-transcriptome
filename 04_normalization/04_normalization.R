#!/usr/bin/env Rscript

#************************************
# DATA NORMALIZATION
#
# DESCRIPTION: This scripts uses edgeR to TMM normalize
# gene counts
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(edgeR)

# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))

# (Wrangled/cleaned) count files
counts_b1 <- read.table(file = file.path(path_data, "03_gene_filtration", "b1_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)
counts_b2 <- read.table(file = file.path(path_data, "03_gene_filtration", "b2_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Design (keys with meta data)
design_data <- read.table(file = file.path(path_data, "02_mapping_summaries", "design_data_with_mapping_summaries.tsv"), as.is = TRUE, header = TRUE)

# Normalization
#-------------------------------------------------------------------------------
# Make DGE object
create_norm_factor <- function(counts){
  # Create DGEList object for edgeR
  y = DGEList(counts = counts)
  
  # Calculate and extract normalization factors
  y <- calcNormFactors(object = y, method = "TMM")
  
  norm_factor <- select(y$samples, norm.factors, lib.size) %>% 
    dplyr::rename(norm_factor = norm.factors, 
           lib_size = lib.size) 
  

  return(norm_factor)
  
}

# Calculate normalization factors
norm_factor_counts_b1 <- create_norm_factor(counts_b1)
norm_factor_counts_b2 <- create_norm_factor(counts_b2)
 
# Add to design
design_with_norm_factors <- norm_factor_counts_b1 %>%
  bind_rows(norm_factor_counts_b2) %>% 
  rownames_to_column(var = "sample") %>% 
  #mutate(sample = rownames(norm_factor)) %>% 
  full_join(design_data, by = "sample") 


# Output 
#-------------------------------------------------------------------------------
# Write design data (i.e. keys with meta data) 
write.table(x = design_with_norm_factors, file = file.path(path_data, "04_normalization", "design_finished.tsv"), quote = FALSE, sep = "\t")


# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "04_normalization", "session_info", "04_session_info.txt"))



