#!/usr/bin/env Rscript

#************************************
# DECONVOLUTION ANALYSIS
#
# DESCRIPTION: This prepares data for cell deconvolution analysis with CIBERSORTx
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

# Load ENSG to gene symbol data
gencode_gtf <- rtracklayer::import("PROJECT_PATH/reference/gencode.v28.primary_assembly.annotation.gtf")

# Load count data (for discovery data)
counts_b2 <- read.table(file = file.path(path_data, "03_gene_filtration", "b2_filtered_counts.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# Load LM22 signature matrix (for gene overlap check up)
LM22 <- read.table(file = file.path(path_data, "08_cell_type_deconvolution", "LM22.txt"), as.is = TRUE, header = TRUE, sep = "\t")

# Design (keys with meta data)
#design_data <- read.table(file = file.path(path_data, "02_mapping_summaries/design_data_v2.tsv"), as.is = TRUE, header = TRUE) # CHECK WHY WE USE THIS? CAN WE USE FINISHED INSTEAD?
design_data <- read.table(file = file.path(path_data, "04_normalization", "design_finished.tsv"), as.is = TRUE, header = TRUE)


# Data wrangling: Preparation
#-------------------------------------------------------------------------------
# Obtain ENSG and gene symbol/id from gtf file
ENSG_gene_symbol <- gencode_gtf %>% 
  as.data.frame() %>%
  select(gene_id, gene_name) %>%
  dplyr::rename(gene = gene_id) %>%
  unique()

# Change to gene symbol and remove redundancy
counts_b2_gene_symbol <- counts_b2 %>% 
  rownames_to_column(var = "gene") %>%
  left_join(ENSG_gene_symbol, "gene") %>% 
  select(-gene) 

design_data_b2 <- design_data %>% 
  filter(batch == "2")

# Data wrangling: Redundancy evaluation
#-------------------------------------------------------------------------------
# Check overlap of redundant genes with signature matrix
LM22_genes <- LM22 %>% 
  select("Gene.symbol") 

# Find redundant genes
redundant_genes <- counts_b2_gene_symbol %>% 
  select(gene_name) %>% 
  group_by(gene_name) %>% 
  summarise(gene_count = n()) %>% 
  filter(gene_count > 1) %>% 
  select(gene_name) 

# Check if any overlap 
inner_join(LM22_genes, redundant_genes, by = c("Gene.symbol" = "gene_name")) # No overlap

# Check how many of 547 ref genes we have in our data
inner_join(LM22_genes, select(counts_b2_gene_symbol, gene_name), by = c("Gene.symbol" = "gene_name")) # 506

# Data wrangling: Cleaning
#-------------------------------------------------------------------------------
# Remove redundant genes (as they are not in LM22 anyways)
counts_b2_CIBERSORTx <- counts_b2_gene_symbol %>% 
  anti_join(redundant_genes, by = "gene_name") %>% 
  column_to_rownames("gene_name")

# Data wrangling: CPM NORMALIZATION
#-------------------------------------------------------------------------------
# Obtain total sample sizes for normalization
total_count <- setNames(design_data_b2$N_noRiboFeature, design_data_b2$sample)
total_count <- total_count[colnames(counts_b2_CIBERSORTx)]

# Convert to DGEList
dge <- DGEList(counts = counts_b2_CIBERSORTx, lib.size = total_count)

# Convert to CPM
cpm_counts <- cpm(dge, log = FALSE)

# Make column name 
cpm_counts_b2_CIBERSORTx <- cpm_counts %>%
  as.data.frame() %>%
  rownames_to_column("Gene")

# Output
#-------------------------------------------------------------------------------
# Write wrangled data to file
write.table(x = cpm_counts_b2_CIBERSORTx, file = file.path(path_data, "08_cell_type_deconvolution", "mixture_file.tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "08_cell_type_deconvolution", "session_info", "08a_session_info.txt"))



