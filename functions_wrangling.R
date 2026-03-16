#!/usr/bin/env Rscript

#************************************
# FUNCTIONS
#
# DESCRIPTION: Functions primarily for
# different operation of data wrangling
# called by other scripts in project
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)

# Calculate mapping summaries (NB: these scripts are not included, but calculations
# are seen below)
#-------------------------------------------------------------------------------
# Given count data, this function calculates additional summaries
# based on total number of mapped reads, unmapped reads, feature mapped reads etc.
calculate_sample_summaries <- function(counts, summary_data, ribo_genes){
  
  summary_data_updated <- summary_data %>%
    mutate(N_total = colSums(counts),
           N_feature = colSums(counts[-c(1:4),]),
           N_ribo = colSums(counts[ribo_genes,]),
           N_noRiboFeature = N_feature - N_ribo,
           N_mapped = N_noFeature + N_ambiguous + N_feature,
           ribo_fraction = N_ribo/N_noRiboFeature)
  
  return(summary_data_updated)
}


# Data filtration
#-------------------------------------------------------------------------------
# Function that based on ribo gene list removes these genes from given count data
remove_ribo_transcripts <- function(counts, ribo_genes_df){
  
  counts_no_ribo <- counts %>%
    anti_join(ribo_genes_df, by = c("gene" = "ribo_genes"))
  
  return(counts_no_ribo)
}


# Removes genes of low expression from counts data using the definition that
# considers only genes with a CPM > 0.5 in more than two samples in at least
# n_min pregnant people to be of "adequate expression" (inspired by org. code
# made by Line skotte, which was however made only with batch 2 pregnant people)
remove_lowly_expressed_genes_manual <- function(counts, design, n_min){
  
  # reformat design
  design_df <- select(design, sample, N_noRiboFeature, id, batch)
  m <- data.frame(sample = colnames(counts)[2:length(colnames(counts))]) %>%
    left_join(design_df, by = "sample")
  
  # Scale data with feature counts for each sample
  e <- t(t(select(counts, -gene))/(m$N_noRiboFeature/1e6)) # temporary normalization for filtering
  
  # Set thresholds
  dp <- 0.5 # read depth at transcript
  sp <- 3 # number of samples with observation
  
  # Identify lowly expressed genes
  get_filt <- function(gene){
    sum(sapply(unique(m$id), # NOTE: REMOVED BATCH ARGUMENT
               function(ind){sum((e[gene,m$id==ind]>dp))>=sp}
    )
    )
  }
  
  f <- sapply(1:(dim(e)[1]), function(ge){get_filt(ge)})
  
  # Keep genes where more than 0.5 CPM in at least 3 samples/time points in at least n_min pregnant people
  keep <- f>=n_min # 15 for discovery
  keep_genes <- data.frame(gene = pull(counts, gene)[keep])
  
  # Update counts
  counts_no_low_exprs <- counts %>%
    right_join(keep_genes, by = "gene")
  
  return(counts_no_low_exprs)
  
}



# Function for gene count normalization. Uses a linear model to estimate the effect
# of different covariates that could effect the technical variability in the gene
# count data. The residuals for each gene model are extracted to be used as the 
# normalized counts
extract_gene_count_residuals <- function(gene_number, counts, design){
  
  # Extract relevant data from design
  design_needed <- design_data %>%
    select(sample, id, ribo_fraction, lib_size, norm_factor) %>%
    mutate(id = as.factor(id),
           id = forcats::fct_inorder(id)) # Just one instance
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = counts[gene_number, ], # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_needed, by = "sample") %>%
    select(-sample)
  
  #### Extract residuals ####
  # We extract residuals to focus on the variation in gene expression not
  # explained by the covariates
  mod1 <- lm(count ~ log(lib_size) + norm_factor + log(ribo_fraction) + id, data=gene_data)
  resids <- mod1$residuals
  
  return(resids)
}


