#!/usr/bin/env Rscript

#************************************
# FUNCTIONS
#
# DESCRIPTION: Functions scripts used to
# maked plots called by other scripts 
# in project 
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(ggfortify)


# Plotting PCA of gene counts with option of coloring with group
#-------------------------------------------------------------------------------
make_PCA_plot <- function(counts, design, group){
  
  pca_matrix <- counts %>%
    # make the "gene" column become the rownames of the table
    column_to_rownames("gene") %>%
    # coerce to a matrix
    as.matrix() %>%
    # transpose the matrix so that rows = samples and columns = genes
    t()
  
  # Make tmp normalization for PCR
  pca_matrix_norm <- pca_matrix/(rowSums(pca_matrix)/1e6)
  
  # Perform the PCA
  pca <- prcomp(pca_matrix)
  
  # Make plot
  gene_pca_plot <- autoplot(pca,
                            data = design,
                            colour = group)
  
  return(gene_pca_plot)
  
}







