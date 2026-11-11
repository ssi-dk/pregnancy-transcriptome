#!/usr/bin/env Rscript

#************************************
# FUNCTIONS
#
# DESCRIPTION: Functions primarily for different operation of data wrangling
# called by other scripts in project
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(sessioninfo)


# Function to write sessioninfo of a script
#-------------------------------------------------------------------------------
save_session_info <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  
  con <- file(path, open = "wt")
  on.exit(close(con))
  
  writeLines(
    capture.output(sessioninfo::session_info(pkgs = "loaded")),
    con
  )
}

# Calculate mapping summaries
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

# Removes genes of low expression from count data with edgeR's filterByExpr()
# based on a grouping the samples where one group consists of all "pre partum"
# samples, i.e. samples taken before the women gave birth, and the second group
# the post partum samples, i.e. the samples taken after the women gave birth
remove_lowly_expressed_genes_edgeR <- function(counts, design){
  
  # Make DGElist object with edgeR
  y = DGEList(counts = counts)
  
  # Make group according to whether the sample is pregnant or not pregnant (pre vs. post partum)
  ga_week_adj <- select(design, sample, ga_week_adj)
  
  group <- data.frame(sample = colnames(counts)[2:length(colnames(counts))]) %>%
    left_join(ga_week_adj, by = "sample") %>%
    mutate(week_binary = case_when(ga_week_adj == "pp" ~ 1, # Basically pregnant or not
                                   ga_week_adj != "pp" ~ 2)) %>%
    pull(week_binary) %>%
    factor()
  
  # Obtain genes that are not of low expression with edgeR
  keep <- filterByExpr(y, group = group)
  
  # Update object
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # Update counts to keep only genes of adequate expression
  counts_no_low_exprs <- counts %>%
    right_join(y$genes, by = "gene")
  
  return(counts_no_low_exprs)
  
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

# Creates a test to compare a genes association to the ribo ratio (number of
# reads mapped to ribo transcript over number of total mapped genes per sample)
# to create a proxy for ribo association. It then removes the genes which have
# a low (<1) test statistic from the count data (inspired by org. code from
# Line Skotte)
remove_ribo_associated_genes <- function(counts, design, output_path){
  
  # Extract meta data only for samples in count data and correct order to follow count data
  design_selected <- data.frame(sample = colnames(select(counts, -gene))) %>%
    left_join(design, by = "sample")
  rownames(design_selected) <- design_selected$sample
  
  # Reformat data
  counts_with_names <- counts %>%
    column_to_rownames(var = "gene")
  counts_norm_by_feature_2 <- t(t(counts_with_names)/design_selected$N_feature) # Note: possibly like this because we reordered design with design_selected
  
  ribo_test <- sapply(X = 1:(nrow(counts_norm_by_feature_2)),
                      FUN = function(i){
                        summary(lm(counts_norm_by_feature_2[i,] ~ design_selected$N_ribo/design_selected$N_mapped))$coef[2,3:4] # Extract test and p-value for N_ribo coeffecients
                        # Note: returns coefficients for both N_ribo and interaction between N_ribo and N_mapped, we want N_ribo to get the estimated change of this when
                        # N_mapped is held constant
                      }
  )
  
  # Reformat data
  colnames(ribo_test) <- rownames(counts)
  rownames(ribo_test) <- c("t-value", "p-value")
  ribo_test_t <- t(ribo_test)
  
  # Keep only genes with test statistics less than 1
  keep <- (ribo_test_t[,1]<1)
  keep_genes <- data.frame(gene = pull(counts, gene)[keep])
  
  # Write file with ribo associate genes and make qqplot of p-values
  ribo_asso_genes <- data.frame(gene = pull(counts, gene)[!keep])
  
  write.table(x = ribo_asso_genes,
              file = paste0(output_path, "ribo_associated_genes.tsv"),
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  png(file = paste0(output_path, "qqplot_ribo_asso_p_value.png"))
  qqp(ribo_test_t[,2], ci=FALSE)
  dev.off()
  
  # Update counts
  counts_no_ribo_asso <- counts %>%
    right_join(keep_genes, by = "gene")
  
  return(counts_no_ribo_asso)
  
}



# Function for gene count normalization. Uses a linear model to estimate the effect
# of different covariates that could effect the technical variability in the gene
# count data. The residuals for each gene model are extracted to be used as the 
# normalized counts
extract_gene_count_residuals <- function(gene_number, counts, design){
  
  # Extract relevant data from design
  design_needed <- design %>%
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


