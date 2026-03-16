#!/usr/bin/env Rscript

#************************************
# FUNCTIONS FOR MODELING
#
# DESCRIPTION: This script contains the function for the linear mixed effects
# models used in the project to compare gene expression during pregnancy to post
# partum and cell type proportions during pregnancy compared to post-partum.
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************


# Imports and initialization
#-------------------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(performance)

# Linear mixed effects model with only id as random intercept
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size and a sample
# normalization factor. Individual is random intercept in model.
# Run for one gene at a time
run_linear_mixed_effects_model <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  mod0 <-       lme4::lmer(count ~ log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  
  # Obtain coeffecients and p-values
  coef <- try(summary(mod1)$coef[,c(1,2)], silent=TRUE)
  if (inherits(coef, "try-error")) {
    coef <- NA
  }
  
  # Compare models and extract pvalue
  anova_results <- anova(mod0, mod1)
  
  pval <- anova_results$`Pr(>Chisq)`[2]
  if(is.null(pval)){
    pval <- NA
  }
  
  # Gather all in df
  model_results_df <- as.data.frame(coef) %>%
    rownames_to_column(var = "coef_effect") %>%
    mutate(gene = rownames(counts)[gene_number],
           singularity = isSingular(mod1),
           pvalue = pval)
  
  return(model_results_df)
}



# Linear mixed effects model for cell type proportions with only id as random intercept 
#-------------------------------------------------------------------------------
# Uses cell type proportions and meta data including individual and gestational 
# age in weeks (adjusted)

run_linear_mixed_effects_model_cell_props <- function(cell, props_design){
  
  print(cell)
  # Prepare data for model containing counts and covariates across all sample for given gene
  prop_data <- props_design %>% 
    select(id, ga_week_adj) %>% 
    mutate(selected_cell = pull(props_design, cell))
  
  if (sum(pull(props_design, cell)) > 0){ # Only if cell prop has non been computed to zero
    # Run model
    mod1 <-       lme4::lmer(selected_cell ~ ga_week_adj + (1 | id), data=prop_data)
    mod0 <-       lme4::lmer(selected_cell ~ (1 | id), data=prop_data)
    
    # Obtain coeffecients and p-values
    coef <- try(summary(mod1)$coef[,c(1,2)], silent=TRUE)
    if (inherits(coef, "try-error")) {
      coef <- NA
    }
    
    # Compare models and extract pvalue
    anova_results <- anova(mod0, mod1)
    
    pval <- anova_results$`Pr(>Chisq)`[2]
    if(is.null(pval)){
      pval <- NA
    }
    
    # Gather all in df
    model_results_df <- as.data.frame(coef) %>%
      rownames_to_column(var = "coef_effect") %>%
      mutate(model_cell_type = cell,
             singularity = isSingular(mod1),
             pvalue = pval)
    
    return(model_results_df)
  }

}














