#!/usr/bin/env Rscript

#************************************
# LINEAR MIXED EFFECTS MODEL FUNCTION
#
# DESCRIPTION: This script contains the overall functions used in the project
#
# Author: Freja Dahl Hede, FREH@dksund.dk
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


# Linear mixed effects model with only id AND batch as random intercept 
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size and a sample
# normalization factor. Individual and batch are random intercepts in model.
# Run for one gene at a time
run_linear_mixed_effects_model_combined <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id) + (1 | batch), data=gene_data)
  mod0 <-       lme4::lmer(count ~ log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id) + (1 | batch), data=gene_data)
  
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


# Linear mixed effects model with only id as random intercept and sex included as factor
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size, a sample
# normalization factor, and sex. Individual is random intercept in model.
# Run for one gene at a time
run_linear_mixed_effects_model_with_sex <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + sex_child + (1 | id), data=gene_data)
  mod0 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  
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


# Linear mixed effects model with only id as random intercept and PARITY included as factor
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size, a sample
# normalization factor, and parity. Individual is random intercept in model.
# Run for one gene at a time
run_linear_mixed_effects_model_with_sex <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + parity + (1 | id), data=gene_data)
  mod0 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  
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


# Linear mixed effects model with only id as random intercept and MOTHER'S BMI included as factor
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size, a sample
# normalization factor, and pre preg bmi. Individual is random intercept in model.
# Run for one gene at a time
run_linear_mixed_effects_model_with_sex <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + pre_preg_bmi + (1 | id), data=gene_data)
  mod0 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  
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



# Linear mixed effects model with only id as random intercept and MOTHER'S AGE included as factor
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size, a sample
# normalization factor, and the age of the mother. Individual is random intercept in model.
# Run for one gene at a time
run_linear_mixed_effects_model_with_sex <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + mom_age_birth + (1 | id), data=gene_data)
  mod0 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  
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



# Linear mixed effects model with only id as random intercept and MOTHER'S AGE included as factor
#-------------------------------------------------------------------------------
# Uses gene counts and meta data including individual, gestational age in weeks (adjusted),
# the original fraction of ribosomal RNA in sample, library/sample size, a sample
# normalization factor, and the weight of the child at birth. Individual is random intercept in model.
# Run for one gene at a time
run_linear_mixed_effects_model_with_sex <- function(gene_number, counts, design_data){
  
  # Prepare data for model containing counts and covariates across all sample for given gene
  gene_data <- data.frame(count = log(counts[gene_number, ] + 1), # Log transform --> Try scaling also
                          sample = colnames(counts)) %>%
    full_join(design_data, by = "sample") %>%
    select(-sample)
  
  # Run model
  mod1 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + birth_weight_g + (1 | id), data=gene_data)
  mod0 <-       lme4::lmer(count ~ ga_week_adj + log(lib_size) + norm_factor + log(ribo_fraction) + (1 | id), data=gene_data)
  
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




