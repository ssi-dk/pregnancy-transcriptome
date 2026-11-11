#!/usr/bin/env Rscript

#************************************
# FUNCTIONS
#
# DESCRIPTION: Functions scripts used to maked plots called by other scripts 
# in project 
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(ggfortify)

# Plotting qqplot (log10)
#-------------------------------------------------------------------------------
# Makes custumizable qqplot with log10 transformation (Made by Line Skotte)
qqp <- function(x, ci=TRUE, add=FALSE, ylab="Observed log10(p-value)", xlab="Expected log10(p-value)", maxLogP, ...){
  
  # Prep
  x<-x[!is.na(x)]
  if(!missing(maxLogP))
    x[x<10^-maxLogP]<-10^-maxLogP
  N<-length(x)
  x<-sort(x)
  
  # Plotting
  e <- -log((1:N-0.5)/N,10)
  if(add)
    points(e,-log(x,10),...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
    abline(0,1,col='firebrick',lwd=2)
  }
  
  # Add confidence intervals
  if(ci){
    c95<-qbeta(0.95,1:N,N-(1:N)+1)
    c05<-qbeta(0.05,1:N,N-(1:N)+1)
    lines(e,-log(c95,10))
    lines(e,-log(c05,10))
  }
}



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
  #pca_matrix_norm <- pca_matrix/(rowSums(pca_matrix)/1e6)
  
  # Perform the PCA
  pca <- prcomp(pca_matrix)
  
  # Make plot
  gene_pca_plot <- autoplot(pca,
                            data = design,
                            colour = group)
  
  # Extract data behind plot
  variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  
  # Extract PCA coordinates
  pca_source_data <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    mutate(
      PC1_variance_percent = variance_explained[1],
      PC2_variance_percent = variance_explained[2]
    ) %>% 
    left_join(select(design, "sample", group), by = "sample")
  
  return(list(
    plot = gene_pca_plot, 
    data = pca_source_data))
  
}

# Log fold estimates across pregnancy with confidence intervals for one gene plot
#-------------------------------------------------------------------------------
plot_gene_prediction <- function(data, selected_gene, gene_color){
  
  tmp_data <- data %>% 
    filter(gene_name == selected_gene)
  
  plot_title <- tmp_data %>% pull(gene_name_p_value) %>% unique() 
  
  gene_plot <- tmp_data %>% 
    ggplot() +
    geom_segment(aes(x = week_number, xend = week_end, y = Estimate, yend = Estimate, color =  "Estimate"), linewidth=2) + 
    geom_segment(aes(x = week_number, xend = week_end, y = CI_lower, yend = CI_lower, color = "95% CI"), linewidth = 0.5) + 
    geom_segment(aes(x = week_number, xend = week_end, y = CI_upper, yend = CI_upper, color = "grey"), linewidth = 0.5) + 
    scale_color_manual(values = c("Estimate" = gene_color, "95% CI" = "grey")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(title = plot_title, x = "Pregnancy week", y = "Log(fold change) to PP level", color = "Color") +
    theme(plot.title = element_text(size = 20),
          axis.text = element_text(size = 15),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15),
          legend.position = "bottom")
  
  return(gene_plot)
  
} 








