#!/usr/bin/env Rscript

#************************************
# EVALUATE AND PLOTS STRING RESULTS
#
# DESCRIPTION: This script uses the results from the enrichment analysis 
# performed in STRING and makes visualisations to compare the modules from WGCNA.
# Note that the data files are obtained via the analysis page, download section
# from the STRING permalinks.
#
# Author: Freja Dahl Hede, FREH@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Package inst.
if (!("vissE" %in% installed.packages())) { # required to make one of the clusterProfiler plots
  # Install this package if it isn't installed yet
  BiocManager::install("vissE")
}

# Libraru
library(tidyverse)
library(vissE)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(GOSemSim)
library(ggpubr)
library(scales)

# Path
path <- '/PATH/data/'

# Load list of significant genes
gene_list <- read.table(file = paste0(path, "07_STRING_networks/significant_gene_list.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

# List of files with STRING output files for all modules (NB: Filename must be initiated with module name)
module_files <- list.files(path=paste0(path, "07_STRING_networks/STRING_modules"), pattern='*.tsv', full.names = TRUE)

# File with number of (protein annotated) genes in each module for GeneRatio calculation
n_in_modules_file <- read.table(file = paste0(path, "07_STRING_networks/STRING_modules/protein_annotations/", "n_proteins_in_modules.tsv"), header = FALSE, sep = " ", col.names = c("n", "module"))

# STRING enrichment genes for all significant genes 
string_gsc <- read.table(file = paste0(path, "07_STRING_networks/STRING_modules/all_sig_genes/enrichment.all.tsv"), as.is = TRUE, header = TRUE, sep = "\t", comment.char = "")

# File with number of (protein annotated) genes in each module for GeneRatio calculation
n_genes_file <- read.table(file = paste0(path, "07_STRING_networks/STRING_modules/all_sig_genes/protein_annotations/", "n_proteins_in_modules.tsv"), header = FALSE, sep = " ", col.names = c("n", "module"))

# Data convertion (STRING --> clusterProfiler) for all significant genes
#-------------------------------------------------------------------------------
# Calculate number of genes overlapping with entire geneset (GO Biological Process in below example)
n <- dplyr::select(tibble(string_gsc), matching.proteins.in.your.network..labels.) %>% 
  separate_longer_delim(matching.proteins.in.your.network..labels., delim = ",") %>% 
  unique() %>% 
  nrow()

# Create result df for enrichGO object
results_df <- data.frame("ID" =  string_gsc$term.ID,
                         "Description" = string_gsc$term.description,
                         "GeneRatio" = paste0(string_gsc$observed.gene.count, "/", n), # To avoid "/0" 
                         "BgRatio" = NA, # unsure what N and possibly also M is
                         "RichFactor" = as.numeric(NA), # = count / N, 
                         "FoldEnrichment" = as.numeric(NA),
                         "zScore" = as.numeric(NA), 
                         "pvalue" = as.numeric(NA),
                         "p.adjust" = as.numeric(string_gsc$false.discovery.rate) ,
                         "qvalue" = as.numeric(NA), #qscore = -log(p.adjust, base=10), https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
                         "geneID" = string_gsc$matching.proteins.in.your.network..labels.,
                         "Count" = as.integer(string_gsc$observed.gene.count))

# Create enrichResult object 
eR_object <- new("enrichResult",
                 result = results_df,
                 #pvalueCutoff = 0.05, # dotplot doesn't work with this
                 pAdjustMethod = "BH",
                 organism = "Homo Sapiens",
                 ontology = "BP",
                 gene = pull(gene_list, gene_name),
                 keytype = "SYMBOL")   

# FUNCTIONS 
#-------------------------------------------------------------------------------
# Make function that load all STRING output file, extracts data for selected 
# gene set (eg. GO BP), and converts to enrichResult input format
prepare_STRING_eR_data <- function(STRING_file, term, n){
  
  # Make results dataframe
  results_df_sig_genes <- data.frame("ID" =  STRING_file_term$term.ID,
                                     "Description" = STRING_file_term$term.description,
                                     "GeneRatio" = paste0(STRING_file_term$observed.gene.count, "/", n),
                                     "BgRatio" = NA, # unsure what N and possibly also M is
                                     "RichFactor" = as.numeric(NA), # = count / N, 
                                     "FoldEnrichment" = as.numeric(NA),
                                     "zScore" = as.numeric(NA), 
                                     "pvalue" = as.numeric(NA),
                                     "p.adjust" = STRING_file_term$false.discovery.rate,
                                     "qvalue" = as.numeric(NA), #qscore = -log(p.adjust, base=10), https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
                                     "geneID" = STRING_file_term$matching.proteins.in.your.network..labels.,
                                     "Count" = STRING_file_term$observed.gene.count)
  
  # Create compareCluster results object
  eR_object_sig_genes <- new("enrichResult",
                             result = results_df_sig_genes) 
  
  return(eR_object_sig_genes)
}  

# Make function that load all STRING output files, extracts data for selected 
# gene set (eg. GO BP), and converts to compareCluster input format
prepare_STRING_cC_data <- function(STRING_files_list, term, n_file){
  # Make empty vector for loop filling of data
  STRING_Cluster <- c()
  STRING_ID <- c()
  STRING_Description <- c()
  STRING_GeneRatio <- c()
  STRING_p.adjust <- c()
  STRING_geneID <- c()
  STRING_Count <- c()
  
  # Loop through files and extract relevant data
  for (file_name in STRING_files_list){
    # Extract module color
    module_color <- gsub(".*/|\\.[^.]+", "", file_name)
    
    # Read STRING enrichment results file for module and filter for selected term
    STRING_file <- read.table(file = file_name, as.is = TRUE, header = TRUE, sep = "\t", comment.char = "") %>% dplyr::filter(X.category == term)
    
    # Calculate number of annotated proteins
    n <- filter(n_file, module == module_color) %>% pull(n) 
    
    # Append needed data for results df
    STRING_Cluster <- c(STRING_Cluster, rep(module_color, nrow(STRING_file)))
    STRING_ID <- c(STRING_ID, STRING_file$term.ID)
    STRING_Description <- c(STRING_Description, STRING_file$term.description)
    if (nrow(STRING_file > 0)){
      STRING_GeneRatio <- c(STRING_GeneRatio, paste0(STRING_file$observed.gene.count, "/", n))} # To avoid "/0" 
    STRING_p.adjust <- c(STRING_p.adjust,  as.numeric(STRING_file$false.discovery.rate))
    STRING_geneID <- c(STRING_geneID, STRING_file$matching.proteins.in.your.network..labels.)
    STRING_Count <- c(STRING_Count, as.integer(STRING_file$observed.gene.count))
    
  }
  
  # Make results dataframe
  results_df_modules <- data.frame("Cluster" = STRING_Cluster, 
                                   "ID" =  STRING_ID,
                                   "Description" = STRING_Description,
                                   "GeneRatio" = STRING_GeneRatio,
                                   "BgRatio" = NA, # unsure what N and possibly also M is
                                   "RichFactor" = as.numeric(NA), # = count / N, 
                                   "FoldEnrichment" = as.numeric(NA),
                                   "zScore" = as.numeric(NA), 
                                   "pvalue" = as.numeric(NA),
                                   "p.adjust" = STRING_p.adjust,
                                   "qvalue" = as.numeric(NA), #qscore = -log(p.adjust, base=10), https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
                                   "geneID" = STRING_geneID,
                                   "Count" = STRING_Count)
  
  # Create compareCluster results object
  cC_object_modules <- new("compareClusterResult",
                           compareClusterResult = results_df_modules) 
  
  return(cC_object_modules)
}

# Enrichment analysis: All significant genes, selected terms
#-------------------------------------------------------------------------------
# Create enrichResults object

# Extract n
n_sig_genes <- n_genes_file %>% pull(n)

# Letters
letter <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L")

# Make empty list for plotting
term_barplot <- list()

# Letter counter
l <- 1

# Create selected term bar plots
for (loop_term in c("GO Process", "GO Function", "GO Component", "Reactome", "KEGG")){
 # Create clusterProfiler object 
  eR_object_sig_genes_term <- prepare_STRING_eR_data(string_gsc, loop_term, n_sig_genes)
  
  print("Created eR object")
  # Make each plot in list for grid plotting
  term_barplot[[loop_term]] <- barplot(eR_object_sig_genes_term, showCategory = 12) + 
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          plot.title = element_text(size = 18),
          legend.position = "right", 
          legend.text = element_text(size=6), 
          legend.title = element_text(size=9),
          plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
    labs(title = letter[l], 
         subtitle = loop_term) +
    guides(fill=guide_legend(title="Adj. p-value")) + #, ncol = 1)) +
    scale_y_discrete(labels = label_wrap(40))
  
  l <- l + 1
  
}

do.call(ggarrange, c(term_barplot, ncol = 3, nrow = 2, align="hv"))#, common.legend = TRUE, legend="right"))


# Enrichment analysis: All significant genes, selected terms BUT common scale
#-------------------------------------------------------------------------------
# Create enrichResults object
n_terms <- 12

# Extract n
n_sig_genes <- n_genes_file %>% pull(n)

# Create full dataframe to join with
full_term_data <- head(prepare_STRING_eR_data(string_gsc, "GO Function", n_sig_genes)@result, n=0)

# Create full dataset for facet wrap plotting
for (loop_term in c("GO Process", "GO Function", "GO Component", "Reactome", "KEGG")){
  
  # Create clusterProfiler object with term column 
  eR_object_sig_genes_term <- prepare_STRING_eR_data(string_gsc, loop_term, n_sig_genes)@result 
  
  # Select only top enriched terms
  eR_object_sig_genes_term_wrangled <- eR_object_sig_genes_term %>% 
    arrange(p.adjust) %>% 
    head(n = n_terms) %>% 
    mutate(term = loop_term) 

  
  # Add term info to data
  full_term_data <- bind_rows(full_term_data, eR_object_sig_genes_term_wrangled)
  
}

total_plot <- full_term_data %>% 
  group_by(term) %>% 
  mutate(Description = fct_reorder(Description, p.adjust, .desc = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = Count, y = Description, fill = p.adjust)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title = element_text(size = 13),
        plot.title = element_text(size = 18),
        legend.position = "right", 
        legend.text = element_text(size=9), 
        legend.title = element_text(size=11),
        plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) +
  #labs(title = letter[l], 
  #     subtitle = loop_term) +
  #guides(fill=guide_legend(title="Adj. p-value")) + #, ncol = 1)) +
  scale_y_discrete(labels = label_wrap(43)) +
  scale_fill_gradient(low="firebrick3", high="royalblue", name = "Adj.\nP-value") +
  facet_wrap(~term, scales = "free", ncol = 3)

total_plot

# Enrichment analysis: All significant genes, all terms
#-------------------------------------------------------------------------------
# Make empty list for plotting
term_barplot <- list()

# Letter counter
l <- 1
# Create selected term bar plots
for (loop_term in c("GO Process", "GO Function", "GO Component", "COMPARTMENTS", "KEGG", "Reactome", "DISEASES", "TISSUES")){
  # Create clusterProfiler object 
  eR_object_sig_genes_term <- prepare_STRING_eR_data(string_gsc, loop_term, n_sig_genes)
  
  print("Created eR object")
  # Make each plot in list for grid plotting
  term_barplot[[loop_term]] <- barplot(eR_object_sig_genes_term, showCategory = 12) + 
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 9),
          plot.title = element_text(size = 18),
          legend.position = "right", 
          legend.text = element_text(size=6), 
          legend.title = element_text(size=9)) +
    labs(title = letter[l], 
         subtitle = loop_term) +
    guides(fill=guide_legend(title="Adj. p-value")) + #, ncol = 1)) +
    scale_y_discrete(labels = label_wrap(40))
  
  l <- l + 1
  
}

do.call(ggarrange, c(term_barplot, ncol = 3, nrow = 3, align="hv"))#, common.legend = TRUE, legend="right"))

# Data convertion (STRING --> clusterProfiler) for module genes
#-------------------------------------------------------------------------------
# Select term for visualization
select_term <- "GO Process"

# Make function that load all STRING output files, extracts data for selected 
# gene set (eg. GO BP), and converts to compareCluster input format
prepare_STRING_cC_data <- function(STRING_files_list, term){
  # Make empty vector for loop filling of data
  STRING_Cluster <- c()
  STRING_ID <- c()
  STRING_Description <- c()
  STRING_GeneRatio <- c()
  STRING_p.adjust <- c()
  STRING_geneID <- c()
  STRING_Count <- c()
  
  # Loop through files and extract relevant data
  for (file_name in STRING_files_list){
    # Extract module color
    module_color <- gsub(".*/|\\.[^.]+", "", file_name)

    # Read STRING enrichment results file for module and filter for selected term
    STRING_file <- read.table(file = file_name, as.is = TRUE, header = TRUE, sep = "\t", comment.char = "") %>% dplyr::filter(X.category == term)
    
    n <- filter(n_in_modules_file, module == module_color) %>% pull(n) 

    # Append needed data for results df
    STRING_Cluster <- c(STRING_Cluster, rep(module_color, nrow(STRING_file)))
    STRING_ID <- c(STRING_ID, STRING_file$term.ID)
    STRING_Description <- c(STRING_Description, STRING_file$term.description)
    if (nrow(STRING_file > 0)){
    STRING_GeneRatio <- c(STRING_GeneRatio, paste0(STRING_file$observed.gene.count, "/", n))} # To avoid "/0" 
    STRING_p.adjust <- c(STRING_p.adjust,  as.numeric(STRING_file$false.discovery.rate))
    STRING_geneID <- c(STRING_geneID, STRING_file$matching.proteins.in.your.network..labels.)
    STRING_Count <- c(STRING_Count, as.integer(STRING_file$observed.gene.count))
    
  }
  
  # Make results dataframe
  results_df_modules <- data.frame("Cluster" = STRING_Cluster, 
                                   "ID" =  STRING_ID,
                                   "Description" = STRING_Description,
                                   "GeneRatio" = STRING_GeneRatio,
                                   "BgRatio" = NA, # unsure what N and possibly also M is
                                   "RichFactor" = as.numeric(NA), # = count / N, 
                                   "FoldEnrichment" = as.numeric(NA),
                                   "zScore" = as.numeric(NA), 
                                   "pvalue" = as.numeric(NA),
                                   "p.adjust" = STRING_p.adjust,
                                   "qvalue" = as.numeric(NA), #qscore = -log(p.adjust, base=10), https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
                                   "geneID" = STRING_geneID,
                                   "Count" = STRING_Count)
  
  
  
  # Create compareCluster results object
  cC_object_modules <- new("compareClusterResult",
                           compareClusterResult = results_df_modules) 
  
  return(cC_object_modules)
}

# Create object for GO BP
BP_cC_object_modules <-  prepare_STRING_cC_data(files, select_term)


# Plots
#-------------------------------------------------------------------------------
# Create all module dot plots
for (loop_term in c("GO Process", "GO Function", "GO Component", "COMPARTMENTS", "KEGG", "Reactome", "DISEASES", "TISSUES")){
  print(loop_term)
  term_cC_object_modules <-  prepare_STRING_cC_data(files, loop_term)
  
  module_dotplot <- dotplot(term_cC_object_modules) + 
    theme(axis.text.y = element_text(size = 8)) +
    labs(x = "Module (and number of annotated genes in module)", title = loop_term) +
    scale_y_discrete(labels = label_wrap(45))
  
}





