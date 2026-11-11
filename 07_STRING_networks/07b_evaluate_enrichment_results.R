#!/usr/bin/env Rscript

#************************************
# EVALUATE AND PLOTS STRING RESULTS
#
# DESCRIPTION: This script uses the results from the enrichment analysis 
# performed in STRING and makes visualisations to compare the modules from WGCNA.
# Note that the data files are obtained via the analysis page, download section
# from the STRING permalinks.
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(enrichplot)
library(scales)
library(patchwork)


# Call projects paths
source("PROJECT_PATH/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))

# Load list of significant genes
gene_list <- read.table(file = file.path(path_results, "05_log_linear_mixed_effects_piecewise_model", "all_pvalues_with_gene_name.tsv"), as.is = TRUE, header = TRUE)

# List of files with STRING output files for all modules (NB: Filename must be initiated with module name)
module_files <- list.files(path = file.path(path_data, "07_enrichment_analyses", "STRING_modules"), pattern='*.tsv', full.names = TRUE)

# STRING enrichment genes for all significant genes 
string_gsc <- read.table(file = file.path(path_data, "07_enrichment_analyses", "enrichment.all.tsv"), as.is = TRUE, header = TRUE, sep = "\t", comment.char = "")

# File with number of (protein annotated) genes in each module for GeneRatio calculation and for all pregnancy associated genes
n_genes_file <- read.table(file = file.path(path_data, "07_enrichment_analyses", "n_proteins_all_sig_genes.tsv"), header = FALSE, sep = " ", col.names = c("n", "module"))
n_in_modules_file <- read.table(file = file.path(path_data, "07_enrichment_analyses", "STRING_modules", "protein_annotations", "n_proteins_in_modules.tsv"), header = FALSE, sep = " ", col.names = c("n", "module"))


#### Barplot of 12 most enriched terms
#-------------------------------------------------------------------------------
# Create full dataset for facet wrap plotting
n_terms <- 12
select_terms <- c("GO Process", "GO Function", "GO Component", "Reactome", "KEGG")

top_data_all_sig <- string_gsc %>% 
  dplyr::select(X.category, term.description, observed.gene.count, false.discovery.rate) %>% 
  filter(X.category %in% select_terms) %>% 
  group_by(X.category) %>% 
  slice_min(order_by = false.discovery.rate, n = n_terms, with_ties = FALSE) %>% 
  ungroup()

  
top_terms_preg_genes <- top_data_all_sig %>% 
  group_by(X.category) %>% 
  mutate(term.description = fct_reorder(term.description, false.discovery.rate, .desc = TRUE)) %>% 
  ungroup() %>% 
  ggplot(aes(x = observed.gene.count, y = term.description, fill = false.discovery.rate)) + 
  geom_bar(stat = "identity", alpha = 0.8) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
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
  scale_y_discrete(labels = label_wrap(30)) +
  scale_fill_gradient(low="firebrick3", high="royalblue", name = "Adj.\nP-value") +
  facet_wrap(~X.category, scales = "free", ncol = 3) +
  labs(x = "Annotated Gene Count", y = "Description")

top_terms_preg_genes


# Dotplot of most enriched terms for each co-expression module
#-------------------------------------------------------------------------------
##### Data convertion function #####
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
    module_color <- gsub("_[^_]+$", "", basename(file_name))
    
    # Read STRING enrichment results file for module and filter for selected term
    STRING_file <- read.table(file = file_name, as.is = TRUE, header = TRUE, sep = "\t", comment.char = "") %>% dplyr::filter(X.category == term)
    
    n <- filter(n_in_modules_file, module == module_color) %>% pull(n) 
    
    # Append needed data for results df
    STRING_Cluster <- c(STRING_Cluster, rep(module_color, nrow(STRING_file)))
    STRING_ID <- c(STRING_ID, STRING_file$term.ID)
    STRING_Description <- c(STRING_Description, STRING_file$term.description)
    if (nrow(STRING_file) > 0){
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


##### Prepare data #####
# Preselect which terms to plot
selected_loop_terms <- c("GO Process", "GO Function", "GO Component", "KEGG", "Reactome")

# Prepare data for plotting
# Consistent module order
module_order <- unique(
  gsub("_[^_]+$", "", basename(module_files))
)

# Prepare enrichment data
module_enrichment_list <- lapply(selected_loop_terms, function(term) {
  obj <- prepare_STRING_cC_data(module_files, term)
  as.data.frame(obj)
})

names(module_enrichment_list) <- selected_loop_terms

# Preselect for fixed sizing
full_module_enrichment_data <- bind_rows(module_enrichment_list)

##### Create plot limits #####
# Common limits across all plots
max_count <- full_module_enrichment_data %>% 
  mutate(GeneRatio_numeric = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio))) %>% 
  pull(GeneRatio_numeric) %>% 
  max(na.rm = TRUE)

# Ensure actual values
filter_pvalue <- full_module_enrichment_data %>% 
  filter(is.finite(p.adjust) & p.adjust >= 0 & !is.na(p.adjust)) %>% 
  dplyr::select(p.adjust)

# Extract min and max
fdr_limits <- range(
  filter_pvalue,
  na.rm = TRUE)

##### Create x-axis labels #####
# Module labels including number of proteins
module_n <- n_in_modules_file %>%
  distinct(module, .keep_all = TRUE)

x_labels <- setNames(
  paste0(module_order, 
         "\n(", 
         n_in_modules_file$n[match(module_order, n_in_modules_file$module)],
        ")"),
  module_order)


##### Create individual plots #####
# Set variables
plots <- list()
plot_heights <- numeric(length(selected_loop_terms))
n_show = 5 # Number of most enriched terms to display in plot


# Loop through selected terms
for (i in seq_along(selected_loop_terms)) {
  
  # Define data for selected term
  term <- selected_loop_terms[i]
  df <- module_enrichment_list[[term]]
  
  # Identify the 10 most significant terms per module
  top_terms <- df %>%
    group_by(Cluster) %>%
    slice_min(
      order_by = p.adjust,
      n = n_show,
      with_ties = FALSE
    ) %>%
    ungroup() %>%
    distinct(ID)
  
  df_wrangled <- df %>%
    semi_join(top_terms, by = "ID") %>% 
    mutate(Cluster = factor(Cluster, levels = module_order),
           Term = paste(ID, Description, sep = "___"),
           GeneRatio_numeric = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
  
  # Add to plot data
  if (i == 1){
    plot_data_df <- df_wrangled %>% 
      #mutate(Term = term) %>% 
      dplyr::select(Term, Cluster, Description, GeneRatio_numeric, p.adjust) }
  
  plot_data_df <- df_wrangled %>% 
    #mutate(Term = term) %>% 
    dplyr::select(Term, Cluster, Description, GeneRatio_numeric, p.adjust) %>% 
    bind_rows(plot_data_df)
  
  # Order terms by significance
  term_order <- df_wrangled %>%
    group_by(Term) %>%
    dplyr::summarise(
      best_FDR = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(best_FDR) %>%
    pull(Term)
  
  df_wrangled$Term <- factor(
    df_wrangled$Term,
    levels = rev(term_order)
  )
  
  # Wrapped enrichment term labels
  y_labels <- setNames(
    df_wrangled$Description,
    as.character(df_wrangled$Term)
  )
  
  p <- ggplot(df_wrangled, aes(x = Cluster, y = Term)) +
    geom_point(
      aes(size = GeneRatio_numeric, fill = p.adjust), shape = 21, color = "black", stroke = 0.2) +
    scale_x_discrete(
      limits = module_order,
      labels = x_labels,
      drop = FALSE) +
    scale_y_discrete(
      labels = y_labels) +
      #labels = label_wrap(45)) +
    scale_size_continuous(
      name = "Gene count",
      limits = c(0, max_count),
      range = c(1, 6)) +
    scale_fill_viridis_c(
      name = "Adjusted p-value (FDR)",
      limits = fdr_limits,
      direction = -1,
      option = "plasma") +
    labs(title = term,
         x = "Module (number of annotated genes in module)",
         y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 9),
          axis.title.x = element_text(size = 9, margin = margin(t = 2)),
          plot.title = element_text(size = 9, face = "bold", margin = margin(b = 1, t = 1)),
          legend.position = "right",
          plot.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "pt"))
  
  # Only display x-axis labels on bottom plot
  if (i < length(selected_loop_terms)) {
    p <- p +
      labs(x = NULL) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  plots[[i]] <- p
  
  # Calculate height based on wrapped labels
  unique_labels <- unique(y_labels)
  n_lines <- sum(str_count(unique_labels, "\n") + 1)
  plot_heights[i] <- max(2, n_lines * 0.22 + 0.7)
}

# Combine panels
combined_plot <- wrap_plots(
  plots,
  ncol = 1,
  heights = plot_heights,
  guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 12, face = "bold"),
        legend.position = "right")

combined_plot



# Output
#-------------------------------------------------------------------------------
# Save barplot with top enriched terms for pregnancy associated genes
png(file.path(path_results, "07_enrichment_analyses", "Barplot_top_terms_preg_genes.png"), width = 1800, height = 1000, units = 'px', res=120)
print(top_terms_preg_genes)
dev.off()

# Plot data
write.table(x = top_data_all_sig, file = file.path(path_results, "07_enrichment_analyses", "plot_tables", "top_term_barplot.tsv"), quote = FALSE, sep = "\t")

# Save module combined plot
ggsave(
  file.path(path_results, "07_enrichment_analyses", "STRING_enrichment_module_dotplots_combined.pdf"),
  combined_plot,
  width = 12,
  height = sum(plot_heights) * 0.50,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE)

# Plot data
write.table(x = plot_data_df, file = file.path(path_results, "07_enrichment_analyses", "plot_tables", "top_term_module_dotplots.tsv"), quote = FALSE, sep = "\t")


# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "07_enrichment_analyses", "session_info", "07b_session_info.txt"))






