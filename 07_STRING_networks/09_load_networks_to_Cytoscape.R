#!/usr/bin/env Rscript

#************************************
# NETWORK VISUALIZATIONS
#
# Create network visualizations for enriched genes in modules
#
# Author: aidn@ssi.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
#  Initialize
# Make sure the required packages are installed. If any package is missing you 
# can always install it with the install.packages("packageName") command and load it with library(packageName).


if(!"tidyverse" %in% installed.packages()){
  install.packages('tidyverse')
}
if(!"dplyr" %in% installed.packages()){
  install.packages('dplyr')
}
if(!"stringr" %in% installed.packages()){
  install.packages("stringr")
}
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}
if(!"RColorBrewer" %in% installed.packages()){
  install.packages("RColorBrewer")
}

# Libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(RCy3)
library(RColorBrewer)
help(package=RCy3)

# Creat queries and networks
#-------------------------------------------------------------------------------
#  ------------------------------------------------------------------------------
#' @title String Protein Query
#'
#' @description Retrieves a string Network for a set of Proteins
#' @details Enter protein names or identifiers to query the STRING database for protein-protein interactions.
#' @param cutoff (optional) <BoundedDouble (0.0<=0.4<=1.0)>: The confidence score reflects the cumulated evidence that this interaction exists. Only interactions with scores greater than this cutoff will be returned
#' @param includesViruses (optional) true|false: By default, a query will search for identifiers in both the protein and virus databases. By changing this to 'false', only the protein database will besearched
#' @param limit (optional) <BoundedInteger (0<=10<=10000)>: The maximum number of proteins to return in addition to the query set
#' @param network (optional) network=current Specifies a network by name, or by SUID 
#' @param newNetName (required) newNetName=<String>: Name for the network to be created
#' @param query (required) <String>: Comma separated list of protein names or identifiers
#' @param species (optional) <String>: Species name. This should be the actual taxonomic name (e.g. homo sapiens, not human)
#' @param taxonID (optional) <int>: The species taxonomy ID. See the NCBI taxonomy home page for IDs. If both species and taxonID are set to a different species, the taxonID has priority.
#' @return None
#' @examples \dontrun{
#' 
#' }
#' @export
stringProteinQuery <-  
  function(cutoff=0.7, #change cutoff based on needs
           includesViruses="true",
           limit=0,
           network=NULL,
           newNetName=NULL,
           query=NULL,
           species="Homo sapiens",
           taxonID=9606) {
    if(is.null(query)){
      write("Please provide protein identifiers for the query",stderr())
    } else {
      cmd <- 
        paste0('string protein query cutoff=',
               cutoff,
               ' includesViruses=',
               includesViruses,
               ' limit=',
               limit,
               " network=",
               network,
               ' newNetName="',
               newNetName,
               '" query="',
               query,
               '" species="',
               species,
               '" taxonID=',
               taxonID
        )
      res <- commandsPOST(cmd)
      invisible(res)
    }
  }

# Check connection and close previous session if necessary
# Make sure you have the latest version of Cytoscape installed and you have already started Cytoscape. 
# If you want to close previous session uncomment the second line in the next block.

cytoscapePing()
# closeSession(TRUE)

# Path
path <- 'path/to/data'
#directory with Estimates per week
subpath1 <- '05_log_linear_mixed_effects_piecewise_model/'

subpath2 <- '06_WGCNA/'

subpath3 <- '09_STRING_networks/'

model_coef_b2 <- read.table(file = paste0(path, subpath1, "model_coef_b2.tsv"), as.is = TRUE, header = TRUE, sep = "\t")

# Model pvalues
pvals_df <- read.table(file = paste0(path, subpath1, "significant_genes_table.tsv"), as.is = TRUE, header = TRUE, row.names = 1)

pvals_df_no_version <- pvals_df %>%
  mutate(gene = str_remove(gene, "\\..*"))%>%
  select(gene, gene_name, pvalue_combined) %>%
  mutate(pvalue_combined = as.numeric(pvalue_combined))

write.csv(pvals_df_no_version, file = paste0(path, subpath3, "pvals_df_no_version.csv"))

max_p_value <- pvals_df_no_version %>%
  select(pvalue_combined) %>%
  max()

min_p_value <- pvals_df_no_version %>%
  select(pvalue_combined) %>%
  min()

# Calculate the mean of pvalue_combined
mean_p_value <- pvals_df_no_version %>%
  summarise(mean_pvalue_combined = mean(pvalue_combined)) %>%
  pull(mean_pvalue_combined)


median_p_value <- pvals_df_no_version %>%
  summarise(median_pvalue_combined = median(pvalue_combined)) %>%
  pull(median_pvalue_combined)


# WGCNA module file
b1_model_estimates_modules <- read.table(file = paste0(path, subpath2, "VST_model_normalization_signed_updated/gene_module_model_estimates.tsv"), as.is = TRUE, header = TRUE, sep = "\t") %>%
  select(gene, gene_name, color) %>%
  unique() 

# filter coefficient list down to 720 genes only.
# Use semi_join to filter model_coef_b2 based on the gene column from b1_model_estimates_modules
filtered_model_coef_b2 <- semi_join(model_coef_b2, b1_model_estimates_modules, by = "gene")

coef_effects_to_remove <- c("(Intercept)", "log(ribo_fraction)", "norm_factor", "log(lib_size)")

filtered_model_coef_b2_no_version <- filtered_model_coef_b2 %>%
  mutate(gene = str_remove(gene, "\\..*"))%>%
  mutate(coef_effect = str_remove(coef_effect, "ga_week_adjw")) %>%
  filter(!coef_effect %in% coef_effects_to_remove)

write.csv(filtered_model_coef_b2_no_version, file = paste0(path, subpath3, "filteredmodel_coef_b2_no_version.csv"))

### Make sure you have the required Apps installed inside Cytoscape
# Remove anything after the period in the gene names
b1_model_estimates_modules_no_version <- b1_model_estimates_modules %>%
  mutate(gene = str_remove(gene, "\\..*"))


# Group the genes by color and create a list of lists
genes_by_color <- b1_model_estimates_modules_no_version %>%
  group_by(color) %>%
  summarise(genes = list(gene)) %>%
  deframe()

# Display the list of lists
print(genes_by_color)


## Create STRING networks for each module
network_names <- c()
string_network_names <- c()
# set flag for visual style
flag = 0
for (color in names(genes_by_color)) {
  # Access the list of genes for the current color
  genes_list <- genes_by_color[[color]]
  
  # Print or process the genes list
  print(paste("Color:", color))
  print(genes_list)

  # create a comma separated string variable to store the protein names for the query
  query_string <- paste(genes_list, collapse=",") 
  Sys.sleep(0.5)
  # run string protein query using the function to retrieve interactions between the proteins in the module
  # alternatively uncomment the 2 lines of code below this command and run using that -- But don't do both!!!
  stringProteinQuery(includesViruses="false", limit=0, newNetName = color, query=query_string)
  Sys.sleep(1)
  # string.cmd <- paste('string protein query cutoff=0.4 includesViruses=false limit=0 taxonID=9606 newNetName=\"',disease_name,'\" query=\"',query_string, '\"',sep='')
  # commandsRun(string.cmd)

  ## ## #  VISUAL STYLE CREATION ## ## # 
  if(flag==0){
    #  Create a style to apply to all networks
    # Sometimes the app creates one style for each network and not a generally applicable style (e.g. STRING style v1.5), 
    # so you have to make sure one visual style exists that can be applied to all networks
    copyVisualStyle(from.style="STRING - black", to.style="Pregnancy style")
    flag=flag+1
    Sys.sleep(0.5)
    # Remove string glass balls and fill colors
    deleteStyleMapping(style.name="Pregnancy style", visual.prop="NODE_CUSTOMGRAPHICS_1")
    Sys.sleep(0.5)
    deleteStyleMapping(style.name="Pregnancy style", visual.prop="NODE_CUSTOMGRAPHICS_2")
    Sys.sleep(0.5)
    deleteStyleMapping(style.name="Pregnancy style", visual.prop="NODE_CUSTOMGRAPHICS_3")
    Sys.sleep(0.5)
    deleteStyleMapping(style.name="Pregnancy style", visual.prop="NODE_FILL_COLOR")
    Sys.sleep(0.5)
    # set default node color to gray
    setNodeColorDefault('#C0C0C0', style.name ="Pregnancy style")
    Sys.sleep(0.5)
    lockNodeDimensions(TRUE, style.name = "Pregnancy style")
    Sys.sleep(0.5)
    # Set node size to 30
    setNodeSizeDefault(new.size=30, style.name = "Pregnancy style")
    Sys.sleep(0.5)
    setNodeLabelPositionDefault(new.nodeAnchor='NE', new.graphicAnchor='SW', new.justification='c',new.xOffset=0,new.yOffset=0, style.name = "Pregnancy style")
    Sys.sleep(0.5)
    updateStyleMapping('Pregnancy style',mapVisualProperty('node label','display name','p'))
    Sys.sleep(0.5)
    
    # Load Omics Visualizer table once
    # The same table that we used previously to load the network will also be used to visualize data with Omics Visualizer
    table_file <- paste0(path, subpath3, "filteredmodel_coef_b2_no_version.csv")
    ov.cmd1 <- paste0('ov load file="',table_file,'" newTableName="All_data" firstRowAsColumnNames=true startLoadRow=1')
    commandsRun(ov.cmd1)
    
  }
  
  
  #####  END OF VISUAL STYLE CREATION ##### 
  
  # Apply the visual style you created above
  setVisualStyle("Pregnancy style")
  Sys.sleep(1)
  
  data <- read.delim(paste0(path, subpath3, "pvals_df_no_version.csv"), stringsAsFactors = FALSE, sep=",")
  # the key column from the attribute file.
  # that will be used to map the data to the preloaded network
  loadTableData(data,data.key.column="gene",table.key.column="query term")
  updateStyleMapping('Pregnancy style',mapVisualProperty('node fill color', 
                                                         'pvalue_combined','c', 
                                                         c(min_p_value,median_p_value, max_p_value),
                                                         c('#252525','#969696','#F7FCF5')))
  
  # connect OV table to network
  ov.cmd2 <-'ov connect mappingColNet="query term" mappingColTable="gene"'
  commandsRun(ov.cmd2)
  Sys.sleep(1)
  ov.cmd5 <- 'ov viz apply inner continuous colorMin="#2166AC" colorMax="#B2182B" colorMid="#ffffff" attributes="Estimate" rangeMax=2.8 rangeMin=-2.8'
  commandsPOST(ov.cmd5)
  Sys.sleep(1)
}

### Save session

saveSession(paste(getwd(),"pregnancy_networks_jul2025_2.cys", sep = '/'))

updateStyleMapping('Pregnancy style',mapVisualProperty('node label background transparency', 
                                                       'pvalue_combined','c', 
                                                       c(min_p_value,median_p_value, max_p_value),
                                                       c('#252525','#969696','#F7FCF5')))
setVisualPropertyDefault(list(visualProperty = "NODE_LABEL_BACKGROUND_COLOR", value = '#F7FCF5'),'Pregnancy style')
setVisualPropertyDefault(list(visualProperty = "NODE_LABEL_BACKGROUND_TRANSPARENCY", value = '200'),'Pregnancy style')
setVisualPropertyDefault(list(visualProperty = "NODE_LABEL_BACKGROUND_SHAPE", value = 'ROUND_RECTANGLE'),'Pregnancy style')

### Apply Layout
# These steps can be taken as a first effort to layout the network, but generally if the layout is applied manually the visualization is better.
# Cose seems to work best amongst the automated layouts tested.

network.list <- getNetworkList()
for (i in network.list){
  if (grepl("STRING network",i)){
    setCurrentNetwork(i)
    # apply organic layout - set as default through the GUI
    layoutNetwork()
    # apply CoSE layout
   #layoutNetwork("cose")
  }
}

# Output
#-------------------------------------------------------------------------------

##Export figures for all networks
for (i in network.list){
  if (grepl("STRING network",i)){
    setCurrentNetwork(i)
  exportPNG(filename=paste0(path, subpath3, i, '.png'))
  }}

