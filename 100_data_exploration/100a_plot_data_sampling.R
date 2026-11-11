#!/usr/bin/env Rscript

#************************************
# DATA STATISTICS PLOTS
#
# DESCRIPTION: This script creates relevant data descriptive plots for the manuscript
#
# Author: Freja Dahl Hede, FREH@dksund.dk
#************************************

# Imports and initialization
#-------------------------------------------------------------------------------
# Libraries
library(tidyverse)
library(scales)
library(lubridate)

# Call projects paths
source("/home/projects/ssi_gen1/data/pregnancy/BSiG_RNAseq_final/scripts/paths.R")

# Functions
source(file.path(path_project, "scripts", "functions_wrangling.R"))

# Keys that link library and barcode to samples for both batches
# Design (keys with meta data)
#design <- read.table(file = file.path(path_data, "01_data_wrangling", "design_data.tsv"), as.is = TRUE, header = TRUE)
design <- read.table(file = file.path(path_data, "04_normalization", "design_finished.tsv"), as.is = TRUE, header = TRUE)


# Plots
#-------------------------------------------------------------------------------
# Create colors 
postdelivery_colors <- c("Pregnancy sample" = "cornflowerblue", "Post partum sample" = "darkolivegreen4", "Birth" = "red")
postdelivery_shapes <- c("Pregnancy sample" = 16, "Post partum sample" = 16, "Birth" = 17)

#days_of_gestation = as.numeric(birth_date - (date_due - 40*7)),
#ga = (40*7 - as.numeric(date_due - date)))

data_sampling_plot_batch <- design %>% 
    mutate(batch_name = case_when(batch == "1" ~ "Replication",
                                  batch == "2" ~ "Discovery")) %>%
    mutate(#id = fct_reorder(id, as.numeric(ymd(birth_date)), .desc = TRUE),
           id = fct_reorder(id, days_of_gestation, .desc = TRUE),
           postdelivery = case_when(postdelivery == TRUE ~ "Post partum sample",
                                    postdelivery == FALSE ~ "Pregnancy sample")) %>% 
    ggplot() +
    geom_point(aes(x = ga, y = id, color = postdelivery, shape = postdelivery), size = 2) + # Samples
    geom_line(aes(x = ga, y = id), color = "cornflowerblue", linewidth = 0.2) + # Samples
    geom_point(aes(x = days_of_gestation, y = id, shape = "Birth", color = "Birth"), size = 2, position = position_nudge(y = -0.2)) + # Birth
    scale_x_continuous(breaks = seq(0, 280, by = 7*5), labels = seq(0, 40, by = 5)) +

    scale_color_manual(values = postdelivery_colors, name = "Samples") +  # Set custom colors
    scale_shape_manual(values = postdelivery_shapes, name = "Samples") +  # Set custom shapes
    theme_light() + 
    theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        legend.text = element_text(size=12), 
        legend.position = "bottom") +
    xlab("Gestational age (weeks)") + ylab("Participants") +
    #facet_wrap(~batch_name, nrow = 2, scales = "free_y", shrink = TRUE)  
    facet_grid(batch_name ~ ., scales="free", space = "free", switch="y") 
    

# Output
#-------------------------------------------------------------------------------
# Save plot
png(file.path(path_results, "100_data_exploration", "figures", "data_sampling_plot_split_by_batch.png"), width = 800, height = 1100, units = "px", res=120)
print(data_sampling_plot_batch)
dev.off()


# Session info
#-------------------------------------------------------------------------------
save_session_info(file.path(path_data, "100_data_exploration", "session_info", "100a_session_info.txt"))


