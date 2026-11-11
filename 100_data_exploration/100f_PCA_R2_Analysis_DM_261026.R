library(tidyverse)
library(data.table)
library(DESeq2)
library(ggplot2)
library(pheatmap)

meta_data <- fread("PROJECT_PATH/data/04_normalization/design_finished.tsv")
head(meta_data)

path="PROJECT_PATH/data/03_gene_filtration/"
b1 = fread(paste0(path, "b1_filtered_counts.tsv"))
b2 = fread(paste0(path, "b2_filtered_counts.tsv"))
b1 <- filter(b1, V1 %in% intersect(pull(b1, V1), pull(b2, V1)))
b2 <- filter(b2, V1 %in% intersect(pull(b1, V1), pull(b2, V1)))

b1 <- tibble::column_to_rownames(as.data.frame(b1), var = "V1")
b2 <- tibble::column_to_rownames(as.data.frame(b2), var = "V1")

vst_b1 <- vst(as.matrix(b1))
vst_b2 <- vst(as.matrix(b2))

shared <- intersect(rownames(vst_b1), rownames(vst_b2))
joined <- cbind(
  vst_b1[shared, , drop = FALSE],
  vst_b2[shared, , drop = FALSE]
)
head(joined[,1:10])
dim(joined)

# reorder metadata to match joined columns
meta <- meta_data[match(colnames(joined), meta_data$sample), ]

# final alignment check
stopifnot(all(meta$sample == colnames(joined)))

# Batch correct:
design <- model.matrix(~ ga_adj1, data = meta)
joined_bc <- limma::removeBatchEffect(
  joined,
  batch = meta$batch,
  covariates = cbind(
    norm_factor = meta$norm_factor),
  design = design
)


### Make a meta data for the factors you want to see with correct encoding:
colnames(meta_data)
Sequencing_QC <- c("norm_factor", "lib_size", "N_multimapping", "N_noFeature", "N_total", "N_feature", 
                   "N_ribo",  "N_mapped", "library")
Biological_Char <- c("id", "batch", "postdelivery", "ga_adj1", "days_of_gestation", "natural_initiation", "mom_age_birth", "pre_preg_bmi", "parity", "birth_weight_g", "birth_length_cm", "sex_child")
seq_meta <- select(meta_data, all_of(c('sample', Sequencing_QC))) %>% mutate(library = factor(library))
bio_meta <- select(meta_data, all_of(c('sample', Biological_Char))) %>% 
mutate(gestational_age = ga_adj1, natural_initiation = as.logical(natural_initiation), fetal_sex = ifelse(sex_child == "M", T, F),
       batch = factor(batch), gestational_age = ifelse(postdelivery == TRUE, NA, gestational_age)) %>% select(-c(ga_adj1, sex_child))


# function to get PC-meta R2
pc_meta_r2 <- function(countdf, meta, n_pcs = 10, ngenes = 20000, filter_outliers = FALSE, outlier_mad = 5, plot = TRUE, title = NULL) {
  
  # countdf: features x samples
  # meta: metadata with UID column matching colnames(countdf)
  meta <- meta %>%  as.data.frame() %>% filter(sample %in% colnames(countdf)) %>% arrange(match(sample, colnames(countdf)))
  stopifnot(all(meta$sample == colnames(countdf)))
    
  # PCA on samples
  countdf <- countdf[order(apply(countdf, 1, var, na.rm = TRUE), decreasing = TRUE)[seq_len(ngenes)], ]
  message("Using ", nrow(countdf), " most variable genes")
  pcs <- prcomp(t(countdf), scale. = FALSE)
  pc_df <- as.data.frame(pcs$x[, seq_len(n_pcs), drop = FALSE])
  pc_df$sample <- rownames(pc_df)
  df <- left_join(pc_df, meta, by = "sample")
  pc_cols <- paste0("PC", seq_len(n_pcs))
  meta_cols <- setdiff(colnames(meta), c("sample", "UID"))
  
  res <- expand.grid(PC = pc_cols, metadata = meta_cols, stringsAsFactors = FALSE) %>%
    pmap_dfr(function(PC, metadata) {
        d <- df %>% select(all_of(c(PC, metadata))) %>% drop_na()
        # Remove extreme values of this PC
        if (filter_outliers) { 
            pc_med <- median(d[[PC]], na.rm = TRUE)
            pc_mad <- mad(d[[PC]], na.rm = TRUE)
            if (is.finite(pc_mad) && pc_mad > 0) { d <- d %>% filter(abs(.data[[PC]] - pc_med) / pc_mad <= outlier_mad) }
        }
        if (is.logical(d[[metadata]])) d[[metadata]] <- factor(d[[metadata]])
        if (nrow(d) < 3 || n_distinct(d[[metadata]]) < 2) {return(NULL)}
        #print(head(d))
        fit <- tryCatch(lm(reformulate(metadata, response = PC), data = d), error = function(e) NUL)
        if (is.null(fit)) return(NULL)
        tibble(PC = PC, metadata = metadata, R2 = summary(fit)$r.squared, pval = anova(fit)$`Pr(>F)`[1])
    }) %>% arrange(desc(R2))
  
  if (plot) {
    mat <- res %>% select(metadata, PC, R2) %>% pivot_wider(names_from = PC, values_from = R2, values_fill = 0) %>% as.data.frame()
    rownames(mat) <- mat$metadata
    mat$metadata <- NULL
    mat <- mat[, pc_cols, drop = FALSE]
    mat <- mat[apply(mat, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
    pheatmap(mat, cluster_cols = FALSE, cluster_rows = TRUE, display_numbers = TRUE, number_format = "%.2f", main = title %||% "Metadata variance explained in PCs")
  }
  return(list(res, pc_df))
}

pc_b1_bio <- pc_meta_r2(vst_b1, bio_meta[match(colnames(vst_b1), bio_meta$sample), ], title = 'Replication_Biological')
pc_b2_bio <- pc_meta_r2(vst_b2, bio_meta[match(colnames(vst_b2), bio_meta$sample), ], title = 'Discovery_Biological')
pc_b1_seq <- pc_meta_r2(vst_b1, seq_meta[match(colnames(vst_b1), seq_meta$sample), ], title = 'Replication_Technical')
pc_b2_seq <- pc_meta_r2(vst_b2, seq_meta[match(colnames(vst_b2), seq_meta$sample), ], title = 'Discovery_Technical')
pc_bio <- pc_meta_r2(joined, bio_meta[match(colnames(joined), bio_meta$sample), ], title = 'Combined_Biological')
pc_seq <- pc_meta_r2(joined, seq_meta[match(colnames(joined), seq_meta$sample), ], title = 'Combined_Technical')
pc_bio_bc <- pc_meta_r2(joined_bc, bio_meta[match(colnames(joined_bc), bio_meta$sample), ], title = 'Combined_Batch_Corrected_Biological')
pc_seq_bc <- pc_meta_r2(joined_bc, seq_meta[match(colnames(joined_bc), seq_meta$sample), ], title = 'Combined_Batch_Corrected_Technical')

# set order for covariates based on the joined data: 
plevels = bind_rows(pc_bio[[1]], pc_seq[[1]]) %>% mutate(PC = factor(PC, levels = paste0("PC", 1:10))) %>% arrange(PC, desc(R2)) %>% pull(metadata) %>% unique() %>% rev()
plevels <- gsub("_", " ", plevels)

# Visualisation function:
visualise_heatmap <- function(res, name) {
    df <- res[[1]] %>% arrange(metadata) %>% mutate(PC = factor(PC, levels = paste0("PC", 1:10)), 
                                                    meta_info = factor(gsub("_", " ", metadata), levels = plevels), 
                                                    padj = p.adjust(pval, method = 'BH'))
    name_clean <- gsub("_", ' ', name)
    p <- ggplot(df, aes(x = PC, y = meta_info, fill = R2, label = round(R2, 2))) +
          geom_tile() +
          geom_text(aes(label = round(R2, 2), fontface = ifelse(padj < 0.05, "bold", "plain"))) +
          colorspace::scale_fill_continuous_sequential(
            palette = "Reds",
            limits = c(0, 1)) +
          theme_bw()+
          ylab("Covariate")+
          ggtitle(name_clean)
    print(p)

    name <- paste0("PC_Heatmap_OR_", name , "_261026.pdf")
    pdf(name, width = 6, height = 6)
    print(p)
    dev.off()
}

visualise_heatmap(pc_b1_bio, "Replication_Biological")
visualise_heatmap(pc_b2_bio, "Discovery_Biological")
visualise_heatmap(pc_bio, "Combined_Biological")
visualise_heatmap(pc_bio_bc, "Combined_Batch_Corrected_Biological")

saveRDS(list('Discovery_Biological' = pc_b2_bio,
             "Replication_Biological" = pc_b1_bio, 
             "Combined_Biological" = pc_bio,
             "Combined_Batch_Corrected_Biological" = pc_bio_bc), 
             "PC_R2_Heatmap_OR_Plot_Biological_RawData_261026.RDS")           

openxlsx::write.xlsx(
  list(
    "Discovery" = pc_b2_bio[[1]],
    "Replication" = pc_b1_bio[[1]],
    "Combined" = pc_bio[[1]],
    "Combined_Batch_Corrected" = pc_bio_bc[[1]]
  ),
  "PC_R2_Biological_Heatmap_OR_Data_261026.xlsx")  

visualise_heatmap(pc_b1_seq, "Replication_Technical")
visualise_heatmap(pc_b2_seq, "Discovery_Technical")
visualise_heatmap(pc_seq, "Combined_Technical")
visualise_heatmap(pc_seq_bc, "Combined_Batch_Corrected_Technical")
saveRDS(list('RNASeq_Characterisation_Discovery_Technical' = pc_b2_seq, 
             "RNASeq_Characterisation_Replication_Technical" = pc_b1_seq, 
             "Combined_Technical" = pc_seq,
             "Combined_Batch_Corrected_Technical" = pc_seq_bc), "PC_R2_Heatmap_OR_Plot_Technical_RawData_261026.RDS")

openxlsx::write.xlsx(
  list(
    "Discovery" = pc_b2_seq[[1]],
    "Replication" = pc_b1_seq[[1]],
    "Combined" = pc_seq[[1]],
    "Combined_Batch_Corrected" = pc_seq_bc[[1]]
  ),
  "PC_R2_Technical_Heatmap_OR_Data_261026.xlsx"
)
