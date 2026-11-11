# To obtain gene symbol file for ENSG --> gene symbol conversion for GENCODE release 28

curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.metadata.HGNC.gz > ../../data/05_log_linear_mixed_effects_piecewise_model/gencode.v28.metadata.HGNC.gz
gunzip ../../data/05_log_linear_mixed_effects_piecewise_model/gencode.v28.metadata.HGNC.gz


(Note, the remaining scripts are for additional data accession to make decisions on models, thresholds etc., but are not needed in the downstream analyses)
(Note, in 05b_evaluate_model_results.R there is a section making a log fold change plot for three genes where they are found to be significant in both the discovery and replication batch, but where the p value of the models made on the combined batches is not lower. This was checked to decide that we will utilize a conservative threshold/criteria for when a gene is considered significant)




