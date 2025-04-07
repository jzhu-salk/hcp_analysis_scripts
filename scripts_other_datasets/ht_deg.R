library(org.Hs.eg.db)

ht_mod_mat = model.matrix(design(ht_dds_final), ht_colData)

tumor = colMeans(ht_mod_mat[ht_dds_final$tissue_type == 'Tumor',])
normal = colMeans(ht_mod_mat[ht_dds_final$tissue_type == 'Normal',])

ht_res = results(ht_dds_final, contrast = tumor - normal)
summary(ht_res, alpha = .05)

res_to_df <- function(single_res){
  df = single_res %>% data.frame %>%
    mutate(ENSG = gsub("\\..*", "", row.names(.)),
           Symbol = mapIds(org.Hs.eg.db, keys = ENSG, column="SYMBOL", keytype="ENSEMBL", multiVals="first"),
           Cutoff = dplyr::case_when(
             is.na(padj) | padj > 0.05  ~ 'ns',
             log2FoldChange >  1 ~ 'inc',
             log2FoldChange < -1 ~ 'dec',
             .default = 'nfc'))
  return(df)
}

ht_df = res_to_df(ht_res) 