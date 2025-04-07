library(dplyr)
library(org.Hs.eg.db)

#######################################################
mod_mat = model.matrix(design(dds_final), colData(dds_final))
hNAT = colMeans(mod_mat[dds_final$condition == 'hNAT',])
hCP = colMeans(mod_mat[dds_final$condition == 'hCP',])
hNP = colMeans(mod_mat[dds_final$condition == 'hNP',])

##############################################
res_to_df <- function(single_res){
  df = single_res %>% data.frame %>%
    mutate(ENSG = gsub("\\..*", "", rownames(.)),
           Symbol = mapIds(org.Hs.eg.db, keys = ENSG, column="SYMBOL", keytype="ENSEMBL", multiVals="first"),
           Cutoff = dplyr::case_when(
             is.na(padj) | padj > 0.05  ~ 'ns',
             log2FoldChange >  1 ~ 'inc',
             log2FoldChange < -1 ~ 'dec',
             .default = 'nfc'))
  return(df)
}

# Individual v Rest Comparisons ----
hcp_res = results(dds_final, contrast = hCP - (hNP/2 + hNAT/2))
summary(hcp_res, alpha=.05)
hcp_df = res_to_df(hcp_res)

hnp_res = results(dds_final, contrast = hNP - (hCP/2 + hNAT/2))
summary(hnp_res, alpha=.05)
hnp_df = res_to_df(hnp_res)

hnat_res = results(dds_final, contrast = hNAT - (hNP/2 + hCP/2))
summary(hnat_res, alpha=.05)
hnat_df = res_to_df(hnat_res)
