hcp_gene_mat = gene_mat %>%
  dplyr::select(grep('^hCP', colnames(.)))

hcp_colData = colnames(hcp_gene_mat) %>%
  data.frame(row.names = .)

hcp_dds = DESeqDataSetFromMatrix(countData = hcp_gene_mat,
                                 colData = hcp_colData,
                                 design = ~ 1) 

hcp_dds %>% counts %>% rowSums %>% summary
# remove genes with low reads
hcp_has_reads = rowSums(counts(hcp_dds) > 10 ) > 3

hcp_dds_filtered = hcp_dds[hcp_has_reads,]
hcp_dds_filtered %>% counts %>% rowSums %>% summary
hcp_dds_final = DESeq(hcp_dds_filtered)

# Prepare Normalized Counts
library(org.Hs.eg.db)
hcp_norm_counts = hcp_dds_final %>% counts(normalized=T) %>% data.frame %>%
  mutate(Symbol = mapIds(org.Hs.eg.db, keys = gsub('\\..*', '', rownames(.)), column="SYMBOL", keytype="ENSEMBL", multiVals="first"))
# write.table(hcp_norm_counts, file = 'norm_counts/hcp_norm_counts.tsv', quote = F, sep='\t')

hcp_gsea_counts = hcp_norm_counts %>%
  filter(!is.na(Symbol)) %>%
  remove_rownames %>%
  distinct(Symbol, .keep_all = T) %>%
  column_to_rownames('Symbol') %>%
  mutate(DESCRIPTION = NA) %>%
  relocate(DESCRIPTION) %>%
  rownames_to_column(var = 'NAME')
# write.table(hcp_gsea_counts, file = 'gsea/hcp_gsea_norm_counts.tsv', sep = '\t', row.names = F, quote = F)

# Will need to manually add additional lines to top of this file and rename as gct for GSEA:
# cat > hcp_norm_counts.gct
# #1.2
# 14419   24

# cat hcp_gsea_norm_counts.tsv >> hcp_gsea_norm_counts.gct