library(ggpubr)

# Also need to locally load ssGSEA-gpmodules for functions
# Downloaded from: https://github.com/GSEA-MSigDB/ssGSEA-gpmodule
# source('ssGSEA-gpmodule/src/common.R')
# source('ssGSEA-gpmodule/src/ssGSEA.Library.R')

# Run ssGSEA ----
input_gct = 'gsea/all_samples_gsea_norm_counts.gct'
ductal_acinar_gene_set = 'gsea/acinar_ductal/gene_set/acinar_ductal_gene_set.gmt'

ssGSEA.project.dataset(input.ds = input_gct,
                       output.ds = 'gsea/acinar_ductal/ssgsea_output/acinar_ductal_ssgsea',
                       gene.sets.dbfile.list = ductal_acinar_gene_set,
                       gene.symbol.column = 'Name', gene.set.selection = 'ALL',
                       sample.norm.type = 'none', weight = 0.75, combine.mode = 'combine.add', min.overlap = 1)

# Plot ssGSEA Results ----
scores = read.csv('gsea/acinar_ductal/ssgsea_output/acinar_ductal_ssgsea.gct', sep='\t', skip = 2) %>%
  dplyr::select(-Description) %>%
  column_to_rownames('Name')

nmf_clusters = hcp_consensus_annotations %>% dplyr::select(Basis) %>%
  merge(colData, by = 'row.names', all=T) %>%
  rename(Condition = 'condition') %>%
  column_to_rownames('Row.names') %>%
  arrange(Basis)

# Figure S3B.
pheatmap(mat = scores, cluster_rows = T, cluster_cols = T,  scale = 'row', 
         annotation_col = nmf_clusters, 
         annotation_colors = list(Basis = c('1' = "#0072B2", '2' = "#D55E00", '3' = "#CC79A7"),
                                  Condition = c(hCP='#E69F00', hNAT='#56B4E9', hNP='#009E73')), 
         treeheight_row = 25, treeheight_col = 25,
         main = 'Ductal vs Acinar ssGSEA')
