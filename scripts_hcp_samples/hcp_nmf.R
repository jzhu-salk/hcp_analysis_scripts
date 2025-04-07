library(NMF)

hcp_nmf = nmf(x = hcp_var, rank = seq(2,7), method = 'brunet', rng = 1, nrun = 500, .options = 'P45')
# Figure S4A.
plot(hcp_nmf, what='all')


# Consensus
consensus_num = 3
consensus_num = consensus_num - 1 # off by 1 since started with k=2
hcp_basis = hcp_nmf[["fit"]][[consensus_num]]@fit@H %>% apply(.,2, which.max) %>%
  data.frame(initial_basis = .) %>%
  mutate(basis = case_when(
    initial_basis == 3 ~ '1',
    initial_basis == 2 ~ '2',
    initial_basis == 1 ~ '3'
  ))

hcp_consensus = predict(hcp_nmf[["fit"]][[consensus_num]], what = 'consensus') %>% 
  # data.frame(consensus = .) %>%
  data.frame(initial_consensus = .) %>%
  mutate(consensus = case_when(
    initial_consensus == 3 ~ '1',
    initial_consensus == 2 ~ '2',
    initial_consensus == 1 ~ '3'
  ))


hcp_silhouette = silhouette(hcp_nmf[['fit']][[consensus_num]], what='consensus')

hcp_consensus_annotations = data.frame(row.names = row.names(hcp_consensus),
                                       Basis = hcp_basis$basis,
                                       Consensus = hcp_consensus$consensus) %>%
  merge(data.frame(Silhouette = hcp_silhouette[,3]), by = 'row.names') %>%
  column_to_rownames('Row.names') %>%
  dplyr::select(colnames(.) %>% rev) %>%
  arrange(Basis)

# in order of basis 1, 2, 3
hcp_consensus_order = c('hCP24', 'hCP51', 'hCP36', 'hCP49', 'hCP7', 'hCP17', 'hCP8', 'hCP16', 'hCP27', 'hCP1', 'hCP28', 'hCP44', 
                        'hCP29', 'hCP32', 'hCP18', 'hCP26', 'hCP20', 'hCP3', 'hCP25',
                        'hCP30', 'hCP15', 'hCP48', 'hCP11', 'hCP35')

# paletteMartin for annotation colors
# Figure 4B.
hcp_consensus_heatmap = hcp_nmf[["consensus"]][[consensus_num]][hcp_consensus_order, hcp_consensus_order] %>%
  pheatmap(mat = . , annotation_col = hcp_consensus_annotations,
           main = paste0('hCP NMF Consensus: Rank = ', consensus_num+1),
           border_color = NA, treeheight_row = 0, treeheight_col = 0, annotation_names_col = T,
           cluster_rows = F, cluster_cols = F, 
           annotation_colors = list(Basis = c('1' = "#0072B2", '2' = "#D55E00", '3' = "#CC79A7"),
                                    Etiology = c('A'="#AAAA49", 'D'="#009292",'G'="#ff6db6",'I'="#490092", 'O'="#006ddb"),
                                    Consensus = c('1'="#b66dff", '2'="#6db6ff", '3' = "#920000")))
