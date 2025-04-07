hcp_var = hcp_dds_final %>% getVarianceStabilizedData %>% data.frame %>%
  arrange(as.matrix(.) %>% rowVars) %>%
  tail(500)

hcp_pca = hcp_var %>% t %>% prcomp(center=T, scale.=T)

hcp_percentVar = format(round(100*summary(hcp_pca)$importance[2,], digits=2), digits = 2)

# Figure 4A.
hcp_pca_plot = hcp_pca$x %>% merge(colData, by = 'row.names') %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color=condition), size = 8, shape = 1, stroke = 2) +
  scale_color_manual(values = '#E69F00') +
  geom_text_repel(aes(label=Row.names), nudge_y = -1, size = 5)  +
  labs(title = 'PCA - hCP Patient Samples',
       subtitle = 'Top 500 Variable Genes',
       x = paste0("PC1 - ", hcp_percentVar[1], "%"),
       y = paste0("PC2 - ", hcp_percentVar[2], "%"),
       color = element_blank()) +
  theme_minimal()  +
  theme(
    legend.position = 'none',
  )
plot(hcp_pca_plot)
