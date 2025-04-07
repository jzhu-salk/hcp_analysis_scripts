library(ggplot2)
library(ggrepel)
library(colorBlindness)

# Figure 3A.
var = dds_final %>% getVarianceStabilizedData %>% data.frame %>%
  arrange(as.matrix(.) %>% rowVars) %>%
  tail(500)

pca = var %>% t %>% prcomp(center=T, scale.=T)

percentVar = round(100*summary(pca)$importance[2,], digits=2)

pca_plot = pca$x %>% merge(colData, by = 'row.names') %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(aes(color=condition), size = 8, shape = 1, stroke = 2) +
  scale_color_manual(values = palette.colors(n = 9)[2:4]) +
  geom_text_repel(aes(label=Row.names), nudge_y = -1)  +
  labs(title = 'PCA - All Patient Samples',
       subtitle = 'Top 500 Variable Genes',
       x = paste0("PC1 - ", percentVar[1], "%"),
       y = paste0("PC2 - ", percentVar[2], "%"),
       color = element_blank()) +
  theme_minimal()

plot(pca_plot)
