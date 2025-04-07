library(pheatmap)
library(tidyverse)

# Figure 3B.
pheatmap(mat = var, annotation = colData, main = 'All Conditions\nTop 500 Variable Genes', 
         treeheight_row = 0, 
         # treeheight_col = 0,
         cluster_cols = T, cluster_rows = T, annotation_names_col = F, show_rownames = F, annotation_legend = T,
         annotation_colors = list(condition = c(hCP='#E69F00', hNAT='#56B4E9', hNP='#009E73')))
