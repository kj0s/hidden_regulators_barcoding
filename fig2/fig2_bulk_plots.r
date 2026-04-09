# mapping clusters:

cluster_map <- data.frame(
  seurat_clusters = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"),
  celltype = c("primitive-RBC","Mk-ery_prog","Mk",
               "early_lymphoid_progenitor","HSC",
               "basophil-mast_progenitor","pre-cDC","Ery",
               "Mk","CDP","CD14+Monocyte","CMP-GMP",
               "Ery","Ery","Mk","unknown")
)

df <- left_join(df, cluster_map, by = "seurat_clusters")

library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("/vast/projects/Sisseq/ST223/ST223.annotated_fate.rds")

ST223.annotated_fate@meta.data$rna_umap <- ST223.annotated_fate@reductions$X1d.rna.umap@cell.embeddings[,1]

# plot
ggplot(ST223.annotated_fate@meta.data, aes(x = fate_umap, y = rna_umap, colour = as.factor(fate_cluster))) +
  geom_point(size = 1) +
  scale_colour_manual(values = cluster_colour) +
  labs(x = "Fate UMAP 1D", y = "WNN UMAP 1D", colour = "Fate Cluster")

## check if broken: 
# colnames(obj@meta.data)
# df <- obj@meta.data %>%
#  rename(
#    fate_umap = YOUR_COLUMN,
#    rna_umap = YOUR_OTHER_COLUMN
#  )
