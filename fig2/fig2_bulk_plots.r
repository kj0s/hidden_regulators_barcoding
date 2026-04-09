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

obj <- readRDS("ST223.annotated_fate.rds")

df <- obj@meta.data %>%
  dplyr::select(fate_umap, rna_umap, fate_cluster, seurat_clusters) %>%
  na.omit()

ggplot(df, aes(x = fate_umap, y = rna_umap, color = as.factor(fate_cluster))) +
  geom_point(size = 1) +
  theme_classic() +
  labs(color = "Fate cluster")

ggplot(df, aes(x = fate_umap, y = rna_umap, color = as.factor(seurat_clusters))) +
  geom_point(size = 1) +
  theme_classic() +
  labs(color = "WNN cluster")

set.seed(18)
df$random <- runif(nrow(df))

ggplot(df, aes(x = fate_umap, y = random, color = as.factor(fate_cluster))) +
  geom_point(size = 0.8, alpha = 0.8) +
  theme_classic()

## check if broken: 
colnames(obj@meta.data)
df <- obj@meta.data %>%
  rename(
    fate_umap = YOUR_COLUMN,
    rna_umap = YOUR_OTHER_COLUMN
  )
