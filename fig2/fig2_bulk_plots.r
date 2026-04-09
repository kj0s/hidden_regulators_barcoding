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
