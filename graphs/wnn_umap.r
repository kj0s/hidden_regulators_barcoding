# Load libraries
library(Seurat)
library(ggplot2)

# Load object
ST223.annotated <- readRDS("ST223.annotated.rds")

# Plot WNN UMAP
DimPlot(ST223.annotated, reduction = 'wnn.umap',
        label = TRUE, repel = TRUE, label.size = 4) + NoLegend()
