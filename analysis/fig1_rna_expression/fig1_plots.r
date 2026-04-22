# Final Run Check
# Wed Apr 22 2026
# KJ 
# needed files: ST227_Day3_21_GP.rds
# outputs: RNA_markers.txt

# Installing pckages
BiocManager::install("SingleR")
BiocManager::install("BiocManager")
install.packages("pak")
pak::pkg_install("immunogenomics/presto")

#finding reference dataset in celldex
library(celldex)
library(SingleR)
library(Seurat)
library(dplyr)
surveyReferences()

ST227 <- readRDS(".../ST227_Day3_21_GP.rds")

#Redo clustering and umap for day 3 and 21 only
ST227 <- FindMultiModalNeighbors(
  ST227, reduction.list = list("pca", "sketch_apca"), 
  dims.list = list(1:27, 1:19), modality.weight.name = list("RNA.weight", "ADT.weight")
)

ST227 <- RunUMAP(ST227, nn.name = "weighted.nn", n.neighbors = 30L, reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ST227 <- FindClusters(ST227, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

DefaultAssay(ST227) <- 'sketch_RNA'
RNA_markers <- FindAllMarkers(ST227, only.pos = TRUE)

DefaultAssay(ST227) <- 'sketch_ADT'
ADT_markers <- FindAllMarkers(ST227, only.pos = TRUE)
# note: No features pass logfc.threshold threshold; returning empty data.frame here. 

RNA_markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

write.table(top10, file= "RNA_markers.txt", sep="\t", row.names=FALSE, col.names=TRUE)
write.table(top10, file = "/vast/projects/Sisseq/previous_code_test_KJ/intermediates/RNA_markers.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


DoHeatmap(ST227, features = top10$gene) + NoLegend()
# ADT layer empty here. data not normalised?

#Remove contaminating stromal cells
ST227 <- subset(ST227, idents = c("19", "29"), invert = TRUE)

#Rename clusters
new.cluster.ids <- c("Ery", "cDC2", "Ery", "Ery", "CD14+ Mono", "T cells", "Pre-B cells", "Mast Cells", "Neutrophils", "Megakaryocytes", "Megakaryocytes", "DC3", "Eosinophils", "HSC", "NK cells", "Basophils", "Mast Cells", "cDC1", "CD16+ Mono", "mregDC", "GMP", "cDC2", "NK cells", "Macrophages", "T cells", "cDC2", "Megakaryocytes", "Eosinophils")

# Assign new names
ST227$celltype <- plyr::mapvalues(Idents(ST227), 
                                  from = levels(Idents(ST227)), 
                                  to = new.cluster.ids)
Idents(ST227) <- ST227$celltype

p1 <- DimPlot(ST227, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
ggplot2::ggsave("annotated_umap.pdf", plot = p1, width = 7, height = 6, units = "in") 
p1
dev.off()

#Compare with BM dataset
ST227.sce.fig <- as.SingleCellExperiment(ST227, assay = "sketch_RNA")

# automatically annotate using this prior dataset


ref <- fetchReference("novershtern_hematopoietic", "2024-02-26")
pred.main <- SingleR(test = ST227.sce.fig, ref = ref, labels = ref$label.main)

ST227[["SingleR.labels"]] <- pred.main$labels
Idents(ST227) <- "SingleR.labels"

plotScoreHeatmap(pred.main,
                 clusters = ST227$celltype, order.by.clusters = TRUE,
                 show.labels = TRUE, show.pruned = FALSE,
                 filename="heatmap_SingleR_main.png", width = 7, height = 8
)


p2 <- DimPlot(ST227, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
ggplot2::ggsave(filename = "D3_D21_annotated_umap_KJ.pdf",plot = p2, width = 7, height = 6) 
p2
dev.off()

pred.fine <- SingleR(test = ST227.sce.fig, ref = ref, labels = ref$label.fine)
plotScoreHeatmap(pred.fine,
                 clusters = ST227$seurat_clusters, order.by.clusters = TRUE,
                 show.labels = TRUE, show.pruned = FALSE,
                 filename="heatmap_SingleR_fine.pdf", width = 15,height = 18
)

ST227[["SingleR.labels.fine"]] <- pred.fine$labels
Idents(ST227) <- "SingleR.labels.fine"

p3 <- DimPlot(ST227, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
pdf("D3_D21_annotated_umap_fine.pdf", width = 7, height = 6 )
p3
dev.off()
