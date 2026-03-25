#Annotation
#Use singleR and Immgen
BiocManager::install("SingleR")
BiocManager::install("celldex")
#finding reference dataset in celldex
library(celldex)
library(SingleR)
library(Seurat)
library(dplyr)
surveyReferences()
# installing devtools, devtools::install_github('immunogenomics/presto')
# use pak:pak('immunogenomics/presto')

setwd("/vast/projects/Sisseq/")

untar("human-haematopoiesis-sis-seq.tar.gz")

HHP_sisseq <- readRDS("human-haematopoiesis-sis-seq.rds")

HHP_sisseq <- RunUMAP(HHP_sisseq, nn.name = "weighted.nn", n.neighbors = 30L, reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
HHP_sisseq <- FindClusters(HHP_sisseq, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

DefaultAssay(HHP_sisseq) <- 'sketch_RNA'
RNA_markers <- FindAllMarkers(HHP_sisseq, only.pos = TRUE)

DefaultAssay(HHP_sisseq) <- 'sketch_ADT'
ADT_markers <- FindAllMarkers(HHP_sisseq, only.pos = TRUE)

RNA_markers  %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

write.table(top10, file= "RNA_markers.txt", sep="\t", row.names=FALSE, col.names=TRUE)

DoHeatmap(HHP_sisseq, features = top10$gene) + NoLegend()

#Remove contaminating stromal cells
HHP_sisseq <- subset(HHP_sisseq, idents = c("19", "29"), invert = TRUE)

#Rename clusters
new.cluster.ids <- c("Ery", "cDC2", "Ery", "Ery", "CD14+ Mono", "T cells", "Pre-B cells", "Mast Cells", "Neutrophils", "Megakaryocytes", "Megakaryocytes", "DC3", "Eosinophils", "HSC", "NK cells", "Basophils", "Mast Cells", "cDC1", "CD16+ Mono", "mregDC", "GMP", "cDC2", "NK cells", "Macrophages", "T cells", "cDC2", "Megakaryocytes", "Eosinophils")

# Assign new names
HHP_sisseq$celltype <- plyr::mapvalues(Idents(HHP_sisseq), 
                                  from = levels(Idents(HHP_sisseq)), 
                                  to = new.cluster.ids)
Idents(HHP_sisseq) <- HHP_sisseq$celltype

p1 <- DimPlot(HHP_sisseq, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
ggplot2::ggsave("annotated_umap.pdf", plot = p1, width = 7, height = 6, units = "in") 
p1
dev.off()

#Compare with BM dataset
HHP_sisseq.sce.fig <- as.SingleCellExperiment(HHP_sisseq, assay = "sketch_RNA")

# put the nvstb cells in directory!

ref <- fetchReference("novershtern_hematopoietic", "2024-02-26")
pred.main <- SingleR(test = HHP_sisseq.sce.fig, ref = ref, labels = ref$label.main)

HHP_sisseq[["SingleR.labels"]] <- pred.main$labels
Idents(HHP_sisseq) <- "SingleR.labels"

plotScoreHeatmap(pred.main,
                 clusters = HHP_sisseq$celltype, order.by.clusters = TRUE,
                 show.labels = TRUE, show.pruned = FALSE,
                 filename="heatmap_SingleR_main.png", width = 7, height = 8
)


p2 <- DimPlot(HHP_sisseq, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
ggplot2::ggsave(filename = "D3_D21_annotated_umap_KJ.pdf",plot = p2, width = 7, height = 6) 
p2
dev.off()

pred.fine <- SingleR(test = HHP_sisseq.sce.fig, ref = ref, labels = ref$label.fine)
plotScoreHeatmap(pred.fine,
                 clusters = HHP_sisseq$seurat_clusters, order.by.clusters = TRUE,
                 show.labels = TRUE, show.pruned = FALSE,
                 filename="heatmap_SingleR_fine.pdf", width = 15,height = 18
)

HHP_sisseq[["SingleR.labels.fine"]] <- pred.fine$labels
Idents(HHP_sisseq) <- "SingleR.labels.fine"

p3 <- DimPlot(HHP_sisseq, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
pdf("D3_D21_annotated_umap_fine.pdf", width = 7, height = 6 )
p3
dev.off()
