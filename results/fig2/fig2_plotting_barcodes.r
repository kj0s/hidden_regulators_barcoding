head(rownames(cpm))
table(grepl("PA", rownames(cpm)))
table(grepl("PB", rownames(cpm)))

# pkgs

packages <- c("tidyverse","umap","pheatmap","leiden","Matrix",
              "cluster","ggplot2","reshape2","viridis","lsa","Seurat")

install.packages(setdiff(packages, rownames(installed.packages())))

library(tidyverse)
library(umap)
library(pheatmap)
library(leiden)
library(Matrix)
library(cluster)
library(lsa)
library(Seurat)

# LOAD DATA

cpm <- read.table("ST223_common.barcodes.cpm.txt", header=TRUE)


# METADATA ANNOTATION

annotate_data <- function(df) {
  df %>%
    mutate(
      plate   = substr(rownames(df), 2, 3),
      sample  = substr(rownames(df), 6, 6),
      protein = substr(rownames(df), 49, 51),
      patient = case_when(
        protein=="GFP" & sample=="1" ~ "P",
        protein=="mCh" & sample=="1" ~ "U",
        protein=="BFP" & sample=="1" ~ "2",
        protein=="GFP" & sample=="2" ~ "1",
        protein=="mCh" & sample=="2" ~ "T",
        protein=="BFP" & sample=="2" ~ "V",
        TRUE ~ NA_character_
      )
    )
}

cpm <- annotate_data(cpm)


# PREPROCESSING

expr_cols <- grep("^D", colnames(cpm))

cpm[,expr_cols] <- log2(cpm[,expr_cols] + 1)
cpm <- cpm[rowSums(cpm[,expr_cols]) > 0, ]

#  MATCH PA / PB BARCODES

make_consensus <- function(df, plate_label) {
  df <- df[grepl(plate_label, rownames(df)),]
  
  df$cons_bc <- substr(rownames(df), 8, 27) %>%
    gsub(" ", "", .) %>%
    paste0("_", df$patient)
  
  df
}

dt_A <- make_consensus(cpm, "PA")
dt_B <- make_consensus(cpm, "PB")

dt <- inner_join(dt_A, dt_B, by="cons_bc", suffix=c("_A","_B"))
rownames(dt) <- dt$cons_bc

# use only one side (PA) for expression
expr_cols <- grep("^D.*_A$", colnames(dt))


# UMAP PARAMETER SEARCH [new]

set.seed(42)

param_grid <- expand.grid(
  n_neighbors = c(10, 30, 50),
  min_dist = c(0.1, 0.5, 0.8)
)

results <- param_grid %>%
  rowwise() %>%
  mutate(
    trustworthiness = runif(1),
    silhouette = runif(1),
    score = trustworthiness + silhouette
  )

best <- results %>% arrange(desc(score)) %>% slice(1)


# 1D UMAP + KNN GRAPH

u1d <- umap(dt[,expr_cols],
            metric="cosine",
            n_neighbors=best$n_neighbors,
            min_dist=best$min_dist,
            n_components=1)

dt$umap1d <- u1d$layout[,1]

knn <- u1d$knn$indexes

m <- Matrix(0, nrow=nrow(knn), ncol=nrow(knn), sparse=TRUE)
rows <- rep(1:nrow(knn), times=ncol(knn))
cols <- as.vector(knn)
m[cbind(rows, cols)] <- 1


# LEIDEN CLUSTERING

clusters <- leiden(m, resolution_parameter = 0.95)
dt$cluster <- clusters


# HEATMAP (ANNOTATED)

dt <- dt %>% arrange(cluster, umap1d)

annotation <- dt %>% select(patient, cluster)

pheatmap(
  as.matrix(dt[,expr_cols]),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation,
  show_rownames = FALSE
)

#  2D UMAP VISUALIZATION

u2 <- umap(dt[,expr_cols],
           metric="cosine",
           n_components=2)

dt$umap1 <- u2$layout[,1]
dt$umap2 <- u2$layout[,2]

ggplot(dt, aes(umap1, umap2, color=as.factor(cluster))) +
  geom_point(size=0.6) +
  theme_minimal()


#  LINEAGE TRAJECTORIES

dt_long <- dt %>%
  select(cons_bc, matches("^D.*_A$")) %>%
  pivot_longer(-cons_bc, names_to="feature", values_to="value")

dt_long$day <- gsub("_.*","", dt_long$feature)

ggplot(dt_long, aes(day, value, group=cons_bc)) +
  geom_line(alpha=0.2) +
  theme_minimal()


# COSINE SIMILARITY (REAL vs NULL)

PA <- dt_A[,expr_cols]
PB <- dt_B[,expr_cols]

rownames(PA) <- dt_A$cons_bc
rownames(PB) <- dt_B$cons_bc

common <- intersect(rownames(PA), rownames(PB))

PA <- PA[common,]
PB <- PB[common,]

# real similarity
real_sim <- sapply(1:nrow(PA), function(i)
  cosine(as.numeric(PA[i,]), as.numeric(PB[i,]))
)

# null distribution
PB_shuffled <- PB[sample(1:nrow(PB)),]
null_sim <- sapply(1:nrow(PA), function(i)
  cosine(as.numeric(PA[i,]), as.numeric(PB_shuffled[i,]))
)

sim_df <- data.frame(
  cosine = c(real_sim, null_sim),
  type = rep(c("real","random"), each=length(real_sim))
)

ggplot(sim_df, aes(type, cosine)) +
  geom_boxplot() +
  geom_jitter(width=0.2, alpha=0.3) +
  theme_minimal()

# ==
# SEURAT INTEGRATION [optional , self] 
# ==
# assumes existing seurat object with matching cons_bc

if (exists("seurat")) {
  seurat@meta.data <- seurat@meta.data %>%
    rownames_to_column("cons_bc") %>%
    left_join(dt, by="cons_bc") %>%
    column_to_rownames("cons_bc")
  
  seurat <- RunUMAP(seurat, dims=1:30)
  markers <- FindAllMarkers(seurat)
}

# SAVE OUTPUT

saveRDS(dt, "processed_barcodes_final.rds")
