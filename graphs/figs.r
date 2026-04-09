# =========================
# 1. Install packages (run once)
# =========================
install.packages(c("dplyr", "ggplot2", "reshape2", "skmeans", "umap"))

# =========================
# 2. Load libraries
# =========================
library(dplyr)
library(ggplot2)
library(reshape2)
library(skmeans)
library(umap)

# =========================
# 3. Load your data
# =========================
common.cpm <- read.table("ST223_common.barcodes.cpm.txt", header = TRUE)
common.cell <- read.table("ST223_common.barcodes.cell.txt", header = TRUE)

# =========================
# 4. Merge replicate plates (cleaned)
# =========================
Samples <- unique(substr(names(common.cell), 1, nchar(names(common.cell)) - 3))

avg <- data.frame(row.names = row.names(common.cell))
new.names <- c()

for (s in Samples) {
  cols <- grep(s, names(common.cell))
  if (length(cols) < 2) next
  
  e <- common.cell[, cols]
  avg <- cbind(avg, rowMeans(e))
  new.names <- c(new.names, s)
}

colnames(avg) <- new.names
avg.common.cell.barcodes <- avg

# =========================
# 5. Log transform
# =========================
data_mat <- log2(avg.common.cell.barcodes + 1)
data_mat <- as.matrix(data_mat)

# =========================
# 6. Choose cluster number (Elbow plot)
# =========================
wssplot <- function(data, nc=15){
  wss <- numeric(nc)
  for (i in 1:nc){
    wss[i] <- sum(kmeans(data, centers=i)$withinss)
  }
  plot(1:nc, wss, type="b", pch=19,
       xlab="Number of clusters", ylab="Within-cluster sum of squares")
}

wssplot(data_mat, nc = 15)

# =========================
# 7. Clustering
# =========================
set.seed(123)
clusters <- skmeans(data_mat, 10)

avg.common.cell.barcodes$cluster <- as.factor(clusters$cluster)

# =========================
# 8. UMAP (FATE UMAP)  → FIGURE B
# =========================
umap_res <- umap(data_mat)

avg.common.cell.barcodes$UMAP1 <- umap_res$layout[,1]
avg.common.cell.barcodes$UMAP2 <- umap_res$layout[,2]

p_fate <- ggplot(avg.common.cell.barcodes,
                 aes(x = UMAP1, y = UMAP2, colour = cluster)) +
  geom_point(size = 1) +
  theme_minimal() +
  ggtitle("Fate UMAP (barcode clustering)")

print(p_fate)

# =========================
# 9. Biomass overlay → FIGURE C
# =========================
long <- melt(avg.common.cell.barcodes,
             id.vars = c("cluster", "UMAP1", "UMAP2"))

long$value <- log2(long$value + 1)

p_biomass <- ggplot(long,
                    aes(x = UMAP1, y = UMAP2, colour = value)) +
  geom_point(size = 0.5) +
  scale_colour_viridis_c() +
  theme_minimal() +
  ggtitle("Fate UMAP with biomass overlay (log2)")

print(p_biomass)

# =========================
# 10. Timepoint extraction (if labels exist)
# =========================
long$day <- case_when(
  grepl("D7", long$variable) ~ "D7",
  grepl("D10", long$variable) ~ "D10",
  grepl("D14", long$variable) ~ "D14",
  grepl("D17", long$variable) ~ "D17",
  grepl("D21", long$variable) ~ "D21",
  TRUE ~ NA_character_
)

long <- long %>% filter(!is.na(day))

# =========================
# 11. Line plots → temporal trends
# =========================
line_data <- long %>%
  group_by(cluster, day) %>%
  summarise(mean_value = mean(value), .groups = "drop")

p_line <- ggplot(line_data,
                 aes(x = day, y = mean_value,
                     group = cluster, colour = cluster)) +
  geom_line(linewidth = 1) +
  theme_minimal() +
  ggtitle("Lineage dynamics over time")

print(p_line)

# =========================
# 12. Save outputs
# =========================
ggsave("Fate_UMAP_clusters.png", p_fate, width = 6, height = 5)
ggsave("Fate_UMAP_biomass.png", p_biomass, width = 6, height = 5)
ggsave("Lineage_trends.png", p_line, width = 6, height = 5)
