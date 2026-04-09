purpose of code is to identify which genes drive the diff fate biases. this clusters cells and annotates their cell types using seurat.


load preprocessed single cell data, and calc mitochondrial gene percentage. filter out low quality cells.

## rna and protein integration
there is multimodal data integration here, youre comparing based on more than one metric.

process the rna data separately, principle component analysis on top genes to reduce dimensions, and then do the same for adt protein data.

## how to combine the cell specific weights for modalities? RNA.weight?
process both to create a weighted neighbour graph by:
1. build separate nearest neighbour graphs [30 nearest]
2. end up w 2 diff neighbourhood definitions for each cell
3. function learns cell specific modality weight that combine the two modalities
4. this function then weighs the two cells according to their learned weights, creates merged knn.
OUTPUTS: wsnn graph (weighted shared nearest neighbor) ; RNA.weight metadata, how much each cell relies on RNA vs. protein

findallmarkers() -> runs test on every gene in every cluster to ask which genes are uniquely HIGH in this cluster compared to all other cells?

process:
1. Split cells into two groups:
    a. Group A: Cells in this cluster (e.g., B cells)
    b. Group B: All other cells (not B cells)
  - For each gene, run a statistical test (default: Wilcoxon rank-sum test):

2. Compare expression level of gene X in Group A vs Group B
    a. Calculate p-value: "Is this difference statistically significant?"
    b. Calculate log2 fold-change: "How much higher is it in Group A?"
  - Filter results using your parameters:

only.pos = TRUE → Only keep genes that are HIGHER in the cluster (positive fold-change)
min.pct = 0.25 → Gene must be expressed in at least 25% of cells in that cluster
logfc.threshold = 0.25 → Fold-change must be at least 1.2x higher
=> Output a ranked list of genes per cluster











