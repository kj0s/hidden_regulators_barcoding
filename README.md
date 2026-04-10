# Overview

This project integrates three independent data layers to characterize cell identity, origin, and lineage:

1. **RNA expression** (stored in a Seurat object) — defines the transcriptional state of each cell.  
2. **HTO / SNP donor IDs** — identify the sample or donor of origin.  
3. **SPLINTR lineage barcodes** — track clonal lineage relationships between cells.  

By combining these layers, we can link cell state, biological origin, and lineage history in a single framework.

---

# Methodology

The RNA data is processed using a standard Seurat workflow:
NormalizeData → FindVariableFeatures → ScaleData → RunPCA  → FindNeighbors → FindClusters → RunUMAP


This pipeline produces:

- A **feature space** based on gene expression  
- A **graph structure** connecting cells by similarity  
- **Clusters** representing groups of transcriptionally similar cells  

Each cluster is interpreted as a transcriptional state, which we use as a proxy for cell fate.

## Data filtering and annotation

Before integration, we:

- Remove low-quality cells  
- Exclude ambiguous assignments  
- Filter out unwanted donors  

We then assign donor identities to each cell, ensuring every cell is linked to its biological origin.

## Data integration

Next, we connect the different data layers:

- RNA profiles are linked to lineage barcodes  
- Barcodes are mapped to donor identities  

This results in a unified table where each cell is associated with its expression profile, cluster, donor, and lineage barcode.

---

# Results

The final dataset is organized as a single table with the following structure:

cell
 ├── RNA expression (state)
 ├── cluster (state group)
 ├── donor (sample origin)
 └── barcode (lineage identity)

 
---

# Cell Fate Definitions

We define cell fate using two complementary approaches:

## Phenotypic fate

Cells are grouped into clusters based on gene expression.  
These clusters represent transcriptional states and serve as a proxy for phenotypic cell fate.

## Lineage fate

Lineage fate is derived from shared barcodes across cells, plates, and donors.  
This allows us to track clonal behavior over time and across samples.

We specifically assess:

1. Whether a lineage persists  
2. Whether it expands  
3. Where it appears (across plates or donors)  

---

# Identifying Cell Fates

For each barcode, we:

- Identify all cells carrying that barcode  
- Examine their cluster assignments  
- Determine their donor identities  
- Analyze their distribution across plates  

This approach links lineage information with transcriptional state and sample origin, enabling a more complete view of cell fate.
