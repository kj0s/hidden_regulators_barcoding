## to finish :)
nalyses and clusters cell barcodes by doing barcode based clustering, id the marker genes, com[ares it with transcriptome genes.

4 normalised barcode dataseets: all_barcodes.cpm, all_barcodes.cell, common.cpm, common.cell

extracts unique sample id from col name, avg vals across technical replicates, create a new merged dataset nmed avg.common.cell.barcodes

then log2 transform the data, use wssplot -> function determines the optimal no. of clusters by plotting within sluster sm of squares 

  meaning:  wssplot measures the sum of squared distances between each data point and its assigned cluster centroid; we want to minimise this so we get titght and compact clusters;
  as k increases, wss decreases, introduce concept of ELBOW: the rate of decrease shifts from fast to slow, meaning adding more clusters provides little benefit. visualising this means we plot the wss graph and identify the elbow, then move forward.

  perform k-means clustring = 1: partition the data into 10 specific clusters, unsupervised ml technique. idea is to find 10 distinct centroids of data, useful for this data because orienttion based on previously calculated cosine similarity matters more than the magnitude of any of the vectors. this forms 10 unique clusters.

the barcode clusters are integrated into a single seurat rna object, unique ids combining patient id and barcodes made, clusters joined w left join.

then the cluster identities re set on the seurat objects, the 10 clusters are visualised on a UMAP. this finds differently expressed marker genes per cluster (top 5, top 10.) then it makes a heatmap of the top 10 markers. 

a different approach is used here to compare: seurat clustering. finds marker genes for seurat derived clusers, and generates a comparative heatmap 
  
