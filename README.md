What determines whether an individual HSC chooses to self-renew or differentiate? And is this decision HERITABLE - does one HSC's choice influence what its daughter cells do? -> clonal fate bias or lineage commitment.

method: 
1. Add unique DNA barcode to EACH HSC
2. Let them grow in vivo (natural conditions)
3. Sample blood over time
4. Sequence barcode → count which cell types
   came from which HSC

## fig2_plotting_barcodes.r 
  characterise what each barcode does: Over 21 days, which cell types does this HSC generate, and how does output change over time?
  creates a 1D umap, pseudo temporal developmental trajectory for each HSC. it also arranges barcodes from fate pattern A to Bp; do HSCs fall into a distinct fate category or is there a spectrum?

## fig3_cluster_seurat_annotation.r

understand what GENES cause these diff fate biases. If HSC A becomes erythroid-biased and HSC B becomes cDC-biased,what genes are different between them? 

barcode data tells us the outcome [ which cell types were made ] and rna sequencing tells us the mechanism aka which genes were expressed. combining them means we can link genes and fate. 

## working with cosine similarity vs corellation distance based metrics: A recent comparison has
suggested that correlation-based distances may outperform other distance metrics when used with k-means or as the basis for Gaussian kernels (Kim et al, 2018). changing the code for fig 2 and 3 to reflect this change. 
