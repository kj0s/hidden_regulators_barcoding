# hidden_reg_self

file 2 works on id cells that contain specific barcodes, match the reference barcodes for gfp, bfp, mcherry proteins, and filters out cells with undefined or unknown barcodes. 
it removes cells based on :
 1. low read count, factor here is less than 1k reads
 2. technical replicates w low correlation aka pearson x < 0.6 [strength and direction of the linear relationship between two continuous variables]
 3. barcodes that appear in only one technical replicate & not confirmed across replicates. [technical replicates are repeated measurements of the same sample.]

it finds VIABLE CELLS -> 
1. if they have matching well characterised barcodes
2. have sufficient depth for sequencing
3. show reproducable results in technical replicates
4. appear in multiple samples, { link to PA and PB?}

then normalises the data using the CPM method for downstream analysis.

