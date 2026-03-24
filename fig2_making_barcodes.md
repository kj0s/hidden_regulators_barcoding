plotting bulk barocde plots!

we are mainly using flourescent protein barcoding

the initial parts are colour gradients for fate biases:
> six colours for diff groups : P 1 2 U V T
> anything with D is a colour for a timepoint

two sister wells are differentiated using Red PA and Blue PB


~line 145 is code start~

read barcode data w all timepoints, 
has an all_timepoints function
  this is just changing sample naming : GFP mcherry etc are changed to the diff groups I U 2 1 etc. 
  does so by adding a PATIENT IDENTITY COLUMN. rules used:
    subsrings made of rows: 
      extract substring from second to third character = PA/PB
      6th character is the protein sample number [1 or 2] {what does this mean?}
      49th to 51st is protein type = gfp, mch or bfp.
  case_when() creates a lookup table which matches protein sample to the patient ids.
  this is all because the experiment has 6 patients, but only 3 flourescent protein markers, 2 samples. we use these as a [component key]. this means they track which genetic bg [protein] in which sample corresonds to which group.

  this is applied to two datasets: cell count data, cpm normalised data, and means we grp barcodes by patient id.
  output: normalised dataset with the names and ids matched. 


  ~ line plot ~
  cell abundance per popn over 5 timepoints, line plots show growth trajectories by plate and population. 

 ~ umap clustering 187 to 290~

the data is loaded and log transformed-> 40 columns as 8 populations across 5 timezones. this is to be changed in your new data.
then apply log2 transformation.

  code from dt_B to dt:
  we split the data into PA vs PB, extract consensus barcode [positions 8-27 row name] and this becomes the composite id.
  then inner join in R to keep values present in both wells aka REPORDUCIBLE CLONES. u only end up analysing clones reliably detected in both duplicates.

then we prepare an expression matrix, run 1Dimentional umap: 
  this is 70 dimentional barcode expression data, it projects it into 1D using cosine similarity. there is a pseudo temporal axis where similar barocodes are nearby

  each barcode has a 70 dimentional fate profile since there are 70 columns. umap first starts in 70D space and finds the 30 nearest neighbours [n_neighbours = 30]. the metric used here is [cosine] and this measures the angles between the fate vectors. 2 barcodes w proportionaly similar outputs [70% A, 30% B] are considered close even if one has 100 cells when the other has 1000. 
  
  ## WHATS A PROPORTIONALLY SIMILAR OUTPUT?

  then, we give UMAP peramaters, [min_dist = 0.8; bandwidth = 50; spread = 10;] min_dist is how much dist barcodes need in 1d to be 'separated'. bandwidth is how quickly probabilities decay with distance, and spread is the expansion factor in the embedding space. 

  ## did these parameters come from trial and error? 

  UMAP initially places all the barcodes at random places on a line. then iteratively adjusts 1d positions to the following:
    barcodes that were neighbours in 70D should stay close in 1D and vice versa. using stochastic gradient descent. 

  each barcode gets a single numeric value representing its position on the pseudo temporal axis. the axis now shows differentiation trajectory. left may be skewed to one fate, middle to one, right to another. 

why 1d? it shows a DEVELOPMENTAL TRAJECTORY for barcodes + Correlates transcriptomic state (RNA UMAP) with fate phenotype (barcode UMAP)

the code then creates individual pdf pages for each barcode showing the stacked area plot of fate bias over 5 timepoints, split by plate

heatmapping:
code sorts barcodes by umap1d position (THEN by cluster), so the heatmap shows them in order of developmental progression for both wells A and B [line 246, 270]

~ lines 730 - 747 ~
create new 1D umap, specific for transcriptome data aka 40 cols of popn data. store this as a separate umap1d coord for each barcode based on gene expression patterns. links fate trajectory to transcriptional state.

fate_umap: Empirical observation of where cells actually went (cells per population per timepoint)
rna_umap: Predicted developmental state based on gene expression

the comparison helps us answer: 
Does transcriptional profiling correctly predict what cells will actually do?
Are there transcriptionally identical cells that make different fate choices?
Which genes drive fate vs. which are just correlated markers?


what does this code answer?
 Do individual clones (barcodes) maintain consistent fate biases?
How reproducible are fate decisions between sister wells?
Which genes drive specific fate outcomes?
How do cell populations expand/contract over time?
Are transcriptional states predictive of barcode fate?













  
