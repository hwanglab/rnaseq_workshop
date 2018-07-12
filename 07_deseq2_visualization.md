# Differential Expression

In this tutorial, you will:

* Make use of the raw counts you generated in the previous session using htseq-count
* Use DESeq2, a bioconductor package designed specifically for differential expression of count-based RNA-seq data
* Visualize differentially expressed genes (features)
* R programming language

## Normalization for DE
 
RNA-seq is often used to compare one tissue type to another (liver vs. spleen), or experimental group vs. control group sample. There are many specific genes transcribed in liver but not in the spleen or vice versa. This is called "a difference in library composition".

Let us overview an R script that performs RNA-Seq differential analyses.

```bash

cd $RNA_RC_DIR
less runDESeq2FromFeatureCount2.r

```

### Loading matrix
* L34: We load the count matrix,
* L37: Assign which sample belongs to which group.
* L39: Retain only genes where half of samples have at least 10 reads,

### DESeq2
DESeq2 (or edgeR) handles both library size and composition issues using log normalization and negative bionomial distribution model. Refer to DESeq2 paper for more details.

* L41: DESeq2 function rlog() normalizes the count matrix,
	* Take log of each value in the count matrix
	* Compute the average log value across samples for each gene 
	* Substract the log value by the average value
	
* L62: DESeq()
	* estimating size factors and dispersions
	* negative bionomial distribution model to fit
	* significance testing

Open DESeq2 report file or we will analyze this via heatmap below.

### Hierachical Clustering
* L42: Computes a distance between each pair of sample
* L47: Perform hierachical clustering

### PCA Plot
Principal Component Analysis (PCA) is a method of dimensionality reduction/feature extraction that transforms the data from a d-dimensional space into a new coordinate system of dimension p, where p<=d.
Here, we use only two PCs. 

### Heatmap
* L70: Load DESeq2 output file
* L71: Focusing significant differentially-expressed genes
* L74: Load log2(TPM2) file
* L82: Generate a heatmap from a log2(TPM2) counts only reported in the refined DESeq2 output table.

Let us run the R script to generate all table and fiugres for DE analysis
```bash
Rscript ./runDESeq2FromFeatureCount2.r
```
Check files in the output directory, `expr_output`

### ERCC DE Analysis

Let us compare the observed versus expected differential expression estimates for the ERCC spike-in RNAs:

```bash

cd $RNA_RC_DIR
Rscript ./Tutorial_ERCC_DE.R ERCC_Controls_Analysis.txt expr_output/DESeq2_UHR_HBR.tsv

```
Tutorial_ERCC_DE.pdf

### Resource
- https://cran.r-project.org/doc/contrib/Torfs+Brauer-Short-R-Intro.pdf
