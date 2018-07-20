# Alignment Free Expression Estimation (Kallisto)
In this session, we will quantify **transcript**-level expression analysis using *Kallisto*. For more information on Kallisto, refer to the [Kallisto project page](https://pachterlab.github.io/kallisto/about.html) and [Kallisto manual page](https://pachterlab.github.io/kallisto/manual.html).

## Obtain transcript sequences in FASTA format
Note that we already have Fasta sequences for the reference genome sequence from earlier in the RNA-seq tutorial. However, Kallisto works directly on target *cDNA/transcript* sequences. Remember also that we have transcript models for genes on chromosome 22. These transcript models were downloaded from Ensembl in GTF format. This GTF contains a description of the coordinates of exons that make up each transcript but it does not contain the transcript sequences themselves. So currently we do not have the transcript sequences needed by Kallisto. There are many places we could obtain these transcript sequences from. For example, we could download them directly in Fasta format from the [Ensembl FTP site](http://www.ensembl.org/info/data/ftp/index.html).

To allow us to compare Kallisto results to expression results from DESeq2, we will create a custom Fasta file that corresponds to the transcripts we used for the DESeq2 analysis.

#### How can we obtain these transcript sequences in Fasta format?

We will use the script `gtf_to_fasta` from the `tophat` package to generate a Fasta sequence from our GTF file. This approach is convenient because it will also include the sequences for the ERCC spike in controls. This allows us to generate Kallisto abundance estimates for those features as well.

```bash
cd $RNA_HOME/refs
../tools/gtf_to_fasta $RNA_REF_GTF $RNA_REF_FASTA chr22_ERCC92_transcripts.fa
```

Use `less` to view the output file `chr22_ERCC92_transcripts.fa`. Each genomic sequence is a transcript after intronic regions are spliced out and only exonic regions are combined together. Each gene has many isoforms.

Note that this file has messy transcript names. Use the following hairball perl one-liner to tidy up the header line for each Fasta sequence

```bash
cd $RNA_HOME/refs
cat chr22_ERCC92_transcripts.fa | \
perl -ne 'if ($_ =~/^\>\d+\s+\w+\s+(ERCC\S+)[\+\-]/){print ">$1\n"}elsif($_ =~ /\d+\s+(ENST\d+)/){print ">$1\n"}else{print $_}' \
> chr22_ERCC92_transcripts.clean.fa

wc -l chr22_ERCC92_transcripts*.fa

```

View the resulting 'clean' file using `less chr22_ERCC92_transcripts.clean.fa`. To view the end of this file use `tail chr22_ERCC92_transcripts.clean.fa`. Note that we have one Fasta record for each Ensembl transcript on chromosome 22 and we have an additional FASTA records for each ERCC spike-in sequence.

Create a list of all transcript IDs for later use:

```bash
cd $RNA_HOME/refs
cat chr22_ERCC92_transcripts.clean.fa | \
	grep ">" | \
	perl -ne '$_ =~ s/\>//; print $_' | \
	sort | \
	uniq > \
	transcript_id_list.txt
```

## Build a Kallisto transcriptome index
Remember that Kallisto does not perform *alignment* or use the original reference genome sequence. Instead, it performs *pseudoalignment* to determine the *compatibility* of reads with targets (transcript sequences in this case). However, similar to alignment algorithms like Tophat or STAR, Kallisto still requires an **index** on the transcript sequences to assess this compatibility efficiently and quickly.

To generate this index, execute:
```bash
cd $RNA_HOME/refs
mkdir kallisto
cd kallisto
kallisto index --index=chr22_ERCC92_transcripts_kallisto_index ../chr22_ERCC92_transcripts.clean.fa
```

## Generate abundance estimates for all samples using Kallisto
As we did with `htseq-count` and `DESeq2`, we will generate transcript abundances for each of our demonstration samples using `Kallisto`.

Run this series of commands:
```bash
echo $RNA_DATA_TRIM_DIR

cd $RNA_HOME/expression
if [ ! -d slueth/input ]; then
	mkdir -p slueth/input
fi

cd slueth/input

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep1_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep2_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep3_ERCC-Mix1 --threads=4 --plaintext $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep1_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep2_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep3_ERCC-Mix2 --threads=4 --plaintext $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

```

### Expectation-maximization algorithm to model transcript quantification used in kallisto

The principle of EM is nicely illustrated by Lior Pachter in his transcript quantification review. Suppose, as shown on the image below, there are three transcripts (green, red, and blue). There are five reads associated with these transcripts. One read (*d*) is unique to the red transcript, while others correspond to two (*b, c, e*) or three (*a*) transcripts. The EM is an iterative procedure. In the first round, transcript abundances are initialized as equal (0.33 each as there are three transcripts). During expectation, reads are apportioned across transcripts based on these abundances. Next, during the maximization step, transcript abundances are re-calculated as follows:

For the red transcript we sum up fraction of each read as `0.33 + 0 + 0.5 + 1 + 0.5` for reads *a, b, c, d, and e,* respectively. We now divide this by the sum of read allocations for each transcript as `2.33 + 1.33 + 1.33` for red, green, and blue transcripts respectively.

All three transcript calculations will look like this:
- Red transcript normalized abundance = 2.33/4.99 = 0.47
- Green transcript normalized abundance = 1.33/4.99 = 0.27
- Blue transcript normalized abundance = 1.33/4.99 = 0.27

During next expectation stage, reads are re-apportioned across transcripts with the updated priors and the procedure is repeated until convergence:

![em](images/em.png)

### Merge kallisto output files
Create a single TSV file that has the TPM abundance estimates for all six samples.

```bash
cd $RNA_HOME/expression/sleuth/input

paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv

ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv

cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2

mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv

rm -f header.tsv
```

Take a look at the final kallisto result file we created:

```bash
head transcript_tpms_all_samples.tsv
tail transcript_tpms_all_samples.tsv
```

The read count file can be analyzed in DESeq2 for DE analysis in general.

## Perform DE analysis of Kallisto expression estimates using Sleuth

We will now use Sleuth to perform a differential expression analysis on the full chr22 dataset produced above. Sleuth is a companion tool that starts with the output of Kallisto, performs DE analysis, and helps you visualize the results.

Regenerate the Kallisto results using the HDF5 format and 100 rounds of bootstrapping (both required for Sleuth to work).

```bash
cd slueth/input

kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep1_ERCC-Mix1 --threads=4 $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep2_ERCC-Mix1 --threads=4 $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=UHR_Rep3_ERCC-Mix1 --threads=4 $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep1_ERCC-Mix2 --threads=4 $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep2_ERCC-Mix2 --threads=4 $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz

kallisto quant -b 100 --index=$RNA_HOME/refs/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=HBR_Rep3_ERCC-Mix2 --threads=4 $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```


## Install R library module
Sleuth is an R package so the following steps will occur in an `R` session. The following section is an adaptation of the [sleuth getting started tutorial](https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html).

The `sleuth` R module may be not available in your R library, so et us install the module in R.

```r
source("http://bioconductor.org/biocLite.R")
biocLite("devtools")    # only if devtools not yet installed
biocLite("remotes")
biocLite("pachterlab/sleuth")
library(sleuth) #check if the module was installed successfully
```

We would like perform PC analysis and additionally estimate the expressions on one gene (e.g., *SYNGR1-203*/ENST00000328933).

Copy a sleuth R script and perform a differential expression analysis on the RNA transcript.

```bash
cd $RNA_HOME/expression/slueth
mkdir results
cp /Informatics_workshop/rnaseq/08_kallisto/Tutorial_KallistoSleuth.R ./
Rscript Tutorial_KallistoSleuth.R
```

```r
#load sleuth library
suppressMessages({
  library("sleuth")
})

#set input and output dirs
datapath = "input"
resultdir = 'results'

#create a sample to condition metadata description
kal_dirs <- list.dirs(datapath,recursive=FALSE)
print(kal_dirs)
sample = c("HBR_Rep1_ERCC-Mix2", "HBR_Rep2_ERCC-Mix2", "HBR_Rep3_ERCC-Mix2", "UHR_Rep1_ERCC-Mix1", "UHR_Rep2_ERCC-Mix1", "UHR_Rep3_ERCC-Mix1")
condition = c("HBR", "HBR", "HBR", "UHR", "UHR", "UHR")
s2c = data.frame(sample,condition)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

#run sleuth on the data
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

#summarize the sleuth results and view 20 most significant DE transcripts
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

#plot an example DE transcript result
p1 = plot_bootstrap(so, "ENST00000328933", units = "est_counts", color_by = "condition")
p2 = plot_pca(so, color_by = 'condition')

#Print out the plots created above and store in a single PDF file
pdf(file=file.path(resultdir,"SleuthResults.pdf"))
print(p1)
print(p2)
dev.off()
```

<!---
## Compare transcript and gene abundance estimates from Kallisto to isoform abundance estimates from StringTie and counts from HtSeq-Count
How similar are the results we obtained from each approach?

We can compare the expression value for each Ensembl transcript from chromosome 22 as well as the ERCC spike in controls.

To do this comparison, we need to gather the expression estimates for each of our replicates from each approach. The Kallisto transcript results were neatly organized into a single file above. For Kallisto gene expression estimates, we will simply sum the TPM values for transcripts of the same gene. Though it is 'apples-to-oranges', we can also compare Kallisto expression estimates to the raw read counts from HtSeq-Count (but only at the gene level in this case). The following R script will pull together the various expression matrix files we created in previous steps and create some visualizations to compare them in the gene level.

First, create the gene version of the Kallisto TPM matrix

```bash

cd $RNA_HOME/expression/kallisto
wget https://raw.githubusercontent.com/griffithlab/rnaseq_tutorial/master/scripts/kallisto_gene_matrix.pl
chmod +x kallisto_gene_matrix.pl
./kallisto_gene_matrix.pl --gtf_file=$RNA_HOME/refs/chr22_with_ERCC92.gtf  --kallisto_transcript_matrix_in=transcript_tpms_all_samples.tsv --kallisto_transcript_matrix_out=gene_tpms_all_samples.tsv

```

Now load files and summarize results from each approach in R
```bash

cd $RNA_HOME/expression
Rscript Tutorial_comparisons.R

```

--->
