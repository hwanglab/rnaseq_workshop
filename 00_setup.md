# Setup and Requirements
## Account Information
1. This is our second Bioinformatics workshop, 2018.
1. We will cover basic RNA-seq analyses
1. No prerequisite is required.
1. For those of you registered for this workshop, we already created a Linux account and in advance sent to each of you credential information to access a Linux network at LRI (Lerner Research Institute).

## SSH client
1. Download [MobaXterm](https://download.mobatek.net/1082018070240950/MobaXterm_Portable_v10.8.zip) and install to your CCF laptop (It does not need an administrator account). 

1. HPC host names/ip address
    ```bash
    10.66.64.37
    10.66.64.38
    10.66.64.39
    10.66.64.40
    ```
## Register environment variables
1. Create a working directory, `~/projects/bioinfo_2018July`
1. Open ~/.bashrc and let us add an alias,

	```bash
	alias l="ls -lt"
	alias binf="cd ~/projects/bioinfo_2018July"
	```

1. Also, add the following environment varilables,
	```bash
	export RNA_HOME=~/projects/bioinfo_2018July/rnaseq
	export RNA_EXT_DATA_DIR=$RNA_HOME/RNA_data
	export RNA_DATA_DIR=$RNA_HOME/data
	export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
	export RNA_REFS_DIR=$RNA_HOME/refs
	export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
	export RNA_REF_FASTA=$RNA_REF_INDEX.fa
	export RNA_REF_GTF=$RNA_REF_INDEX.gtf
	export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2
	```
	
1. Relogin and check if the new environment variables are effective! 

## Download the workshop material
1. `binf`
1. `wget xxxx/bioinfo_2018July.tar.gz`
1. `untar bioinfo_2018July.tar.gz`

## Required programs in Linux
The program should be in the path!
1. fastqc
1. trim_galore
1. cutadapt
1. samtools
1. hisat2
1. picard
1. igv.sh
1. htseq-count
1. R (DESeq2,pheatmap,gplots,ggplot2,ggrepel,data.table,stringr,GenomicRanges)
1. PDF viewer
1. python 2.7+
1. perl
1. vim

## Get familiar with
1. tmux
1. copy and paste between Windows and linux terminal
1. vim text editor


# Reference
Some contents/material used in this workshop are borrowed from
https://github.com/griffithlab/rnaseq_tutorial/wiki