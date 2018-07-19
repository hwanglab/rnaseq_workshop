# RNA-seq Data

The test data consists of two commercially available RNA samples: [Universal Human Reference](https://github.com/griffithlab/rnaseq_tutorial/wiki/ResourceFiles/UHR.pdf) (UHR) and [Human Brain Reference](https://github.com/griffithlab/rnaseq_tutorial/wiki/ResourceFiles/HBR.pdf) (HBR). The UHR is total RNA isolated from a diverse set of 10 cancer cell lines. The HBR is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old. In this workshop, we consider UHR as a **control** group and HBR for an **experimental** group.

## ERCC (within a sample and between samples)
We added an aliquot of the [ERCC ExFold RNA Spike-In Control Mixes](https://github.com/griffithlab/rnaseq_tutorial/wiki/ResourceFiles/ERCC.pdf) to each sample. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies). This range allows us to test the degree to which the RNA-seq assay (including all laboratory and analysis steps) accurately reflects the relative abundance of transcript species within a sample. There are two 'mixes' of these transcripts to allow an assessment of differential expression output between samples if you put one mix in each of your two comparisons. In our case, Mix1 was added to the UHR sample, and Mix2 was added to the HBR sample.

![ercc_mix](images/ERCC-f1.jpg)
**Transcript molar ratios in ERCC Spike-In Mixes**

https://www.thermofisher.com/order/catalog/product/4456740

## Two triplicates
We have 3 complete experimental replicates for each sample. This allows us to assess the technical variability of our overall process of producing RNA-seq data in the lab.

For all libraries, we prepared low-throughput (Set A) [TruSeq Stranded Total RNA Sample Prep Kit](http://products.illumina.com/products/truseq_stranded_total_rna_sample_prep_kit.html) libraries with Ribo-Zero Gold to remove both cytoplasmic and mitochondrial rRNA.

Triplicate, indexed libraries were made starting with 100ng Agilent/Stratagene Universal Human Reference total RNA and 100ng Ambion Human Brain Reference total RNA. The Universal Human Reference replicates received 2 ul of 1:1000 ERCC Mix 1. The Human Brain Reference replicates received 1:1000 ERCC Mix 2. The libraries were quantified with KAPA Library Quantification qPCR and adjusted to the appropriate concentration for sequencing. The triplicate, indexed libraries were then pooled prior to sequencing. Each pool of three replicate libraries was sequenced across 2 lanes of a HiSeq 2000 using paired-end sequence chemistry with 100bp read lengths.

So to summarize we have:

- UHR + ERCC Spike-In Mix1, Replicate 1
- UHR + ERCC Spike-In Mix1, Replicate 2
- UHR + ERCC Spike-In Mix1, Replicate 3
- HBR + ERCC Spike-In Mix2, Replicate 1
- HBR + ERCC Spike-In Mix2, Replicate 2
- HBR + ERCC Spike-In Mix2, Replicate 3

Each data set has a corresponding pair of FastQ files (read 1 and read 2 of paired-end reads).

The reads are paired-end 101-mers generated on an Illumina HiSeq instrument. The test data has been pre-filtered for reads that appear to map to chromosome 22.

Go to $RNA_DATA_DIR. Unpack the test data. You should see 6 sets of paired-end FASTQ files. One for each of our sample replicates above. We have 6 pairs (12 files) because in FASTQ format, read 1 and read 2 of each read pair (fragment) are stored in separate files.

```bash
cd $RNA_DATA_DIR
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
l
```

In the data directory, view the first two read records of a file (in FASTQ format each read corresponds to 4 lines of data)

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
```

Q3.1 Identify the following components of each read: read name, read sequence, and quality string

How many reads are there in the first library? Decompress file on the fly with 'zcat', pipe into 'grep', search for the read name prefix and count a total number of reads.


```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -Pc "^\@HWI"
```

## QC check of RNA-Seq reads
Q3.2 Run FastQC on your RNA-seq FASTQ files. First, create an output directory `raw_fastqc` and store report files to the output directory.

Q3.3 Open *_fastqc.html file in a web browser to view the FastQC report and investigate the source/explanation for over-represented sequences.

## Adapter trimming
Use trim_galore to trim sequence adapter from the read FASTQ files. The output of this step will be trimmed FASTQ files for each data set.

Refer to the [trim_galore manual](https://github.com/FelixKrueger/TrimGalore) for a more detailed explanation

Prepare a subdirectory to store trimmed read files if not exist,
```bash
mkdir $RNA_DATA_TRIM_DIR
``` 

Run the following bash shell script (trim_galore.sh) and run it.
 
```bash
#!/bin/bash -l

# module load trim_galore/0.4.4 # make trim_galore available in path
# module load FastQC/0.11.5

cd $RNA_DATA_DIR # change to RNA-seq read directory

for read_base in `ls *read1.fastq.gz | rev | cut -c 11- | rev | uniq` # list base file names to process
do
    echo "running trim_galore on $read_base ..."
    trim_galore -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -a2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --paired --retain_unpaired --gzip --fastqc -o $RNA_DATA_TRIM_DIR ${read_base}1.fastq.gz ${read_base}2.fastq.gz
    
    echo "renaming qc read_1 ..." 
    mv ${RNA_DATA_TRIM_DIR}/${read_base}1_val_1.fq.gz ${RNA_DATA_TRIM_DIR}/${read_base}1.fastq.gz
    
    echo "renaming qc read_2 ..."
    mv ${RNA_DATA_TRIM_DIR}/${read_base}2_val_2.fq.gz ${RNA_DATA_TRIM_DIR}/${read_base}2.fastq.gz
    
    echo "Done."
done
```

** cutadapt should be in the path

Options used:
- specified two adapters to trim if any
- retain unpaired even when one pair becomes too short after trimming
- gzip output
- run fastqc on trimmed reads
