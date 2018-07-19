# Alignment

## HISAT2 alignment
HISAT2 is a reference genome-based, splice-aware, NGS short-read aligner. It uses a graph-based alignment strategy and has overcome HISAT and TopHat2 as a standard. HISAT2's output is SAM/BAM files for each data set. RNA-seq BAM files can be very useful in many applications, including genomic variants (e.g., RNA-editing), gene fusion, and gene expression analysis, as well as for discovery of *de novo* isoforms. Now that both RNA-seq paired-end reads files and HISAT2 reference index are available, let us perform alignments with HISAT2 to the genome and transcriptome.

Refer to the HISAT2 manual for a more detailed explanation:

* https://ccb.jhu.edu/software/hisat2/manual.shtml


### Setup

First, if not present, begin by making the corresponding output directory for our alignment results.
```bash
l $RNA_ALIGN_DIR # Check directory contents
mkdir -p $RNA_ALIGN_DIR #If not found, create
cd $RNA_ALIGN_DIR #Change to directory
```

### HISAT2 basic usage:

```bash
hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
```

Extra options specified below:

* `-p 4` - Tells HISAT2 to use eight CPUs for bowtie alignments.
* `--RNA-strandness RF` - Specifies strandness of RNAseq library.
  * We will specify RF since the TruSeq strand-specific library was used to make these libraries. See <a href="https://github.com/griffithlab/rnaseq_tutorial/blob/master/manuscript/supplementary_tables/supplementary_table_5.md">here</a> for details.
* `--rg-id $ID` - Specifies a read group ID that is a unique identifier.
* `--rg SM:$SAMPLE_NAME` - Specifies a read group sample name. This, together with rg-id, will allow you to determine which reads came from which sample in the merged bam later on.
* `--rg LB:$LIBRARY_NAME` - Specifies a read group library name. This, together with rg-id, will allow you to determine which reads came from which library in the merged bam later on.
* `--rg PL:ILLUMINA` - Specifies a read group sequencing platform.
* `--rg PU:$PLATFORM_UNIT` - Specifies a read group sequencing platform unit.
  * Typically this consists of FLOWCELL-BARCODE.LANE
* `--dta` - Reports alignments tailored for transcript assemblers.
* `-x </path/to/hisat2/index>` - The HISAT2 index filename prefix built earlier including splice sites and exons (minus the trailing .X.ht2).
* `-1 </path/to/read1.fastq.gz>` - The read 1 FASTQ file.
 * Optionally: gzip(.gz) or bzip2(.bz2) compressed.
* `-2 </path/to/read2.fastq.gz>` - The read 2 FASTQ file.
 * Optionally: gzip(.gz) or bzip2(.bz2) compressed.
* `-S </path/to/output.sam>` - The output SAM format text file of alignments.

## HISAT2 run on the triplicates
Let us run HISAT2 with the following commands:
```bash
hisat2 -p 4 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_TRIM_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools sort -o ./UHR_Rep1.bam
hisat2 -p 4 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_TRIM_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools sort -o ./UHR_Rep2.bam
hisat2 -p 4 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_TRIM_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools sort -o ./UHR_Rep3.bam

hisat2 -p 4 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_TRIM_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools sort -o ./HBR_Rep1.bam
hisat2 -p 4 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_TRIM_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools sort -o ./HBR_Rep2.bam
hisat2 -p 4 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_TRIM_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz | samtools sort -o ./HBR_Rep3.bam
```
 These commands are also in the bash shell script `hisat2.sh`, located in the `$RNA_ALIGN_DIR` directory.

**Important:**

In the above alignments, we are treating each library as an independent data set. If you had multiple lanes of data for a single library, you could align them all together in one HISAT2 command. Similarly, you might combine technical replicates into a single alignment run (perhaps after examining them and removing outliers).

To combine multiple lanes, you would provide all the read1 files as a comma-separated list for the `-1` input argument, then all read2 files as a comma-separated list for the `-2` input argument (where both lists have the same order).
You can also combine individual BAM files to merge a bam file. This is the approach we will take.

### HISAT2 Alignment Summary

HISAT2 generates a summary of the alignments printed to the terminal. Notice the number of total reads, reads aligned and various metrics regarding how the reads aligned to the reference.

### Merge HISAT2 BAM files

Make a single BAM file combining all UHR data and another for all HBR data. This could be done in several ways, such as:
* `samtools merge`
* `bamtools merge`
* using picard-tools (see below)

We chose the third method because it did the best job of merging the bam header information. NOTE: `sambamba` also retains header info.

```bash
cd $RNA_HOME/alignments/hisat2
java -Xmx2g -jar $RNA_HOME/tools/picard.jar MergeSamFiles OUTPUT=UHR.bam INPUT=UHR_Rep1.bam INPUT=UHR_Rep2.bam INPUT=UHR_Rep3.bam
java -Xmx2g -jar $RNA_HOME/tools/picard.jar MergeSamFiles OUTPUT=HBR.bam INPUT=HBR_Rep1.bam INPUT=HBR_Rep2.bam INPUT=HBR_Rep3.bam
```

Count the alignment (BAM) files to make sure all were created successfully (you should have 8 total).

```bash
ls -l *.bam | wc -l
ls -l *.bam
```
#### Q4.1 By merging the triplicate BAM files, do we lose the information from which BAM file each read originates?
#### Q4.2 Consider two extra post-processing steps 1) MAPQ>=1 and 2) collapse duplicated reads within an individual sample

### Up next
[Alignment Visualization](05_postalignment-visualization.md)
