# Alignment Visualization
The Integrative Genomics Viewer (IGV) is a high-performance visualization tool for interactive exploration of large, integrated genomic datasets, developed by Broad Institute. It supports a wide variety of data types, including array-based and next-generation sequence data, and genomic annotations. 

## Indexing BAM (.bai)
Before we can view our alignments in the IGV browser we need to index our BAM files. We will use samtools index for this purpose. For convenience later, index all bam files. Run the following command lines (or run `create_bai.sh`). An index file (.bai) will be generated for each BAM file.

```bash

echo $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR
find *.bam -exec echo samtools index {} \; | sh

```

## Visualize alignments
Start IGV on your laptop.

```bash
binf
cd tools/varan-gie_0.2.9
bash ./igv_varan.sh &
```

Load the UHR.bam & HBR.bam files in IGV. Make sure that set Human (hg38) to reference genome sequence. You can load the necessary files (`UHR.bam` and `HBR.bam`) in IGV directly using 'File' -> 'Load from URL'. You may wish to customize the track names as you load them in to keep them straight. Do this by right-clicking on the alignment track and choosing 'Rename Track'.

Q5.1 Go to an example gene locus on chr22: for example, *EIF3L*, *NDUFA6*, and *RBX1* have nice coverage

Q5.2 Differentially expressed genes,
- For example, *SULT4A1* and *GTSE1* are differentially expressed. Are they up-regulated or down-regulated in the brain (HBR) compared to cancer cell lines (UHR)?

- Mouse over some reads and use the read group (RG) flag to determine which replicate the reads come from. What other details can you learn about each read and its alignment to the reference genome.

### Exercise

Try to find a variant position in the RNAseq data:
- HINT: *DDX17* is a highly expressed gene with several variants in its 3 prime UTR.
- Other highly expressed genes you might explore are: *NUP50*, *CYB5R3*, and *EIF3L* (all have at least one transcribed variant).
- Are these variants previously known (e.g., present in dbSNP or gnomad)?
- Load [gnomad vcf file]($RNA_HOME/refs/gnomad.exomes.r2.0.2.sites.liftover.b38.22.vcf) into the current IGV session.
- How should we interpret the allele frequency of each variant? Remember that we have rather unusual samples here in that they are actually pooled RNAs corresponding to multiple individuals (genotypes).

## Something useful
- Save the session
- Snapshot
- [tutorial](https://github.com/griffithlab/rnaseq_tutorial/wiki/IGV-Tutorial) in depth
