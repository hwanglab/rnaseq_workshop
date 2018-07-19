# Indexing
 
## Create a HISAT2 index

Create a HISAT2 index for chr22 and the ERCC spike-in sequences. HISAT2 can incorporate exons and splice sites into the index file for alignment. First create a splice site file, then an exon file.  Finally, make the aligner FM index.

To learn more about how the HISAT2 indexing strategy is distinct from other NGS aligners, refer to the [HISAT publication](https://www.ncbi.nlm.nih.gov/pubmed/25751142).

Run the following commands. A bash script file (`$RNA_REFS_DIR/hisat2_index.sh`) is also available.

```bash

cd $RNA_REFS_DIR
hisat2_extract_splice_sites.py $RNA_REF_GTF > $RNA_REFS_DIR/splicesites.tsv
hisat2_extract_exons.py $RNA_REF_GTF > $RNA_REFS_DIR/exons.tsv
hisat2-build -p 8 --ss $RNA_REFS_DIR/splicesites.tsv --exon $RNA_REFS_DIR/exons.tsv $RNA_REF_FASTA $RNA_REF_INDEX

```
Check the output file from the command above

Note that all environment variables used in the commands above are defined in `~/bashrc`. Try this to see which values are assigned to the variable.

```bash
echo $RNA_REF_GTF
echo $RNA_REFS_DIR
``` 

** HISAT2 should be in the path
