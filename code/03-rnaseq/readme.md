# RNA-seq
-   Pre-alignment QC [FastQC v0.11.9](https://github.com/s-andrews/FastQC)
-   Alignment [STAR-2.7.3a](https://github.com/alexdobin/STAR); index with [GENCODE v29lift37](https://www.gencodegenes.org/) genome and annotation; note there is a new run of STAR for sQTL
-   Alignment QC [PicardTools 2.21.7](https://github.com/broadinstitute/picard)
-   Compile FastQC and PicardTools metrics [MultiQC v1.9.dev0](https://github.com/ewels/MultiQC)
```
# In picard/
# -d -dd 1: to keep identical sample ID from different folders
python3 -m multiqc -d -dd 1 Walker/ Obrien Werling_final/ hdbr libd -o all_multiqc
```
-   Quantification [Salmon v1.1.0](https://salmon.readthedocs.io/en/latest/); [GENCODE v33lift37](https://www.gencodegenes.org/) decoys-aware index
-   Compile and import quantifications [Tximport 1.14.0](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)

```{R}
txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE, countsFromAbundance="lengthScaledTPM")
write.table(txi$counts,file="gene.noVersion.scaled.counts.tsv",quote=FALSE, sep='\t')
write.table(txi$abundance,file="gene.noVersion.TPM.tsv",quote=FALSE, sep='\t')

txi.tx <- tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE, countsFromAbundance="lengthScaledTPM")
write.table(txi.tx$counts,file="tx.counts.scaled.tsv",quote=FALSE, sep='\t')
write.table(txi.tx$abundance,file="tx.TPM.tsv",quote=FALSE, sep='\t')
```
