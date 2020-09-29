## RNA-seq
- Pre-alignment QC [FastQC v0.11.9](https://github.com/s-andrews/FastQC)
- Alignment [STAR-2.7.3a](https://github.com/alexdobin/STAR); index with [GENCODE v29lift37](https://www.gencodegenes.org/) genome and annotation
- Alignment QC [PicardTools 2.21.7](https://github.com/broadinstitute/picard)
- Compile FastQC and PicardTools metrics [MultiQC v1.9.dev0](https://github.com/ewels/MultiQC)
- Quantification [Salmon v1.1.0](https://salmon.readthedocs.io/en/latest/); [GENCODE v33lift37](https://www.gencodegenes.org/) decoys-aware index
- Compile and import quantifications [Tximport 1.14.0](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)
- 
