[![DOI](https://zenodo.org/badge/268661500.svg)](https://zenodo.org/badge/latestdoi/268661500)

# Developing brain mega xQTL
**Supplemental tables, all QTL summary statistics, and other extended data are available at https://doi.org/10.7303/syn50897018.5**
## 1: Munging data
## 2: Genotype
### 2-1: Pre-imputation
- Run `plinkQC` on data as a sanity check
- First apply PLINK filters, then split by chromosome and sort
    - Walker data is already filtered; split by chromosome and impute
    - For all the other datasets, we applied the same filters that the Walker data used `--hwe 1e-6 --maf 0.01 --mind 0.10 --geno 0.05`
    - Note: for HDBR, we used `--mind 0.3`; for LIBD, we fixed strand flips by running an extra step of [conform-gt](https://faculty.washington.edu/browning/conform-gt.html), which automatically splits the data by chromosome
### 2-2: Post-imputation
- Scripts in `prelim/`: inputs are imputed genotype files downloaded from Michigan Imputation Server; concatenate by chromosomes, index, filter by R2, and take the ***intersection*** of high impute quality variants across datasets
    - Note: except for Walker data, we applied R2>.3 filter during imputation; so here we only applied R2>.3 on Walker imputed data and intersected with the other datasets
- `ancestry.ipynb`: infer data ancestry, make plots
- `IBD.ipynb`: relatedness check
- `Snakefile`
## 3: RNA-seq
-   Pre-alignment QC [FastQC v0.11.9](https://github.com/s-andrews/FastQC)
-   Alignment [STAR-2.7.3a](https://github.com/alexdobin/STAR), index with [GENCODE v29lift37](https://www.gencodegenes.org/) genome and annotation
    - Note: there is a new run of STAR for sQTL
-   Alignment QC [PicardTools 2.21.7](https://github.com/broadinstitute/picard)
-   Compile FastQC and PicardTools metrics [MultiQC v1.9.dev0](https://github.com/ewels/MultiQC)
```
# In picard/
# -d -dd 1: to keep identical sample ID from different folders
python3 -m multiqc -d -dd 1 Walker/ Obrien Werling_final/ hdbr libd -o all_multiqc
```
-   Quantification [Salmon v1.1.0](https://salmon.readthedocs.io/en/latest/), [GENCODE v33lift37](https://www.gencodegenes.org/) decoys-aware index
-   Compile and import quantifications [Tximport 1.14.0](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)

```{R}
txi <- tximport(files, type="salmon", tx2gene=tx2gene, dropInfReps=TRUE, countsFromAbundance="lengthScaledTPM")
write.table(txi$counts,file="gene.noVersion.scaled.counts.tsv",quote=FALSE, sep='\t')
write.table(txi$abundance,file="gene.noVersion.TPM.tsv",quote=FALSE, sep='\t')

txi.tx <- tximport(files, type="salmon", txOut=TRUE, dropInfReps=TRUE, countsFromAbundance="lengthScaledTPM")
write.table(txi.tx$counts,file="tx.counts.scaled.tsv",quote=FALSE, sep='\t')
write.table(txi.tx$abundance,file="tx.TPM.tsv",quote=FALSE, sep='\t')
```
-   Sample swap check: 
    + [VerifyBamID](https://genome.sph.umich.edu/wiki/VerifyBamID) (slow. Use `--smID` to add subject ID to BAM sequence file)
    + `check.ipynb`: called SNP from BAM, merged with imputed genotype (Mike)
## 4: xQTL
### 4-1: cis-eQTL
- `ancestry.ipynb`
- `combat-seq.ipynb`
- `decon.ipynb`: cell type specific and interacting analysis
- `eqtl_analysis.ipynb`: identify optimal #HCP in covariates, gene expression PCA, dTSS, etc.
- `fetal_adult.ipynb`
- `func_enrich.ipynb`: functional enrichment analysis of QTL
- `metadata.ipynb`: plot data age, sex, infer NA sex, etc.
- `module_eigengene.ipynb`
- `paintor.ipynb`: PAINTOR multi-ethnic fine-mapping 
- `pLI.ipynb`
- `sex_specific.ipynb`
- `susie.ipynb`: susie finemapping results
- `tri_egene_biotype.ipynb`
- `tri_h2_supp.ipynb`
- `tri_specific.ipynb`
- `walker_fetal.ipynb`
- `Snakefile`
- `decon.smk`
- `paintor.smk`
### 4-2: cis-isoQTL
- `isoqtl_analysis.ipynb`
- `prep.ipynb`: sex and trimester specific QTL
- `Snakefile`: follows a similar pipeline as cis-eQTL, except that run grouped permutation as GTEx did
### 4-3: cis-sQTL
- `sqtl_analysis.ipynb`
- `e_iso_s.ipynb`
- `qvalue_pi0.ipynb`
- `check.ipynb`: check chunk size
- `Snakefile`
### 4-4: trans
- `gbat.ipynb`
- `Snakefile`
### 4-5: APEX
- `apex_analysis.ipynb`
- `Snakefile`
## 5: Integrative analysis
### 5-1: sLDSC 
- `ldsc_analysis.ipynb`
- `Snakefile`
- `pec.smk`
### 5-2: TWAS-FUSION
- `TWAS.ipynb`
- `LDREF.ipynb`
- `run_focus.sh`
- `Snakefile`
### 5-3: MESC
- `MESC.ipynb`
- `Snakefile`
-  `test.smk`
### 5-4: Colocalization (eCAVIAR)
- `eCAVIAR.ipynb`
- `GRIN2A.ipynb`
- `SP4_gviz.ipynb`
- `sqtlviztools.ipynb`
- `Visualizing_Loci_working.ipynb`
- `celltype.smk`
- `eqtl.smk`
- `isoqtl.smk`
- `mod_ieqtl.smk`
- `sex_tri.smk`
- `sqtl.smk`
- sashimi plot related code
## 6: Further Analyses
### 6-1: eGene/sGene Enrichment
- `fetal_only_egenes.ipynb`: biotype and cell type analysis for fetal-specific eGenes
- `trimester_egenes_sgenes.ipynb`: biotype and cell type analysis for trimester-specific e/sGenes
### 6-2: WGCNA
- `compare_module_enrichment.ipynb`: compare enrichment across networks and across correlated cell types
- `dashboard_generator.ipynb`: generate dashboards using ST6.xlsx
- `dashboards`: folder containing dashboards for each module
