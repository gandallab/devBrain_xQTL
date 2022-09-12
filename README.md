# Fetal brain mega QTL
***Note: please see snakefiles for exact rules and descriptions***

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
- `Snakefile`:
```
rules:
- walker_subj: remove 4 walker subjects that are not in rnaseq
- reheader: reheader to subject ID as in rnaseq; old_name new_name, same order as in original vcf file; or new_names for all samples, same order; can only output same format as input (BCF/VCF, bgzipped or not)
- convert_to_plink, remove_allele_in_id, remove_dup_update_name: remove duplicate position variants (multi-allelic), and convert topmed variant IDs to rsID, keeping indels. Reference: /u/project/gandalm/shared/GenomicDatasets/ABCD_r201_r1/impute/imputeABCD_July2020/results/TOPMED_postimputation-master
- keep_rsid_only: remove variants still in chr:pos:ref:alt ID. ~10% Topmed variants not mapped to rsID. Keeping these for analysis
- variant_qc: filter variants
- concat: concat chr plink files
- plink_to_vcf
- crossmap: hg38 to hg19
- sort_tabix
- vcf_to_plink: plink binary data will be used in ancestry PCA
- checkvcf: sanity check
- pca: merge data with 1000genomes and do PCA
- call ancestry in ancestry.ipynb, and do the following rules for eur, amr, afr: variant_qc_ancestry, concat_ancestry, plink_to_vcf_ancestry, crossmap_ancestry, sort_tabix_ancestry, vcf_to_plink_ancestry, checkvcf_ancestry, pca_ancestry
- rel_check: check plink pi_hat
- remove_rel and ancestry_remove_rel: remove relatives for QTL mapping
- remove_gene_expr_outlier and ancestry_remove_gene_expr_outlier: prepare genotype file for FastQTL
- remove_tx_expr_outlier and ancestry_remove_tx_expr_outlier: prepare genotype file for FastQTL
- add_chr: for STAR 2nd pass
- get_maf: for susie ancestry analysis, get MAF of variants
```
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
- `combat-seq`
- `decon.ipynb`: cell type specific and interacting analysis
- `eqtl_analysis.ipynb`: identify optimal #HCP in covariates, gene expression PCA, dTSS, etc.
- `fetal_adult.ipynb`
- `func_enrich.ipynb`: functional enrichment analysis of QTL
- `locuszoom.ipynb`
- `metadata.ipynb`: plot data age, sex, infer NA sex, etc.
- `paintor.ipynb`: PAINTOR multi-ethnic fine-mapping 
- `pLI.ipynb`
- `sex_specific.ipynb`
- `susie.ipynb`: susie finemapping results
- `tri_specific.ipynb`
- `walker_fetal`

- `Snakefile`
```
prep
    - expr_prep: expression filter, VST normalization, outlier detection, ComBat, create bed file
    - ancestry_expr: subset ancestries
    - cov: generate covariates files
    - ancestry_cov
nominal pass
    - fastqtl_nominal
    - ancestry_fastqtl_nominal
    - merge_nominal
    - ancestry_merge_nominal
    - call_nominal: identify significant nominal associations
    - ancestry_call_nominal
permutation pass
    - fastqtl_perm
    - ancestry_fastqtl_perm
    - merge_perm
    - ancestry_merge_perm
    - call_perm
    - ancestry_call_perm
conditional pass
    - permutations_all_threshold: compute a npval threshold for all expressed genes
    - conditional
    - get_top_variants
susie finemapping
    - vcf_to_dosage
    - make_susie_meta
    - make_susie_expr
    - run_susie
    - merge_susie
    - sort_susie
susie finemapping ancsetry
    (as above)
    - locuszoom: also see susie_analysis.ipynb
Ancestry eQTL effect size
    - make_effect_size_scatter_eur_amr, eur_afr, afr_amr: all nominal associations
    - scatter_eur_afr, eur_amr, afr_amr: x-axis pop1 eGene-eQTL, y-axis sig or non-sig
Torus: functional enrichment
    - make_annot_ensembl: use ensembl regulatory build to annotate our variants
    - run_vep: ensembl variant effect predictor
    - select_consequences: select rsID and consequence/annotation columns
    - add_vep: add VEP to annotation
    - merge_annot_ensembl_vep: merge chr
    - (fastqtl_calculate_sebeta): too slow for all nominal association, run GTEx's FastQTL to get sebeta
    - (merge_sebeta)
    - gtex_fastqtl
    - run_torus
PAINTOR: multi-ethnic fine-mapping
    - make_eur_coord: as in sLDSC for ALL, make variant coord file for EUR. ALL does not cover all shared variants between EUR, AMR, AFR
    - make_eur_annot: need variant annot for shared variants between EUR, AMR, AFR
    - merge_eur_annot:
```
- `decon.smk`
```
rules:
cell type specific
    - ct_cov
    - ct_fastqtl_nominal
    - ct_merge_nominal
    - ct_call_nominal
    - ct_fastqtl_perm
    - ct_merge_perm
    - ct_call_perm
cell type/group interaction
    - make_decon_dosage
    - snps_to_test
    - fix_decon_dosage
    - run_decon_qtl
```
- `paintor.smk`
```
rules:
```
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
```
rules:
STAR
    - index: generate STAR index with Gencode v33, with annotation GTF as recommended
    - See step 1-5 for running STAR, 1st and 2nd pass, WASP filter; Leafcutter bam2junc
Leafcutter
    - cluster: note some output files are not specified here
    - remove_chr
    - write_chr_blacklist
    - pheno_prep: note some output files are not specified here
    - bgzip_tabix
    - concat
    - pheno_process
    - cov
ALL
    - fastqtl_nominal
    - merge_nominal
    - call_nominal
    - fastqtl_perm
    - merge_perm
    - call_perm
Separate ancestry
    - ancestry_pheno
    - ancestry_cov
    - ancestry_fastqtl_nominal
    - ancestry_merge_nominal
    - ancestry_call_nominal
    - ancestry_fastqtl_perm
    - ancestry_merge_perm
    - ancestry_call_perm
Intron annotation
    - prepare_leafviz_annot
    - annotate_intron
    - summarise_annot
GTEx way
    - prepare_leafviz_annot_exons: add gene_id in all_exons file
    - map_clusters_to_genes
    - make_grp_and_bed_file
    - fastqtl_grp_perm
    - gtex_grp_perm_write_chunk_log_list
    - gtex_grp_perm_merge_chunk
Conditional
    - permutations_all_threshold: compute a npval threshold for all tested introns
    - conditional: run QTLtools conditional pass
    - get_top_variants: get top variants per rank
susie finemapping
    - vcf_to_dosage
    - make_susie_meta
    - make_susie_expr
    - run_susie
    - merge_susie
    - sort_susie
Torus: functional enrichment
    - (fastqtl_calculate_sebeta)
    - (merge_beta)
    - gtex_fastqtl
    - gtex_write_chunk_log_list
    - gtex_merge_chunk
    - make_torus_input
    - run_torus
```
### 4-4: trans
- `gbat.ipynb`
- `Snakefile`
```
rules:
- filter (input BAM already passed the filters)
- gtf_filter_mappability: remove exons in  GTF that overlap with ENCODE low mappability regions (map score < 1)
- merge_gtf: merge chr
- count: run featureCounts
- rdata: generate rdata with all sample featureCounts, quantile normalize TPM (gene and sample), standardize, filter genes by expression
- make_geno
- cvBLUP: predict gene expression
- move
- pearsonR: R2 of predicted and observed gene expression
- cor: use cal_cor_covar.txt to include covariates and to combine with supervised SVA; not using script cal_cor.R
- pval
- qval
- sig
```
### 4-5: APEX
- `apex_analysis.ipynb`
- `Snakefile`:
```
rules:
GENESIS
    - From split_chr_prep_vcf to pca_plots_by_group: refer to ABCD_GWAS Snakefile, running GENESIS for ancestry PCA, and ancestry-aware kinship estimation
APEX OLS mapping
    - factor: generating covarites file from known factors and expression factor analysis implemented in apex. (Note: in current version of apex, if a kin matirix is included, eFA will have to be modeled as fixed effects. If no kin/grm included, eFA still can only be modeled as fixed effects. To model inferred factors as random effects, use --epcs $num_factor and --cov $file_with_known_factor_only)
    - cis_ols
APEX LMM mapping with kin
    - make_apex_kin_mat: convert pcrelate RData to kinship sparse matrix for apex, refer to apex documentation for format details
    - cis_lmm_kin: cis-eQTL mapping, mixed ancestry data, lmm with kin as random effects
    - cis_lmm_dtss_kin: specify an alpha for dtss weight
APEX LMM mapping with grm
    - make_apex_grm_mat 
    - cis_lmm_grm
    - cis_lmm_dtss_grm
APEX trans
    - trans_ols
```
## 5: Integrative analysis
### 5-1: sLDSC 
- `ldsc_analysis.ipynb`
- `Snakefile`
```
rules:
    - make_coord_file
Top eQTL
    - make_set_file_mixed_top_eqtl
    - make_annot_mixed_top_eqtl
    - partition_h2_mixed_top_eqtl
Top isoQTL
    (as above, note: use top QTL per gene and GTEx grouped permutation results)
SuSiE variants
    (as above)
Top sQTL
    (as above, note: use top QTL per gene and GTEx grouped permutation results)
eQTL maxCPP
    - make_annot_eqtl_maxCPP: make annotation files with maxCPP as continuous annotation
    - ldsc_eqtl_maxCPP
    - partition_h2_eqtl_maxCPP
    - partition_h2_gtex_brain_cortext_maxCPP
isoQTL_maxCPP
sQTL_maxCPP
CT specific maxCPP
jointly e, iso, sQTL maxCPP
trimester specific maxCPP
```
### 5-2: TWAS-FUSION
- `TWAS.ipynb`
- `LDREF.ipynb`
- `run_focus.sh`
- `Snakefile`
```
rules:
- make_plink
- compute_weights (rn)
- concat_hsq (rn)
- write_wgtlist (rn)
- summary_wgt (rn)
- make_pos_file (rn)
- assoc (rn)
- concat_extract (rn)
- chr_sig_rn
- pos_process_rn
```
### 5-3: MESC
- `MESC.ipynb`
- `Snakefile`
```
rules:
Gene set analysis (using individual level expression and genotype)
    - plink_covar: generate covar file in plink format
    - expr_rel: remove relatatives in expression
    - geno_split: generate chr genotype files
    - effect_egene: estimate eQTL effect sizes
    - gene_set: generate egene set file
    - score_egene: estimate gene set expression scores
    - h2med_egene: estimate h2med
Overall gene expression analysis
    - score_all_gene: estimate overall gene expression score
    - h2med_all_gene: estimate h2med
Overall isoform expression
    - plink_covar_iso
    - expr_rel_iso
    - score_all_iso
    - h2med_all_iso
Overall splicing
    - plink_covar_intron
    - expr_rel_intron
    - score_all_intron
    - h2med_all_intron
Trimester expression genes
    - score_all_gene_tri(test): used (1) trimester specific covaraites and (2/test) EUR corrected file separated into trimesters
    - h2_med_all_gene_tri(test)
Sex specific
```
### 5-4: Colocalization (eCAVIAR)
- `eCAVIAR.ipynb`
- `GRIN2A.ipynb`
- `locuszoom.ipynb`
- `sqtlviztools.ipynb`
- `Visualizing_Loci_working.ipynb`
- `celltype.smk`
- `eqtl.smk`
- `isoqtl.smk`
- `sex_tri.smk`
- `sqtl.smk`