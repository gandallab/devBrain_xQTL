# Fetal brain mega QTL
## Pipeline
#### 1: Munging data
#### 2: Genotype
##### Pre-imputation
- Run plinkQC on data as sanity check
- First apply PLINK filters, then split by chromosome and sort
- Walker data is already filtered; split by chromosome and imputed
- For all the other datasets, we applied the same filters that the Walker data used `--hwe 1e6 --maf 0.01 --mind 0.10 --geno 0.05`
- Note: for HDBR, we used `--mind 0.3`; for LIBD, we fixed strand flips by running an extra step of [conform-gt](https://faculty.washington.edu/browning/conform-gt.html), which automatically splits the data by chromosome
##### Post-imputation
- Scripts in `prelim/`: inputs are imputed genotype files downloaded from Michigan Imputation Server; concatenate by chromosomes, index, filter by R2, and take the **intersection** of high impute quality variants across datasets
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
```
#### 3: RNA-seq
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
-   Sample swap check: 
    + [VerifyBamID](https://genome.sph.umich.edu/wiki/VerifyBamID) (slow. Use `--smID` to add subject ID to BAM sequence file)
    + `check.ipynb`: called SNP from BAM, merged with imputed genotype
#### 4: xQTL
##### cis-eQTL
- `metadata.ipynb`: plot data demographics, age, sex, infer NA sex, etc.
- `eqtl_analysis.ipynb`: identify optimal #HCP in covariates, gene expression PCA, dTSS, etc.
- `susie_analysis.ipynb`: susie finemapping results
- `cell_specific_analysis.ipynb`: cell type/group specific and interaction results
- `functional_enrichment.ipynb`: functional enrichment analysis of QTL
- `PAINTOR.ipynb`: PAINTOR multi-ethnic fine-mapping 
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
    ...as above
cell type/group interaction
    - make_decon_dosage
    - snps_to_test
    - fix_decon_dosage
    - run_decon_qtl
cell type/group specific
    - bgzip_tabix
    - cg_cov
    - cg_fastqtl_nominal
    - cg_call_nominal
Ancestry eQTL effect size
    - make_effect_size_scatter_eur_amr, eur_afr, afr_amr
Torus: functional enrichment
    - make_annot
    - merge_annot
    - fastqtl_calculate_sebeta
    - run_torus
PAINTOR: multi-ethnic fine-mapping
    - make_eur_coord: as in sLDSC for ALL, make variant coord file for EUR. ALL does not cover all shared variants between EUR, AMR, AFR
    - make_eur_annot: need variant annot for shared variants between EUR, AMR, AFR
    - merge_eur_annot:
```
##### cis-isoQTL
- `isoqtl_analysis.ipynb`
- `Snakefile`: follows a similar pipeline as cis-eQTL
##### cis-sQTL
- `sqtl_analysis.ipynb`
- `Snakefile`
```
rules:
STAR
    - index: generate STAR index with Gencode v33, with annotation GTF as recommended
    - See step 1-5 for running STAR, 1st and 2nd pass, WASP filter; Leafcutter bam2junc
Leafcutter combined ancestry
    - cluster: note some output files are not specified here
    - remove_chr: 
    - write_chr_blacklist
    - pheno_prep: note some output files are not specified here
    - bgzip_tabix
    - concat
    - pheno_process
    - cov
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
```
##### trans
- `gbat.ipynb`
- `Snakefile`
```
rules:
- filter (input BAM already passed the filters)
- gtf_filter_mappability: remove exons in  GTF that overlap with ENCODE low mappability regions
- merge_gtf: merge chr
- count: run featureCounts
- rdata: generate rdata with all sample featureCounts, CPM filter, quantile normalize TPM (gene and sample), standardize
- make_geno
- cvBLUP: predict gene expression
- move
- pearsonR: R2 of predicted and observed gene expression
- cor: use cal_cor_covar.txt to include covariates and to combine with supervised SVA; not using script cal_cor.R
- pval
- qval
- sig
```
##### APEX
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
#### 5: Integrative
##### sLDSC 
- `ldsc_analysis.ipynb`: make partitioned h2 plots
- `Snakefile`
```
rules:
- make_coord_file
- make_set_file_mixed_top_eqtl
- make_annot_mixed_top_eqtl
- partition_h2_mixed_top_eqtl
- make_set_file_mixed_top_isoqtl
- make_annot_mixed_top_isoqtl
- partition_h2_mixed_top_isoqtl
```
##### TWAS-FUSION
- `TWAS_analysis.ipynb`
- `Snakefile`: tested running with and without rank normalization of gene expression
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
```
------
## Data and results
#### Quantifications
##### Gene (`working/gene/`)
- [x] Estimated counts (`"countsFromAbundance="lengthScaledTPM"`): `gene.noVersion.scaled.counts.tsv`
- [x] TPM: `gene.noVersion.TPM.tsv`
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: `gene.counts.scaled.normalized.bed.gz`
##### Isoform (`working/transcript/`)
- [x] Estimated counts (`"countsFromAbundance="lengthScaledTPM"`): `tx.counts.scaled.tsv`
- [x] TPM: `tx.TPM.tsv`
- [x] Filtered for expression, normalized, variance-stabilized, transformed, ComBat: `tx.counts.scaled.normalized.bed.gz`
##### Splicing (`working/splicing/`)
- [x] Intron excision ratios: filtered, normalized, ComBat: `leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.bed.gz`

#### Genotype (`working/geno/`)
***Please note some TOPMED variants are not mapped to rsID. They are still in hg38 coordinates as ID (chr:pos). Their coordinates are hg19 in the final genotype file.***
- AMR, AFR, EUR genotype plink files used for FUSION
  
#### Covariates (`working/covariates/`)

#### cis-eQTL (`output/cis_eQTL_nominal_permutation_conditional_finemapping/`)
* Combined (`shared_90HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_info.txt`. Filter for significant: `qval < .05`
- [ ] Permutation pass: full list of eQTL with p-value passing nominal threshold: 
- [x] Conditional pass, top eQTL per rank: `conditional_top_variants.txt`. To remove variants with backward P-value that is not below the threshold of this feature: `filter(V20 == 1)`
- [x] SuSiE finemapping: all variants in non-low purity CS: `mixed_ciseqtl_90hcp_perm_purity_filtered.txt.gz`
* EUR (`EUR_50HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_info.txt`. Filter for significant: `qval < .05`
* AMR (`AMR_15HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_info.txt`. Filter for significant: `qval < .05`
* AFR (`AFR_25HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_info.txt`. Filter for significant: `qval < .05`

#### cis-isoQTL (`output/cis_isoQTL_nominal_permutation_conditional_finemapping/`)
* Combined (`shared_70HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`
- [ ] Permutation pass: full list of isoQTL with p-value passing nominal threshold: 
- [x] Conditional pass, top isoQTL per rank: `conditional_top_variants.txt`. To remove variants with backward P-value that is not below the threshold of this feature: `filter(V20 == 1)`
- [x] SuSiE finemapping: all variants in non-low purity CS: `mixed_cisisoqtl_70hcp_perm_purity_filtered.txt.gz`
* EUR (`EUR_60HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`
* AMR (`AMR_15HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`
* AFR (`AFR_20HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`

#### cis-sQTL (`output/cis_sQTL_nominal_permutation_conditional_finemapping/`)
* Combined (`shared_35HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`
- [ ] Permutation pass: full list of sQTL with p-value passing nominal threshold: 
- [ ] Conditional pass: top sQTL per rank: 
- [ ] SuSiE finemapping; all variants in non-low purity CS: 
* EUR (`EUR_25HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`
* AMR (`AMR_10HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`
* AFR (`AFR_12HCP/`)
- [x] Nominal pass: `all_assoc_nominal.txt.gz`.  Filter for significant: `fdr <= .05`
- [x] Permutation pass: `all_assoc_perm_gene_info.txt`. Filter for significant: `qval < .05`

#### trans-eQTL (`output/trans_eQTL/`)

#### FUSION-TWAS weights (`output/TWAS_wights/`)
- EUR, AMR, AFR; with or without rank normalization