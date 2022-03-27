# QTL
## cis-eQTL
- `metadata.ipynb`: plot data age, sex, infer NA sex, etc.
- `eqtl_analysis.ipynb`: identify optimal #HCP in covariates, gene expression PCA, dTSS, etc.
- `susie.ipynb`: susie finemapping results
- `decon.ipynb`: cell type specific and interacting analysis
- `cell_specific.ipynb`: (OUTDATED) cell type/group specific and interaction results
- `func_enrich.ipynb`: functional enrichment analysis of QTL
- `paintor.ipynb`: PAINTOR multi-ethnic fine-mapping 
- `fetal_adult.ipynb`
- `paintor.smk`
- `pLI.ipynb`
- `sex_specific.ipynb`
-  `SV.ipynb`
- `tri_specific.ipynb`
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
## cis-isoQTL
- `isoqtl_analysis.ipynb`
- `prep.ipynb`: sex and trimester specific QTL
- `Snakefile`: follows a similar pipeline as cis-eQTL, except that run grouped permutation as GTEx
## cis-sQTL
- `sqtl_analysis.ipynb`
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
## trans
- `gbat.ipynb`
- `Snakefile`
```
rules:
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
## APEX
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
