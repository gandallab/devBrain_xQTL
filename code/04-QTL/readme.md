# QTL
## cis-eQTL
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
## cis-isoQTL
- `isoqtl_analysis.ipynb`
- `Snakefile`: follows a similar pipeline as cis-eQTL
## cis-sQTL
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
## trans
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