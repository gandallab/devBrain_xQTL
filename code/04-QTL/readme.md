# QTL mapping 
## cis-eQTL
- metadata.ipynb: plot data demographics, age, sex, infer NA sex, etc.
- eqtl_analysis.ipynb: identify optimal #HCP in covariates, gene expression PCA, dTSS, etc.
- susie_analysis: susie finemapping results
- Snakefile
```
rules:
prep:
    expr_prep: expression filter, VST normalization, outlier detection, ComBat, create bed file
    ancestry_expr: subset ancestries
    cov: generate covariates files
    ancestry_cov
nominal pass:
    fastqtl_nominal
    ancestry_fastqtl_nominal
    call_nominal: identify significant nominal associations
    ancestry_call_nominal
permutation pass:
    fastqtl_perm
    ancestry_fastqtl_perm
    call_perm
    ancestry_call_perm
conditional pass:
    permutations_all_threshold: compute a npval threshold for all expressed genes
    conditional
    get_top_variants
susie finemapping:
    vcf_to_dosage
    make_susie_meta
    make_susie_expr
    run_susie
    merge_susie
    sort_susie
```
## cis-isoQTL
- isoqtl_analysis.ipynb
- Snakefile follows a similar pipeline as cis-eQTL
## cis-sQTL
- sqtl_analysis.ipynb
- Snakefile
```
rules:
- index: generate STAR index with Gencode v33, with annotation GTF as recommended
- See bash scripts step 1-5 for running STAR, 1st and 2nd pass, WASP filter; Leafcutter bam2junc
- cluster
- remove_chr
- write_chr_blacklist
- pheno_prep
- bgzip_tabix
- concat
- pheno_process
- cov
- fastqtl_nominal
```
## trans
## APEX
- apex_analysis.ipynb
- Snakefile
```
rules:
- From split_chr_prep_vcf to pca_plots_by_group: run GENESIS for ancestry PCA, and ancestry-aware kinship estimation
- make_apex_kin_mat: convert pcrelate RData to kinship sparse matrix for apex, refer to apex documentation for format details
- factor: generating covarites file from known factors and expression factor analysis implemented in apex. (Note: in current version of apex, if a kin matirix is included, eFA will have to be modeled as fixed effects. If no kin/grm included, eFA still can only be modeled as fixed effects. To model inferred factors as random effects, use --epcs $num_factor and --cov $file_with_known_factor_only)
- cis_lmm_kin: cis-eQTL mapping, mixed ancestry data, lmm with kin as random effects
- cis_lmm_dtss_kin: specify an alpha for dtss weight
- cis_ols
- make_apex_grm_mat
- cis_lmm_grm
- cis_lmm_dtss_grm
```

  
