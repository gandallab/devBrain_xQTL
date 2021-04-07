## QTL mapping 
- metadata.ipynb in cis-eQTL: plot data basics, age, sex, infer NA sex
- eqtl/iso_analysis.ipynb in cis-e/isoQTL: identify optimal #HCP in covariates, gene expression PCA, dTSS, etc.
- cis-eQTL and cis-isoQTL follow a similar pipeline. See Snakefile
```
rules:
- expr_prep: expression filter, VST normalization, outlier detection, ComBat, create bed file
- ancestry_expr: subset ancestries
- cov: generate covariates files
- ancestry_cov: 
- fastqtl_nominal: calculate nominal association for all feature/cis-SNP pairs
- ancestry_fastqtl_nominal: 
- call_nominal: identify significant nominal associations
- ancestry_call_nominal:
- fastqtl_perm: calculate beta approximated permuted p-value, using optimal #HCP identified from nominal mode
- ancestry_fastqtl_perm:
- call_perm: call eGene/isoTx, top QTL
- ancestry_call_perm:
```
- apex. See Snakefile and apex_analysis.ipynb
```
rules:
- From split_chr_prep_vcf to pca_plots_by_group: refer to ABCD_GWAS Snakefile, running GENESIS for ancestry PCA, and ancestry-aware kinship estimation
- make_apex_kin_mat: convert pcrelate RData to kinship sparse matrix for apex, refer to apex documentation for format details
- factor: generating covarites file from known factors and expression factor analysis implemented in apex. (Note: in current version of apex, if a kin matirix is included, eFA will have to be modeled as fixed effects. If no kin/grm included, eFA still can only be modeled as fixed effects. To model inferred factors as random effects, use --epcs $num_factor and --cov $file_with_known_factor_only)
- cis_lmm_kin: cis-eQTL mapping, mixed ancestry data, lmm with kin as random effects
- cis_lmm_dtss_kin: specify an alpha for dtss weight
- cis_ols:
- make_apex_grm_mat:
- cis_lmm_grm:
- cis_lmm_dtss_grm:
```
- cis-sQTL. See bash scripts and Snakefile
```
rules:
- index: generate STAR index with Gencode v33, with annotation GTF as recommended
- See step 1-5 for running STAR, 1st and 2nd pass, WASP filter; Leafcutter bam2junc
- cluster: note some output files are not specified here
- remove_chr:
- write_chr_blacklist
- pheno_prep: note some output files are not specified here
- bgzip_tabix
- concat
- pheno_process
- cov
- fastqtl_nominal
```
  
