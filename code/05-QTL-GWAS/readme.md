# QTL-GWAS
## sLDSC 
- `ldsc_analysis.ipynb`: make partitioned h2 plots
- `Snakefile`
```
rules:
    - make_coord_file
Top eQTL
    - make_set_file_mixed_top_eqtl
    - make_annot_mixed_top_eqtl
    - partition_h2_mixed_top_eqtl
Top isoQTL
    (as above)
SuSiE variants
    (as above)
Top sQTL
    (as above)
```
## TWAS-FUSION
- `TWAS.ipynb`
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
## MESC
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
```
