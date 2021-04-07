## Genotype 
### Pre-imputation
- Pre-imputation QC code, from starting data to files for imputation
- First apply PLINK filters, and then split by chromosome and sort
- Walker data is already filtered; we split the data by chromosome and imputed it
- For all the other datasets, we applied the same filters that the Walker data have had `--hwe 1e6 --maf 0.01 --mind 0.10 --geno 0.05`
- Note: for HDBR, we used `--mind 0.3`; for LIBD, we fixed strand flips by running an extra step of conform-gt, which automatically splits the data by chromosome
### Post-imputation
- Scripts in `prelim/`: inputs are imputed genotype files downloaded from Michigan Imputation server; concatenate by chromosomes, index, filter by R2, and take the **intersection** of high impute quality variants across datasets
- Note: except for Walker data, we applied R2>.3 filter during imputation; so here we only applied R2>.3 on Walker imputed data and intersected with the other datasets
- Snakefile rules description:
```
rules:
    walker_subj: remove 4 walker subjects that are not in rnaseq
    reheader: reheader to subject ID as in rnaseq; old_name new_name, same order as in original vcf file; or new_names for all samples, same order; can only output same format as input (BCF/VCF, bgzipped or not)
    convert_to_plink, remove_allele_in_id, remove_dup_update_name: remove duplicate position variants (multi-allelic), and convert topmed variant IDs to rsID, keeping indels. Reference: /u/project/gandalm/shared/GenomicDatasets/ABCD_r201_r1/impute/imputeABCD_July2020/results/TOPMED_postimputation-master
    keep_rsid_only: remove variants still in chr:pos:ref:alt ID. ~10% Topmed variants not mapped to rsID. Keeping these for analysis
    variant_qc: filter variants
    concat: concat chr plink files
    plink_to_vcf:
    crossmap: hg38 to hg19
    sort_tabix:
    vcf_to_plink: plink binary data will be used in ancestry PCA
    checkvcf: sanity check
    pca: merge data with 1000genomes and do PCA
    call ancestry in ancestry.ipynb, and do the following rules for eur, amr, afr: variant_qc_ancestry, concat_ancestry, plink_to_vcf_ancestry, crossmap_ancestry, sort_tabix_ancestry, vcf_to_plink_ancestry, checkvcf_ancestry, pca_ancestry
    rel_check: check plink pi_hat
    remove_rel and ancestry_remove_rel: remove relatives for QTL mapping
    remove_gene_expr_outlier and ancestry_remove_gene_expr_outlier: prepare genotype file for FastQTL
    remove_tx_expr_outlier and ancestry_remove_tx_expr_outlier: prepare genotype file for FastQTL
    add_chr: for STAR 2nd pass
```
- ancestry.ipynb: infer data ancestry, make plots
- IBD.ipynb: relatedness check
