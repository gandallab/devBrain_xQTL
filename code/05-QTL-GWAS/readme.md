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
## TWAS-FUSION
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
## MESC
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
## Colocalization (eCAVIAR)
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
