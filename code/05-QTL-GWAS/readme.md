# QTL-GWAS
## sLDSC 
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
## TWAS-FUSION
- `TWAS_analysis.ipynb`
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
```
