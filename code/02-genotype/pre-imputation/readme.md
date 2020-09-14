#### Pre-imputation QC
- Pre-imputation QC code, from starting data to files for imputation
- First apply PLINK filters, and then split by chromosome and sort
- Walker data is already filtered; we split the data by chromosome and imputed it
- For all the other datasets, we applied the same filters that the Walker data have had `--hwe 1e6 --maf 0.01 --mind 0.10 --geno 0.05`
- Note: for HDBR, we used `--mind 0.3`; for LIBD, we fixed strand flips by running an extra step of conform-gt, which automatically splits the data by chromosome
