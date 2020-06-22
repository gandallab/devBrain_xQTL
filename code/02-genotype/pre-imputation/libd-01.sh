#!/bin/bash
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/LIBD_1and2/geno/
#$ -j y
#$ -o ./job.out.impute.prep
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load plink/1.90b3.45
module load htslib/1.9

# apply plink filters
phase1and2_1M=./phase1and2/position.filtered/libd.phase1and2.1M.fixed.vcf.gz
phase1and2_h650=./phase1and2/position.filtered/libd.phase1and2.h650.fixed.vcf.gz
phase2only_1M=./phase2only/position.filtered/libd.phase2only.1M.fixed.final.vcf.gz
phase2only_h650=./phase2only/position.filtered/libd.phase2only.h650.fixed.vcf.gz

plink --vcf ${phase1and2_1M} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode vcf --out ./phase1and2/to_impute/libd.phase1and2.1M.filtered
bgzip ./phase1and2/to_impute/libd.phase1and2.1M.filtered.vcf
tabix -p vcf ./phase1and2/to_impute/libd.phase1and2.1M.filtered.vcf.gz

plink --vcf ${phase1and2_h650} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode vcf --out ./phase1and2/to_impute/libd.phase1and2.h650.filtered
bgzip ./phase1and2/to_impute/libd.phase1and2.h650.filtered.vcf
tabix -p vcf ./phase1and2/to_impute/libd.phase1and2.h650.filtered.vcf.gz

plink --vcf ${phase2only_1M} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode vcf --out ./phase2only/to_impute/libd.phase2only.1M.filtered
bgzip ./phase2only/to_impute/libd.phase2only.1M.filtered.vcf
tabix -p vcf ./phase2only/to_impute/libd.phase2only.1M.filtered.vcf.gz

plink --vcf ${phase2only_h650} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode vcf --out ./phase2only/to_impute/libd.phase2only.h650.filtered
bgzip ./phase2only/to_impute/libd.phase2only.h650.filtered.vcf
tabix -p vcf ./phase2only/to_impute/libd.phase2only.h650.filtered.vcf.gz
