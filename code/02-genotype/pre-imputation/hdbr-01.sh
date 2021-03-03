#!/bin/bash
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/HDBR/geno/
#$ -j y
#$ -o ./to_impute/job.out.impute.prep
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load plink/1.90b3.45
module load htslib/1.9

# filter samples; apply plink filters
list=/u/home/c/cindywen/project-gandalm/isoform_twas/FINAL/hdbrOnlySnpChip.txt
mkdir ./to_impute/
bcftools view -S ${list} HDBR_merged.vcf.gz > ./to_impute/HDBR_sampleFiltered.vcf
plink --vcf ./to_impute/HDBR_sampleFiltered.vcf --maf 0.01 --geno 0.05 --mind 0.3 --hwe 1e-6 --recode vcf --out ./to_impute/HDBR_filtered
bgzip ./to_impute/HDBR_filtered.vcf
tabix -p vcf ./to_impute/HDBR_filtered.vcf.gz
