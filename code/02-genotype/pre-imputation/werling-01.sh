#!/bin/bash
#$ -l h_data=40G,h_rt=4:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Werling2020_eQTL/geno/
#$ -j y
#$ -o ./to_impute/job.out.impute.prep
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load plink/1.90b3.45
module load htslib/1.9

# filter samples; apply plink filters
# had problems with WGS variant ID when working in VCF format, now first work in PLINK format, then convert to VCF when filtering is finished
list=/u/home/c/cindywen/project-gandalm/isoform_twas/FINAL/werling.geno.exclude.samples.txt
mkdir ./to_impute/
bcftools view -S ^${list} PEC_HSB_Yale-UCSF_WGS.vcf.gz > ./to_impute/werling_sampleFiltered.vcf
plink --vcf ./to_impute/werling_sampleFiltered.vcf --allow-extra-chr --make-bed --out ./to_impute/werling_sampleFiltered
plink --bfile ./to_impute/werling_sampleFiltered --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --allow-extra-chr --make-bed --out ./to_impute/werling_filtered
plink --bfile werling_filtered --recode vcf-iid --allow-extra-chr --out werling_filtered
bgzip ./to_impute/werling_filtered.vcf
tabix -p vcf ./to_impute/werling_filtered.vcf.gz
