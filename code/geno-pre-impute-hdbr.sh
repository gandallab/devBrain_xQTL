#!/bin/bash
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -wd /u/project/gandalm/cindywen/isoform_twas/HDBR_data/scripts/
#$ -j y
#$ -o ./../log/job.out.split.chr
#$ -m a
#$ -t 1-22

### split by chr to impute

. /u/local/Modules/default/init/modules.sh
module load vcftools/0.1.14
module load htslib

chr=${SGE_TASK_ID}
input_vcf=/u/project/gandalm/shared/GenomicDatasets/FetalBrain/HDBR/geno/HDBR_merged.vcf.gz
out_folder=/u/project/gandalm/shared/GenomicDatasets/FetalBrain/HDBR/geno/by_chr
mkdir -p ${out_folder}

vcftools --chr ${chr} --gzvcf ${input_vcf} --recode --recode-INFO-all --out ${out_folder}/HDBR.chr${chr}
# required sorted input by Michigan imputation server
vcf-sort ${out_folder}/HDBR.chr${chr}.recode.vcf | bgzip -c > ${out_folder}/HDBR.chr${chr}.recode.vcf.gz
