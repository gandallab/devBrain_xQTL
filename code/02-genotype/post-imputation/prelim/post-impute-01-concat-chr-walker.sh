#!/bin/bash
#$ -l h_data=16G,h_rt=4:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Walker2019_eQTL/geno/imputed/
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9

bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Oz > concat.all.vcf.gz
tabix -p vcf concat.all.vcf.gz
