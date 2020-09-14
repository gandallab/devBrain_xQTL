#!/bin/bash
#$ -l h_data=16G,h_rt=12:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9

#libd
bcftools view -i 'R2>.8' ./LIBD_1and2/geno/phase1and2/imputed/array_1M/concat.all.vcf.gz -Oz > ./LIBD_1and2/geno/phase1and2/imputed/array_1M/concat.R2.8.vcf.gz
bcftools view -i 'R2>.8' ./LIBD_1and2/geno/phase1and2/imputed/array_h650/concat.all.vcf.gz -Oz > ./LIBD_1and2/geno/phase1and2/imputed/array_h650/concat.R2.8.vcf.gz
bcftools view -i 'R2>.8' ./LIBD_1and2/geno/phase2only/imputed/array_1M/concat.all.vcf.gz -Oz > ./LIBD_1and2/geno/phase2only/imputed/array_1M/concat.R2.8.vcf.gz
bcftools view -i 'R2>.8' ./LIBD_1and2/geno/phase2only/imputed/array_h650/concat.all.vcf.gz -Oz > ./LIBD_1and2/geno/phase2only/imputed/array_h650/concat.R2.8.vcf.gz
#tabix -p vcf
touch /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/post-imputation-02-libd.done
