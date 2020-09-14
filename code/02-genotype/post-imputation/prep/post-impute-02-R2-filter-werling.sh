#!/bin/bash
#$ -l h_data=8G,h_rt=4:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9

#werling
bcftools view -i 'R2>.8' ./Werling2020_eQTL/geno/imputed/concat.all.vcf.gz -Oz > ./Werling2020_eQTL/geno/imputed/concat.R2.8.vcf.gz
tabix -p vcf ./Werling2020_eQTL/geno/imputed/concat.R2.8.vcf.gz
touch /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/post-imputation-02-werling.done
