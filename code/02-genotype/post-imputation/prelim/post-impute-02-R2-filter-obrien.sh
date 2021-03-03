#!/bin/bash
#$ -l h_data=16G,h_rt=12:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9

#obrien
bcftools view -i 'R2>.8' ./Obrien2018_eQTL/geno/imputed/concat.all.vcf.gz -Oz > ./Obrien2018_eQTL/geno/imputed/concat.R2.8.vcf.gz
tabix -p vcf ./Obrien2018_eQTL/geno/imputed/concat.R2.8.vcf.gz
touch /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/post-imputation-02-obrien.done
