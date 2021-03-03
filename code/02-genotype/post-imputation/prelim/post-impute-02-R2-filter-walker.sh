#!/bin/bash
#$ -l h_data=20G,h_rt=12:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9

# walker
#bcftools view -i 'R2>.8' ./Walker2019_eQTL/geno/imputed/concat.all.vcf.gz -Oz > ./Walker2019_eQTL/geno/imputed/concat.R2.8.vcf.gz
bcftools view -i 'R2>.3' ./Walker2019_eQTL/geno/imputed/concat.all.vcf.gz -Oz > ./Walker2019_eQTL/geno/imputed/concat.R2.3.vcf.gz
tabix -p vcf ./Walker2019_eQTL/geno/imputed/concat.R2.3.vcf.gz
touch /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/post-imputation-02-walker.3.done
