#!/bin/bash
#$ -l h_data=32G,h_rt=8:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.post.impute.03.8
#$ -j y
#$ -m a

# merge datasets together, intersecting the high imputation quality variants (R2 > .8)
. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9

outfolder=/u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_8
mkdir -p ${outfolder}

#excluding LIBD h650 array and phase2only 1M
walker=./Walker2019_eQTL/geno/imputed/concat.R2.8.vcf.gz
obrien=./Obrien2018_eQTL/geno/imputed/concat.R2.8.vcf.gz
werling=./Werling2020_eQTL/geno/imputed/concat.R2.8.vcf.gz
libd_1and2_1M=./LIBD_1and2/geno/phase1and2/imputed/array_1M/concat.R2.8.vcf.gz
hdbr=./HDBR/geno/imputed/concat.R2.8.vcf.gz
# -c none, -c all, not much difference
# -c none is intersecting by chr:position:ref:alt
# -c all by chr:position
# -c all slightly outputing more intersecting variants, but not much
# also losing ref:alt information

bcftools isec ${walker} ${obrien} ${werling} ${libd_1and2_1M} ${hdbr} -c all -n=5 -p ${outfolder}
