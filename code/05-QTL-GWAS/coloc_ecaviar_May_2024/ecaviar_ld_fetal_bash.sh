#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -j y
#$ -o ../log/job.out.ecaviar.ld.bash
#$ -m a
#$ -t 1-55011

. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624

list=gwas_feature_fetal_all.txt
INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
IFS=$'\t'
params=($INLINE)

locus=${params[0]}
gene=${params[1]}
gwas=${params[2]}
annot=${params[3]}
chr=${params[4]}


cd ../out_1MB/${gwas}/locus${locus}/

plink --bfile /u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
      --r \
      --matrix \
      --extract ${gene}_fetal_${annot}_snps.txt \
      --out ${gene}_fetal_${annot}_gwas

plink --bfile /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/test \
      --r \
      --matrix \
      --extract ${gene}_fetal_${annot}_snps.txt \
      --out ${gene}_fetal_${annot}

