#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -j y
#$ -o ../log/job.out.ecaviar.sig
#$ -m a
#$ -t 1-5062

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

list=gwas_feature_fetal_scz_eqtl.txt
INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
IFS=$'\t' 
params=($INLINE)

locus=${params[0]}
gene=${params[1]}
gwas=${params[2]}
annot=${params[3]}

Rscript scripts/analyze.R \
            --gwas ${gwas} \
            --locus ${locus} \
            --annot fetal_${annot} \
            --gene ${gene}
