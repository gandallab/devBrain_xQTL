#!/bin/bash -l
#$ -cwd
#$ -l h_data=48G,h_rt=4:00:00,highp
#$ -j y
#$ -o ../log/job.out.ecaviar.locus.bash
#$ -m a
#$ -t 1-396

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

list=gwas_fetal_torun.txt
INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
IFS=$'\t'
params=($INLINE)

gwas=${params[0]}
file=${params[1]}
locus=${params[2]}
annot=${params[3]}
mkdir -p ../out_1MB/${gwas}/locus${locus}
Rscript scripts/locus_feature.R \
            --gwas ${gwas} \
            --locus ${locus} \
            --annot ${annot} \
            --gwas_file ${file}

