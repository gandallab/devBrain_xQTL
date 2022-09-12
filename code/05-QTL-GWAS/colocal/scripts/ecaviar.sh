#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=12:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.ecaviar
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load gcc/10.2.0

locus=$1
feature=$2
trait=$3
level=$4


cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_${level}/${trait}/locus_${locus}/

/u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
    -o ${feature}_ecaviar \
    -l ${feature}_1kg_eur.ld \
    -z ${feature}_gwas_zscore.txt \
    -l ${feature}_qtl.ld \
    -z ${feature}_qtl_zscore.txt \
    -f 1 \
    -c 2 \
    -r 0.95
