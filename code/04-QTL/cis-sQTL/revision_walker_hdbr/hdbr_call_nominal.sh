#!/bin/bash
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 5
#$ -cwd
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.call.nominal.hdbr
#$ -m a
#$ -j y

n_hcp=$1
dir=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/HDBR_nominal_${n_hcp}hcp/

#zcat ${dir}chunk*.txt.gz | gzip -c > ${dir}all.chunks.txt.gz

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO
Rscript /u/project/gandalm/cindywen/isoform_twas/eqtl_new/code/scripts/call_nominal.R \
    --outdir ${dir}
#gzip ${dir}all_assoc.txt
