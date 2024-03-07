#!/bin/bash
#$ -l h_data=4G,h_rt=4:00:00,highp
#$ -cwd
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.gtex.fastqtl.merge
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load python/3.9.6
module load R/4.1.0-BIO

cd /u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/

python3 /u/project/gandalm/cindywen/isoform_twas/sqtl_new/fastqtl-gtex/python/merge_chunks.py \
    list_chunks.txt \
    list_log.txt \
    gtex.new \
    -o .
