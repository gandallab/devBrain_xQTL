#!/bin/bash
#$ -l h_data=8G,h_rt=2:00:00,highp
#$ -cwd
#$ -o /u/project/gandalm/cindywen/isoform_twas/eqtl_new/log/job.out.paintor.fixLD
#$ -m a
#$ -j y

gene=$1

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

Rscript scripts/paintor_fixLD.R \
    --gene ${gene}
