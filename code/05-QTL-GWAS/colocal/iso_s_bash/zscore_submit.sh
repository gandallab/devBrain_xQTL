#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=2:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.zscore.iso
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

Rscript /u/project/gandalm/cindywen/isoform_twas/colocal/code/scripts/zscore_iso_s.R \
    --locus $1 \
    --gene $2 \
    --trait $3 \
    --level $4
