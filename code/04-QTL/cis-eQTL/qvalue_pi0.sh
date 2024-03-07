#!/bin/bash -l
#$ -cwd
#$ -l h_data=12G,h_rt=8:00:00,highp,highmem
#$ -pe shared 15
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/eqtl_new/log/job.out.qvalue.pi0.s
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO
Rscript scripts/qvalue_pi0_s.R
