#!/bin/bash
#$ -l h_data=4G,h_rt=4:00:00,highp
#$ -pe shared 4
#$ -cwd
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.qvalue.s.s
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO
Rscript scripts/qvalue_pi0_s_s.R
