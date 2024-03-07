#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=4:00:00,highp
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sLDSC/log/job.out.qqplot
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO
Rscript scripts/qqplot.R
