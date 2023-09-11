#!/bin/bash -l
#$ -cwd
#$ -l h_data=8G,h_rt=8:00:00,highp
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/eqtl_new/log/job.out.all_genes_batch_expr
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO
Rscript all_genes_batch_expr.R
