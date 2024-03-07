#!/bin/bash
#$ -l h_data=8G,h_rt=100:00:00,highp
#$ -cwd
#$ -o /u/project/gandalm/cindywen/isoform_twas/eqtl_new/log/job.out.run.paintor.new
#$ -m a
#$ -j y

gene=$1

. /u/local/Modules/default/init/modules.sh
module load gcc/7.5.0

cd /u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/${gene}_dir/
/u/project/gandalm/shared/apps/PAINTOR_V3.0/PAINTOR \
    -input input.file \
    -LDname LD1,LD2,LD3 \
    -Zhead ZSCORE.P1,ZSCORE.P2,ZSCORE.P3 \
    -annotations TF_binding_site_d,promoter_flanking_region_d,promoter_d,open_chromatin_region_d,enhancer_d,CTCF_binding_site_d \
    -enumerate 2 \
    -in . \
    -out . \
    -Gname Enrichment.Values.2causal \
    -Lname Log.BayesFactor.2causal \
    -RESname results.2causal \
    -set_seed 123
