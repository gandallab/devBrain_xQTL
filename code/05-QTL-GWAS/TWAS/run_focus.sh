#!/bin/bash -l
#$ -cwd
#$ -l h_data=12G,h_rt=4:00:00,highp
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/TWAS/log/job.out.focus
#$ -m a
#$ -t 1-21

. /u/local/Modules/default/init/modules.sh
module load python/3.9.6

# 0. See GWAS sum stats folder

# Had error running the focus in shared folder
# [2022-08-09 20:16:12 - ERROR] import_fusion() missing 2 required positional arguments: 'from_gencode' and 'session'

# pip3 install mygene --user
# module load gcc/11.3.0
# pip3 install rpy2 --user

# need to add DIR column to pos file

# 1. To generate FOCUS db
# ~/.local/bin/focus import ~/project-gandalm/isoform_twas/TWAS/results/gene_all_LDREF_rn/gene_all_FOCUS.pos fusion \
#     --tissue brain \
#     --name all \
#     --assay rnaseq \
#     --output fusion

# 2. Fine-mapping
~/.local/bin/focus finemap \
    /u/project/gandalm/shared/GWAS/SCZ.PGC3.2021/wave3_v3/PGC3_SCZ_wave3.european.autosome.public.v3.FOCUS_cleaned.sumstats.gz \
    /u/project/gandalm/shared/TWAS/code/LDREF/1000G.EUR.${SGE_TASK_ID} \
    ~/project-gandalm/isoform_twas/TWAS/results/gene_all_LDREF_rn/FOCUS/fusion.db \
    --chr ${SGE_TASK_ID} \
    --out ~/project-gandalm/isoform_twas/TWAS/results/gene_all_LDREF_rn/FOCUS/all_${SGE_TASK_ID} \
    --plot \
    --locations 37:EUR
