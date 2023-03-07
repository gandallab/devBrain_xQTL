#!/bin/bash -l 
#$ -cwd
#$ -l h_data=18G,h_rt=12:00:00,highp
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.snakemake.ct
#$ -m a

# source /u/local/apps/anaconda3/2019.03/bin/activate snakemake
module load anaconda3
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
    --snakefile celltype.smk \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00,highp -pe shared {resources.num_cores} -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.pipeline.ct" \
    --jobs 100 \
    --max-jobs-per-second 10 \
    -T 0 \
    -w 60 \
    --default-resources mem_gb=4 time_min=240 num_cores=1 
