#!/bin/bash -l
#$ -cwd
#$ -l h_data=50G,h_rt=8:00:00,highp
#$ -j y
#$ -o ./../log/job.out.snakemake
#$ -m a

source /u/local/apps/anaconda3/2019.03/bin/activate snakemake

snakemake \
    --snakefile pipeline.smk \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00,highp -pe shared {resources.num_cores} -o ./../log/job.out.pipeline" \
    --jobs 100 \
    --printshellcmds \
    --max-jobs-per-second 1 \
    --restart-times 2 \
    --latency-wait 5 \
    --default-resources mem_gb=4 time_min=240 num_cores=1
