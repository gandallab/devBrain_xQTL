#!/bin/bash -l 
#$ -cwd
#$ -l h_data=40G,h_rt=6:00:00,highp
#$ -j y
#$ -o ./../log/job.out.snakemake
#$ -m a

source /u/local/apps/anaconda3/2019.03/bin/activate snakemake

snakemake \
    --snakefile Snakefile \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00,highp -pe shared {resources.num_cores} -o ./../log/job.out.pipeline" \
    --jobs 100 \
    --max-jobs-per-second 5 \
    --restart-times 0 \
    --latency-wait 15 \
    --default-resources mem_gb=4 time_min=240 num_cores=1 

