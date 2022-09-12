#!/bin/bash -l 
#$ -cwd
#$ -l h_data=12G,h_rt=6:00:00
#$ -j y
#$ -o ./../log/job.out.snakemake
#$ -m a

module load anaconda3
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
    --snakefile Snakefile \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00 -pe shared {resources.num_cores} -o ./../log/job.out.pipeline" \
    --jobs 100 \
    --max-jobs-per-second 1 \
    -T 0 \
    -w 15 \
    --default-resources mem_gb=4 time_min=60 num_cores=1 

