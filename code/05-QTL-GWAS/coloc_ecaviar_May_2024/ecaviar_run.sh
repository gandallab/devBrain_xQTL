#!/bin/bash -l 
#$ -cwd
#$ -l h_data=30G,h_rt=24:00:00,highp
#$ -j y
#$ -o ../log/job.out.snakemake.ecaviar.run
#$ -m a

# source /u/local/apps/anaconda3/2019.03/bin/activate snakemake
module load anaconda3
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate snakemake

snakemake \
    --snakefile ecaviar_run.smk \
    --jobscript cluster_jobscript.sh \
    --cluster-sync "qsub -l h_data={resources.mem_gb}G,h_rt=00:{resources.time_min}:00,highp -pe shared {resources.num_cores} -o ../log/job.out.ecaviar.run" \
    --jobs 200 \
    --max-jobs-per-second 10 \
    -T 0 \
    -w 60 \
    --default-resources mem_gb=20 time_min=1440 num_cores=4 \
    --rerun-incomplete 

