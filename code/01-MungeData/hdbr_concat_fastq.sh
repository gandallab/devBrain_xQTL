#!/bin/bash
#$ -l h_data=32G,h_rt=02:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/HDBR/rnaseq/
#$ -j y
#$ -m a
#$ -t 1-173

# Merge fastq files, some samples have multiple regions sequenced

list=/u/home/c/cindywen/project-gandalm/isoform_twas/FINAL/hdbrOnlyFastqPrefix.txt
subjects=/u/home/c/cindywen/project-gandalm/isoform_twas/FINAL/hdbrOnlyIndividual.txt

INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
SUBJECT_ID=`head -n ${SGE_TASK_ID} ${subjects} | tail -n 1`
params=( ${INLINE} )
# count=${#params[@]}
for param in ${params[@]}; do cat ./${param}_1.fastq.gz >> ./concat/${SUBJECT_ID}_1.fastq.gz; done
for param in ${params[@]}; do cat ./${param}_2.fastq.gz >> ./concat/${SUBJECT_ID}_2.fastq.gz; done
touch ./concat/${SUBJECT_ID}.done
