#!/bin/bash
#$ -l h_data=32G,h_rt=02:00:00,highp
#$ -wd /u/project/gandalm/
#$ -j y
#$ -m a
#$ -t 1-4

list=./cindywen/isoform_twas/FINAL/libd.phase2only.fastq.txt
INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
IFS=$'\t' params=($INLINE)
SUBJECT_ID=${params[0]}             #BrXXXX
FASTQ_PREFIX=${params[1]}           #RXXXX
FASTQ_PATH=./shared/GenomicDatasets/FetalBrain/LIBD_phase2/rnaseq/fastq
CONCAT_PATH=./shared/GenomicDatasets/FetalBrain/LIBD_1and2/rnaseq

cat ${FASTQ_PATH}/${FASTQ_PREFIX}_*R1_001.fastq.gz >> ${CONCAT_PATH}/${SUBJECT_ID}_R1_001.fastq.gz
cat ${FASTQ_PATH}/${FASTQ_PREFIX}_*R2_001.fastq.gz >> ${CONCAT_PATH}/${SUBJECT_ID}_R2_001.fastq.gz

touch ${CONCAT_PATH}/${SUBJECT_ID}.done
