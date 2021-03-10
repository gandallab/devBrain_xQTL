#!/bin/bash

INPUT_FILE=$1
INLINE=`head -n ${SGE_TASK_ID} ${INPUT_FILE} | tail -n 1`
IFS=$'\t' params=($INLINE)
SAMPLE_ID=${params[0]}             
IN_FASTQS_R1=${params[1]}          
IN_FASTQS_R2=${params[2]}          
SAMPLEPATH=./hdbr/${SAMPLE_ID}
mkdir -p ${SAMPLEPATH}
./FastQC/fastqc --outdir ${SAMPLEPATH} ${IN_FASTQS_R1} ${IN_FASTQS_R2}
