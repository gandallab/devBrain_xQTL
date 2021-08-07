#!/bin/bash
#$ -l h_data=2G,h_rt=6:00:00
#$ -pe shared 6
#$ -cwd
#$ -j y
#$ -m a
#$ -t 1-120
#$ -o /u/project/gandalm/cindywen/isoform_twas/verifyBamID/log/job.out

. /u/local/Modules/default/init/modules.sh
module load samtools/1.9

INPUT_FILE=input.txt
INLINE=`head -n ${SGE_TASK_ID} ${INPUT_FILE} | tail -n 1`
IFS=$'\t' params=($INLINE)
VCF=${params[0]}
BAM=${params[1]}
STUDY=${params[2]}
SUBJ=${params[3]}

samtools index ${BAM}

mkdir -p ${STUDY}

/u/project/gandalm/shared/apps/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID \
    --vcf ${VCF} \
    --bam ${BAM} \
    --out ${STUDY}/${SUBJ} \
    --verbose \
    --ignoreRG \
    --smID ${SUBJ}
