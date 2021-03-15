#!/bin/bash
#$ -l h_data=4G,h_rt=24:00:00,highp
#$ -pe shared 16
#$ -cwd
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load samtools

STAR_INDEX=/u/project/gandalm/cindywen/isoform_twas/star/gencodev33_STARindex/
INPUT_FILE=$1
data=$2
GTF_FILE=/u/project/gandalm/shared/refGenomes/hg19/Gencode/v33/gencode.v33lift37.annotation.gtf
STAR=/u/project/gandalm/cindywen/isoform_twas/star/STAR-2.7.3a/bin/Linux_x86_64_static/STAR

INLINE=`head -n ${SGE_TASK_ID} ${INPUT_FILE} | tail -n 1`
IFS=$'\t' params=($INLINE)
SAMPLE_ID=${params[0]}
IN_FASTQS_R1=${params[1]}
IN_FASTQS_R2=${params[2]}

SAMPLEPATH=/u/project/gandalm/cindywen/isoform_twas/star/${data}_sqtl_new/${SAMPLE_ID}
mkdir -p ${SAMPLEPATH}

${STAR} --genomeDir $STAR_INDEX \
    --runThreadN 12 \
    --readFilesIn $IN_FASTQS_R1 $IN_FASTQS_R2 \
    --readFilesCommand gunzip -c \
    --outSAMattributes NH HI AS NM MD nM \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --sjdbGTFfile $GTF_FILE \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --outFilterMismatchNmax 999 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outFileNamePrefix $SAMPLEPATH/$SAMPLE_ID.STAR \
    --limitSjdbInsertNsj 1000000
