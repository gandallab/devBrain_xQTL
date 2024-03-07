#!/bin/bash
#$ -l h_data=2G,h_rt=12:00:00
#$ -pe shared 6
#$ -cwd
#$ -j y
#$ -m a
#$ -t 1-654
#$ -o /u/project/gandalm/cindywen/isoform_twas/star/log/job.out.wasp.fix.truncated

. /u/local/Modules/default/init/modules.sh
module load samtools/1.9

INPUT_FILE=allBAMprefix.txt

INLINE=`head -n ${SGE_TASK_ID} ${INPUT_FILE} | tail -n 1`
IFS=$'\t' params=($INLINE)
IN_BAM=${params[0]}.STARAligned.sortedByCoord.out.bam
TEMP_SAM=${params[0]}.STARAligned.sortedByCoord.WASPfiltered.out.sam
OUT_BAM=${params[0]}.STARAligned.sortedByCoord.WASPfiltered.out.bam

# view original bam, keep header, WASP filter, output to sam
samtools view -h ${IN_BAM} | grep -v "vW:i:[2-7]" > ${TEMP_SAM}
# convert sam to bam
samtools view -b ${TEMP_SAM} > ${OUT_BAM}
rm ${TEMP_SAM}
