#!/bin/bash
#$ -l h_data=4G,h_rt=6:00:00,highp
#$ -cwd
#$ -j y
#$ -o ../log/job.out.bam2junc
#$ -m a
#$ -t 1-654

. /u/local/Modules/default/init/modules.sh
module load samtools/1.9
module load htslib/1.9

BAM_LIST=allBAMfiltered.txt
paramline=`head -n ${SGE_TASK_ID} ${BAM_LIST} | tail -n 1`
#params=($paramline)
BAM_FILE=${paramline[0]}

# 1707 in hdbr and walker
if ((${SGE_TASK_ID} == 155)); then 
	BAM_FILE_NAME=1707.1.STARAligned.sortedByCoord.WASPfiltered.out.bam
else
	BAM_FILE_NAME="$(basename ${BAM_FILE})"
fi

OUTDIR=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/junc
mkdir -p ${OUTDIR}

SCRIPT=/u/project/gandalm/cindywen/isoform_twas/sqtl/leafcutter/scripts/bam2junc.sh
bash ${SCRIPT} ${BAM_FILE} ${OUTDIR}/${BAM_FILE_NAME}.junc
echo ${OUTDIR}/${BAM_FILE_NAME}.junc >> juncfiles.txt
