#!/bin/bash
#$ -l h_data=4G,h_rt=1:00:00
#$ -cwd
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.sqtl.extract.test
#$ -m a
#$ -j y

. /u/local/etc/profile.d/sge.sh

trait=$1
locus=$2

assoc=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/gtex.allpairs.txt.gz
list=/u/project/gandalm/cindywen/isoform_twas/colocal/results_sqtl/${trait}/locus_${locus}/locus_intron.txt

intron=`head -n ${SGE_TASK_ID} ${list} | tail -n1`

qsub /u/project/gandalm/cindywen/isoform_twas/colocal/code/scripts/extract_cis_assoc.sh ${intron} ${assoc} /u/project/gandalm/cindywen/isoform_twas/colocal/results_sqtl/${trait}/locus_${locus}/

