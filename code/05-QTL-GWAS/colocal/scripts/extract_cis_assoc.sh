#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.extract.cis
#$ -m a

gene=$1
file=$2
outdir=$3

awk -v a="${gene}" '$1 == a {print}' <(zcat ${file}) > ${outdir}${gene}_all_pairs.txt
