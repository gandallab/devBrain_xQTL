#!/bin/bash
#$ -l h_data=12G,h_rt=12:00:00,highp
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry
#$ -j y
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.mega.checkVCF
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load python/2.7.16

ref=/u/project/gandalm/shared/refGenomes/hg19/Gencode/GRCh37.primary_assembly.genome_noCHR.fa
input=./data.hg19.sorted.vcf.gz
mkdir checkVCF

python2 /u/home/c/cindywen/project-gandalm/isoform_twas/LIBD_data/scripts/checkVCF.py -r ${ref} -o ./checkVCF/data.hg19 ${input}
