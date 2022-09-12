#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=1:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.ld.qtl.s
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624

locus=$1
feature=$2
trait=$3
level=$4

cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_${level}/${trait}/locus_${locus}

plink --vcf /u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeRel.vcf.gz \
      --r \
      --matrix \
      --extract ${feature}_gwas_qtl_snp_set.txt \
      --out ${feature}_qtl
