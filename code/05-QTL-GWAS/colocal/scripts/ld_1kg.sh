#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=00:30:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.ld.1kg
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load plink/1.90b624

locus=$1
feature=$2
trait=$3
level=$4

chr=$(awk -v a="$locus" '$1 == a {print $2}' /u/project/gandalm/cindywen/isoform_twas/colocal/code/tables/${trait}_1Mb.txt)

cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_${level}/${trait}/locus_${locus}/

plink --bfile /u/project/gandalm/shared/LDSCORE/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
      --r \
      --matrix \
      --extract ${feature}_gwas_qtl_snp_set.txt \
      --out ${feature}_1kg_eur

