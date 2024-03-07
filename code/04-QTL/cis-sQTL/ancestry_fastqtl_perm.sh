#!/bin/bash -l
#$ -cwd
#$ -l h_data=6G,h_rt=8:00:00
#$ -pe shared 4
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.perm.pop
#$ -m a
#$ -t 1-1000

ancestry=$1
num_hcp=$2

expr=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/${ancestry}/lc_combat.bed.gz
geno=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/${ancestry}/filtered.hg19.sorted.removeRel.vcf.gz
cov=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/${ancestry}/${num_hcp}HCP_cov.txt
rel=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt

chunk=${SGE_TASK_ID}
outdir=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/${ancestry}_perm_${num_hcp}HCP_1e6/


mkdir -p ${outdir}
/u/project/gandalm/shared/apps/FastQTL/bin/fastQTL.static \
    --vcf ${geno} \
    --bed ${expr} \
    --cov ${cov} \
    --permute 1000 10000 \
    --out ${outdir}chunk${chunk}.txt.gz \
    --chunk ${chunk} 1000 \
    --seed 123 \
    --window 1e6 \
    --normal \
    --exclude-samples ${rel}
touch ${outdir}chunk${chunk}.txt.gz.done
