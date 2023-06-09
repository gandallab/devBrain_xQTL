#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 4
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.walker.nominal
#$ -m a
#$ -t 1-100

num_hcp=$1
chunk=${SGE_TASK_ID}

expr=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/cluster/leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.bed.gz
geno=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeRel.vcf.gz
cov=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/${num_hcp}hcp_cov.txt
outdir=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/Walker_nominal_${num_hcp}hcp/
sample_list=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/code/walker_samples.txt


mkdir -p ${outdir}
/u/project/gandalm/shared/apps/FastQTL/bin/fastQTL.static \
    --vcf ${geno} \
    --bed ${expr} \
    --cov ${cov} \
    --out ${outdir}chunk${chunk}.txt.gz \
    --chunk ${chunk} 100 \
    --seed 123 \
    --window 1e6 \
    --normal \
    --include-samples ${sample_list}
