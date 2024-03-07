#!/bin/bash
#$ -l h_data=4G,h_rt=4:00:00,highp
#$ -cwd
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.gtex.fastqtl
#$ -m a
#$ -t 1-100

expr=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/cluster/leafcutter_perind.counts.nochr.gz.qqnorm_all_fixSubj_combat.bed.gz
geno=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/filtered.hg19.sorted.removeRel.vcf.gz
cov=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/data/40hcp_cov.txt
rel_file=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt
out=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/gtex_chunk${SGE_TASK_ID}

. /u/local/Modules/default/init/modules.sh
module load gcc/10.2.0
/u/project/gandalm/cindywen/isoform_twas/sqtl_new/fastqtl-gtex/bin/fastQTL.static \
    --vcf ${geno} \
    --bed ${expr} \
    --window 1e6 \
    --cov ${cov} \
    --seed 123 \
    --exclude-samples ${rel_file} \
    --chunk ${SGE_TASK_ID} 100 \
    --out ${out}.txt.gz \
    --log ${out}.log \
    --normal
