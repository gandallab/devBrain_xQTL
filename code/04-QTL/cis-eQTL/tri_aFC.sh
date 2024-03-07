#!/bin/bash -l
#$ -cwd
#$ -l h_data=6G,h_rt=4:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/eqtl_new/log/job.out.tri.aFC
#$ -m a
#$ -t 1-22

chr=${SGE_TASK_ID}
tri=$1

# quick and dirty: not excluding relatives?
expr=/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/gene.counts.scaled.normalized.trimester${tri}.bed.gz
geno=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeGeneOutlier.trimester${tri}.vcf.gz
cov=/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/chuanjiao/tri${tri}*txt
qtl=~/project-gandalm/isoform_twas/eqtl_new/results/eur_trimester/sig_pheno_trimester${tri}_aFC.txt
output=~/project-gandalm/isoform_twas/eqtl_new/results/eur_trimester/aFC/tri${tri}_${chr}.txt

. /u/local/Modules/default/init/modules.sh
module load python/3.9.6
python3 /u/project/gandalm/cindywen/isoform_twas/eqtl_new/code/scripts/aFC_DanielVo.py \
    --vcf ${geno} \
    --pheno ${expr} \
    --cov ${cov} \
    --output ${output} \
    --qtl ${qtl} \
    --log_xform 1 \
    --log_base 2 \
    --boot 100 \
    --geno GT \
    --chr ${chr}
