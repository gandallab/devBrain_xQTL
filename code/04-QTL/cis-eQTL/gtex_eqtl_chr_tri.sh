#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=4:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/eqtl_new/log/job.out.gtex.eqtl.chr.tri
#$ -m a
#$ -t 1-22

tri=$1
num_hcp=$2
chunk=${SGE_TASK_ID}

expr=/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/gene.counts.scaled.normalized.trimester${tri}.bed.gz
geno=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/eur/filtered.hg19.sorted.removeGeneOutlier.trimester${tri}.vcf.gz
cov=/u/project/gandalm/cindywen/isoform_twas/eqtl_new/data/eur/chuanjiao/tri${tri}_${num_hcp}HCP.txt
rel_file=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/related.txt
outdir=/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/

. /u/local/Modules/default/init/modules.sh
module load gcc/11.3.0
/u/project/gandalm/cindywen/isoform_twas/sqtl_new/fastqtl-gtex/bin/fastQTL.static \
            --vcf ${geno} \
            --bed ${expr} \
            --window 1e6 \
            --cov ${cov} \
            --seed 123 \
            --exclude-samples ${rel_file} \
            --chunk ${chunk} 22 \
            --out ${outdir}tri${tri}_${num_hcp}_chunk${chunk}.txt.gz \
            --log ${outdir}tri${tri}_${num_hcp}_chunk${chunk}.log \
            --normal

