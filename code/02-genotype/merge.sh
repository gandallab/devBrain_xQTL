#!/bin/bash
#$ -l h_data=4G,h_rt=96:00:00,highp
#$ -wd /u/project/gandalm/cindywen/isoform_twas/genotype/
#$ -V
#$ -j y
#$ -o ./job_out_merge
#$ -m a
#$ -t 1-22

. /u/local/Modules/default/init/modules.sh
module load vcftools/0.1.14
module load htslib

chr=${SGE_TASK_ID}
Walker=./Walker/imputation_results/chr${chr}.dose.vcf.gz
Obrien=./Obrien/imputation_results/chr${chr}.dose.vcf.gz
if [ ${chr} -ge 1 ] && [ ${chr} -le 6 ]; then 
	Werling=./Werling/imputation_results/chr${chr}/chr${chr}.dose.vcf.gz
elif [ ${chr} -ge 7 ] && [ ${chr} -le 9 ]; then
	Werling=./Werling/imputation_results/chr789/chr${chr}.dose.vcf.gz
elif [ ${chr} -ge 10 ] && [ ${chr} -le 12]; then
	Werling=./Werling/imputation_results/chr101112/chr${chr}.dose.vcf.gz
elif [ ${chr} -ge 13 ] && [ ${chr} -le 15 ]; then
	Werling=./Werling/imputation_results/chr131415/chr${chr}.dose.vcf.gz
elif [ ${chr} -ge 16 ] && [ ${chr} -le 18 ]; then
	Werling=./Werling/imputation_results/chr161718/chr${chr}.dose.vcf.gz
elif [ ${chr} -ge 19 ] && [ ${chr} -le 22 ]; then
	Werling=./Werling/imputation_results/chr19202122/chr${chr}.dose.vcf.gz
fi

tabix -p vcf ${Walker}
tabix -p vcf ${Obrien}
tabix -p vcf ${Werling}

out_folder=./merge/
mkdir -p ${out_folder}
vcf-merge ${Walker} ${Obrien} ${Werling} | bgzip -c > ${out_folder}/merge.chr${chr}.vcf.gz
