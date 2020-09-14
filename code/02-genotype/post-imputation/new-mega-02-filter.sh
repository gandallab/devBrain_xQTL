#!/bin/bash
#$ -l h_data=20G,h_rt=8:00:00,highp
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.mega.filter
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load plink/1.90b3.45
module load htslib/1.9

R2_3_input=./merge.reheader.rsID.vcf.gz
output_prefix=./filtered

# In this step, chr notation is changed from chr1 to 1
plink --vcf ${R2_3_input} --hwe 1e-6 --maf 0.01 --geno 0.05 --recode vcf-iid bgz --out ${output_prefix}
tabix -p vcf ${output_prefix}.vcf.gz
