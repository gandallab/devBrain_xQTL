#!/bin/bash
#$ -l h_data=20G,h_rt=8:00:00,highp
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.ancestry.filter
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load plink/1.90b3.45
module load htslib/1.9

input=./merge.reheader.rsID.vcf.gz
mkdir eur
mkdir afr
mkdir amr

output_prefix_eur=./eur/eur.filtered
output_prefix_afr=./afr/afr.filtered
output_prefix_amr=./amr/amr.filtered

# In this step, chr notation is changed from chr1 to 1
plink --vcf ${input} --keep ~/project-gandalm/isoform_twas/eqtl/data/ancestry_list/plink.keep.eur.tsv --hwe 1e-6 --maf 0.01 --geno 0.05 --recode vcf-iid bgz --out ${output_prefix_eur}
tabix -p vcf ${output_prefix_eur}.vcf.gz
plink --vcf ${input} --keep ~/project-gandalm/isoform_twas/eqtl/data/ancestry_list/plink.keep.afr.tsv --hwe 1e-6 --maf 0.01 --geno 0.05 --recode vcf-iid bgz --out ${output_prefix_afr}
tabix -p vcf ${output_prefix_afr}.vcf.gz
plink --vcf ${input} --keep ~/project-gandalm/isoform_twas/eqtl/data/ancestry_list/plink.keep.amr.tsv --hwe 1e-6 --maf 0.01 --geno 0.05 --recode vcf-iid bgz --out ${output_prefix_amr}
tabix -p vcf ${output_prefix_amr}.vcf.gz
