#!/bin/bash
#$ -l h_data=13G,h_rt=4:00:00
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.after-isec-merge
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
#module load plink/1.90b3.45
module load htslib/1.9
module load bcftools/1.9

R2_3_input=./../merge.vcf.gz

# 1. remove 4 walker subjects
list=./walker.remove.subj.txt
bcftools view -S ^${list} ${R2_3_input} -Oz -o ./merge.654.vcf.gz

# 2. reheader to subject ID as in rnaseq
# old_name new_name, same order as in original vcf file
# or new_names for all samples, same order
# -o file.vcf will be bcf format
# can only output same format as input (BCF/VCF, bgzipped or not), so cannot pipe from -Ou
reheader=./geno.reheader2.txt
bcftools reheader -s ${reheader} ./merge.654.vcf.gz -o ./merge.reheader.vcf.gz
tabix -p vcf ./merge.reheader.vcf.gz

# 3. map to dbSNP rsID
# Note not all SNP id are converted to rsID, some are still chr:pos:REF:ALT
dbsnp=/u/project/gandalm/cindywen/isoform_twas/genotype/ref/00-common_all.withChr.vcf.gz
bcftools annotate -a ${dbsnp} -c ID ./merge.reheader.vcf.gz -Oz -o ./merge.reheader.rsID.vcf.gz

