#!/bin/bash
#$ -l h_data=48G,h_rt=8:00:00,highp
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/
#$ -j y
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.mega.03.ancestry
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9
module load python/3.7.0
module load plink/1.90b3.45

for an in eur amr afr
do
  # 1. crossmap to hg19
  input=./${an}/${an}.filtered.vcf.gz
  chain=/u/home/c/cindywen/project-gandalm/isoform_twas/genotype/ref/hg38ToHg19.over.chain.gz
  ref=/u/home/c/cindywen/project-gandalm/isoform_twas/genotype/ref/hg19.fa
  output=./${an}/${an}.hg19.vcf
  python3 ~/.local/bin/CrossMap.py  vcf ${chain} ${input} ${ref} ${output}
  
  bgzip ${output}
  
  # 2. sort
  bcftools sort ${output}.gz -Oz -o ./${an}/${an}.hg19.sorted.vcf.gz
  tabix -p vcf ./${an}/${an}.hg19.sorted.vcf.gz
  
  #3. vcf2plink
  plink --vcf ./${an}/${an}.hg19.sorted.vcf.gz --allow-extra-chr --make-bed --out ./${an}/${an}.hg19
done
