#!/bin/bash
#$ -l h_data=32G,h_rt=8:00:00,highp
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.post.impute.03.isec.merge
#$ -j y
#$ -m a

# merge datasets together
. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9

R2_3_0000=./isec_R2_greater_than_3/0000.vcf.gz
R2_3_0001=./isec_R2_greater_than_3/0001.vcf.gz
R2_3_0002=./isec_R2_greater_than_3/0002.vcf.gz
R2_3_0003=./isec_R2_greater_than_3/0003.vcf.gz
R2_3_0004=./isec_R2_greater_than_3/0004.vcf.gz

R2_8_0000=./isec_R2_greater_than_8/0000.vcf.gz
R2_8_0001=./isec_R2_greater_than_8/0001.vcf.gz
R2_8_0002=./isec_R2_greater_than_8/0002.vcf.gz
R2_8_0003=./isec_R2_greater_than_8/0003.vcf.gz
R2_8_0004=./isec_R2_greater_than_8/0004.vcf.gz
bcftools merge -m id ${R2_3_0000} ${R2_3_0001} ${R2_3_0002} ${R2_3_0003} ${R2_3_0004} -Oz -o isec_R2_greater_than_3/merge.vcf.gz
bcftools merge -m id ${R2_8_0000} ${R2_8_0001} ${R2_8_0002} ${R2_8_0003} ${R2_8_0004} -Oz -o isec_R2_greater_than_8/merge.vcf.gz
