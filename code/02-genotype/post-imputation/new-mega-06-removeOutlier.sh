#!/bin/bash
#$ -l h_data=4G,h_rt=4:00:00,highp
#$ -wd /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/
#$ -o /u/home/c/cindywen/project-gandalm/isoform_twas/genotype/log/job.out.mega.removeOutlier
#$ -j y
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib


# Gene
#outlier_list=~/project-gandalm/isoform_twas/eqtl/data/gene.outlier.txt
# Tx
outlier_list=~/project-gandalm/isoform_twas/isoqtl/data/tx.outlier.txt


#input=data.hg19.sorted.vcf.gz
#output=data.hg19.removeGeneOutlier.vcf.gz
#output=data.hg19.removeTxOutlier.vcf.gz

#input=./eur/eur.hg19.sorted.vcf.gz
#output=./eur/eur.hg19.removeGeneOutlier.vcf.gz
#output=./eur/eur.hg19.removeTxOutlier.vcf.gz

#input=./amr/amr.hg19.sorted.vcf.gz
#output=./amr/amr.hg19.removeGeneOutlier.vcf.gz
#output=./amr/amr.hg19.removeTxOutlier.vcf.gz

input=./afr/afr.hg19.sorted.vcf.gz
#output=./afr/afr.hg19.removeGeneOutlier.vcf.gz
output=./afr/afr.hg19.removeTxOutlier.vcf.gz

bcftools view -S ^${outlier_list} --force-samples -Oz -o ${output} ${input}
tabix -p vcf ${output}
