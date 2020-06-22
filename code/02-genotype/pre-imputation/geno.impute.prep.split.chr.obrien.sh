#!/bin/bash
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Obrien2018_eQTL/geno/to_impute/
#$ -j y
#$ -o ./job.out.split.chr
#$ -m a
#$ -t 1-22

. /u/local/Modules/default/init/modules.sh
module load bcftools/1.9
module load htslib/1.9
module load vcftools/0.1.14

#from filtered vcf file, split by chr, and sort
chr=${SGE_TASK_ID}
bcftools view -r chr${chr} obrien_filtered_withChr.vcf.gz | vcf-sort | bgzip -c > chr${chr}.vcf.gz
