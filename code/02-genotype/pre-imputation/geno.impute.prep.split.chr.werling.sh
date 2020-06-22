#!/bin/bash
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Werling2020_eQTL/geno/to_impute/
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
bcftools view -r ${chr} werling_filtered.vcf.gz | vcf-sort | bgzip -c > chr${chr}.vcf.gz
# array build is hg38, so the following to change chr:
#for i in {1..22}; do zcat chr$i.vcf.gz | awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' > chr$i.withChr.vcf; done
#for i in {1..22}; do bgzip chr$i.withChr.vcf;done
