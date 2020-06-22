#!/bin/bash
#$ -l h_data=4G,h_rt=2:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Obrien2018_eQTL/geno/
#$ -j y
#$ -o ./to_impute/job.out.impute.prep
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load plink/1.90b3.45
module load htslib/1.9

# starting vcf file is imputed and filtered
# use GENOTYPED SNPs only
#bcftools view -e "FILTER=='PASS'" $file -o $out

# apply plink filters
file=./obrien2018_hg38.genotyped.vcf
plink --vcf ${file} --maf 0.01 --geno 0.05 --mind 0.1 --hwe 1e-6 --recode vcf-iid --out ./to_impute/obrien_filtered
# change 1 to chr1
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ./to_impute/obrien_filtered.vcf > ./to_impute/obrien_filtered_withChr.vcf

bgzip ./to_impute/obrien_filtered_withChr.vcf
tabix -p vcf ./to_impute/obrien_filtered_withChr.vcf.gz

