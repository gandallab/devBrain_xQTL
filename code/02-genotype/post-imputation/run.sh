#!/bin/bash -l
#$ -cwd
#$ -l h_data=20G,h_rt=4:00:00
#$ -j y
#$ -o job.out.afr
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load python/3.7.0
module load htslib/1.9
module load bcftools/1.9
module load plink

for i in {1..22}; do plink --bfile ../all_data/isec_R2_greater_than_3/ancestry/merge.reheader.chr${i}_rsid --keep ../all_data/isec_R2_greater_than_3/ancestry/ancestry_list/plink_afr.tsv --make-bed --out afr.pre.filter.${i};done

for i in {1..22}; do echo "afr.pre.filter.${i}" >> afr.mergelist.txt; done

plink --merge-list afr.mergelist.txt --make-bed --out afr.pre.filter

plink --bfile afr.pre.filter --recode vcf-iid bgz --out afr.pre.filter

tabix -p vcf afr.pre.filter.vcf.gz

python3 ~/.local/bin/CrossMap.py vcf /u/project/gandalm/cindywen/isoform_twas/genotype/ref/GRCh38_to_GRCh37.chain.gz afr.pre.filter.vcf.gz /u/project/gandalm/shared/refGenomes/hg19/Gencode/GRCh37.primary_assembly.genome.fa afr.pre.filter.hg19.vcf
# EUR log
# @ 2021-08-26 15:40:56: Read the chain file:  /u/project/gandalm/cindywen/isoform_twas/genotype/ref/GRCh38_to_GRCh37.chain.gz
# @ 2021-08-26 15:40:58: Filter out variants [reference_allele == alternative_allele] ...
# @ 2021-08-26 15:40:58: Updating contig field ...
# @ 2021-08-26 15:40:58: Lifting over ...

# @ 2021-08-26 15:55:41: Total entries: 10959319
# @ 2021-08-26 15:55:41: Failed to map: 1649821
# AFR log
# @ 2021-08-26 20:00:20: Total entries: 10959319
# @ 2021-08-26 20:00:20: Failed to map: 1665805
bgzip afr.pre.filter.hg19.vcf

bcftools sort afr.pre.filter.hg19.vcf.gz -Oz -o afr.pre.filter.hg19.sorted.vcf.gz
tabix -p vcf afr.pre.filter.hg19.sorted.vcf.gz

bcftools +fill-tags afr.pre.filter.hg19.sorted.vcf.gz -- -t MAF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%MAF\n' | bgzip -c > afr.maf.gz
