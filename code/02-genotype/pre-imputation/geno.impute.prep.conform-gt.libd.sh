#!/bin/bash
#$ -l h_data=20G,h_rt=4:00:00,highp
#$ -wd /u/project/gandalm/shared/GenomicDatasets/FetalBrain/LIBD_1and2/geno/
#$ -j y
#$ -o ./job.out.conform-gt
#$ -m a
#$ -t 1-22

. /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load vcftools/0.1.14
module load htslib/1.9

chr=${SGE_TASK_ID}
ref=/u/project/gandalm/shared/refGenomes/1000genomes/chrs/chr${chr}.1kg.phase3.v5a.vcf.gz
# for chr2, use fixed reference file
#ref=/u/home/c/cindywen/project-gandalm/isoform_twas/LIBD_data/chr2.1kg.phase3.v5a.final.vcf.gz

phase1and2_1M=./phase1and2/to_impute/libd.phase1and2.1M.filtered.vcf.gz
phase1and2_h650=./phase1and2/to_impute/libd.phase1and2.h650.filtered.vcf.gz
phase2only_1M=./phase2only/to_impute/libd.phase2only.1M.filtered.vcf.gz
phase2only_h650=./phase2only/to_impute/libd.phase2only.h650.filtered.vcf.gz
phase1and2_prefix_1M=libd.phase1and2.1M.conform-gt.chr${chr}
phase1and2_prefix_h650=libd.phase1and2.h650.conform-gt.chr${chr}
phase2only_prefix_1M=libd.phase2only.1M.conform-gt.chr${chr}
phase2only_prefix_h650=libd.phase2only.h650.conform-gt.chr${chr}

java -jar /u/home/c/cindywen/project-gandalm/isoform_twas/LIBD_data/conform-gt.24May16.cee.jar ref=${ref} gt=${phase1and2_1M} chrom=${chr} out=phase1and2/to_impute/conform-gt/${phase1and2_prefix_1M}
java -jar /u/home/c/cindywen/project-gandalm/isoform_twas/LIBD_data/conform-gt.24May16.cee.jar ref=${ref} gt=${phase1and2_h650} chrom=${chr} out=phase1and2/to_impute/conform-gt/${phase1and2_prefix_h650}
java -jar /u/home/c/cindywen/project-gandalm/isoform_twas/LIBD_data/conform-gt.24May16.cee.jar ref=${ref} gt=${phase2only_1M} chrom=${chr} out=phase2only/to_impute/conform-gt/${phase2only_prefix_1M}
java -jar /u/home/c/cindywen/project-gandalm/isoform_twas/LIBD_data/conform-gt.24May16.cee.jar ref=${ref} gt=${phase2only_h650} chrom=${chr} out=phase2only/to_impute/conform-gt/${phase2only_prefix_h650}
# default: match=ID, both are rsID

vcf-sort phase1and2/to_impute/conform-gt/${phase1and2_prefix_1M}.vcf.gz | bgzip -c > phase1and2/to_impute/conform-gt/${phase1and2_prefix_1M}.sorted.vcf.gz
vcf-sort phase1and2/to_impute/conform-gt/${phase1and2_prefix_h650}.vcf.gz | bgzip -c > phase1and2/to_impute/conform-gt/${phase1and2_prefix_h650}.sorted.vcf.gz
vcf-sort phase2only/to_impute/conform-gt/${phase2only_prefix_1M}.vcf.gz | bgzip -c > phase2only/to_impute/conform-gt/${phase2only_prefix_1M}.sorted.vcf.gz
vcf-sort phase2only/to_impute/conform-gt/${phase2only_prefix_h650}.vcf.gz | bgzip -c > phase2only/to_impute/conform-gt/${phase2only_prefix_h650}.sorted.vcf.gz
