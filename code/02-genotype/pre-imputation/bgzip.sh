#!/bin/bash
#$ -l h_data=4G,h_rt=4:00:00,highp
#$ -cwd
#$ -j y
#$ -o ../../log/job.out.bgzip
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load htslib/1.12

bgzip /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Werling2020_eQTL/geno/to_impute/werling_sampleFiltered.vcf
bgzip /u/project/gandalm/shared/GenomicDatasets/FetalBrain/HDBR/geno/to_impute/HDBR_sampleFiltered.vcf
bgzip /u/project/gandalm/shared/GenomicDatasets/FetalBrain/Obrien2018_eQTL/geno/obrien2018_hg38.genotyped.vcf
