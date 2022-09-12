#!/bin/bash -l
#$ -cwd
#$ -l h_data=4G,h_rt=2:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/colocal/log/job.out.analyze.iso.s
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

locus=$SGE_TASK_ID
level=$1 #isoqtl, sqtl, male_eqtl, female_eqtl, tri1_eqtl, tri2_eqtl
trait=$2
feature=$3 #isoform, intron, egene
cd /u/project/gandalm/cindywen/isoform_twas/colocal/results_${level}/${trait}/locus_${locus}/

cat locus_${feature}.txt | while read gene
do Rscript /u/project/gandalm/cindywen/isoform_twas/colocal/code/scripts/analyze_iso_s.R \
    --locus $locus \
    --gene $gene \
    --trait $trait \
    --level $level
done

# touch analye.done
