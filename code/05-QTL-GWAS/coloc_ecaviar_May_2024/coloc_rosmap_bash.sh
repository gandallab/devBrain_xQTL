#!/bin/bash -l
#$ -cwd
#$ -l h_data=16G,h_rt=4:00:00
#$ -j y
#$ -o ../log/job.out.coloc.rosmap.bash
#$ -m a
#$ -t 1-2149

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO

list=coloc_rosmap_temp.tsv
INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
IFS=$'\t'
params=($INLINE)


locus=${params[0]}
type=${params[1]}


Rscript scripts/coloc_rosmap_sneqtl.R \
            --gwas PGC3_SCZ_wave3.european.autosome.public.v3 \
            --locus ${locus} \
            --gwas_file /u/project/gandalm/shared/GWAS/SCZ.PGC3.2021/wave3_v3/PGC3_SCZ_wave3.european.autosome.public.v3.filtered.tsv.gz \
            --type ${type}
