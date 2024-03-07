#!/bin/bash
#$ -l h_data=4G,h_rt=24:00:00,highp
#$ -pe shared 8
#$ -wd /u/project/gandalm/cindywen/isoform_twas/star/
#$ -j y
#$ -m a

# Custom filtering:
# (1) all MT entries,
# (2) all non-canonical SJ,
# (3) SJ supported by multi-mappers-only,
# (4) SJ supported by too few reads (<=2)

cat ./*_sqtl_new/*/*.STARSJ.out.tab \
| awk '($5 > 0 && $7 > 2 && $6==0 && $1!="chrM" ){out[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]++} END{for(var in out){print var}}' \
| sort -t $'\t' -k1,1 -k2,2n -k3,3n \
> ./gencode_v33lift37_new_sqtl_SJ.filtered.tab
