#!/bin/bash -l
#$ -wd /u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/
#$ -l h_data=4G,h_rt=4:00:00
#$ -pe shared 4
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.prep
#$ -m a

# zcat Walker_nominal_10hcp/all.chunks.txt.gz | awk '{printf "%s-%s %s\n", $1, $2, $4}' | gzip -c > Walker_nominal_10hcp/all.chunks.prep.txt.gz
# zcat HDBR_nominal_20hcp/all.chunks.txt.gz | awk '{printf "%s-%s %s\n", $1, $2, $4}' | gzip -c > HDBR_nominal_20hcp/all.chunks.prep.txt.gz

# duplicated lines exist in file. Here only keep the first occurance of duplicated lines. But do they have the same stats?
zcat mixed_nominal_40hcp_1e6/all.chunks.txt.gz | awk '{printf "%s-%s %s\n", $1, $2, $4}' | awk '!seen[$0]++' | gzip -c > mixed_nominal_40hcp_1e6/all.chunks.prep.txt.gz
