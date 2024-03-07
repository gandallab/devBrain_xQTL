#!/bin/bash -l
#$ -wd /u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/
#$ -l h_data=8G,h_rt=24:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.pair1
#$ -m a

# join -t $' ' -1 1 -2 1 -a 1 -o '1.1 2.2' <(sort -t $'\t' -k1,1 Walker_perm_10hcp/sig_pheno_pairs.txt) <(zcat HDBR_nominal_20hcp/all.chunks.prep.txt.gz | sort -t $'\t' -k1,1) > HDBR_nominal_20hcp/walker_permsig_in_hdbr1.txt

join -t $' ' -1 1 -2 1 -a 1 -o '1.1 2.2' <(sort -t $'\t' -k1,1 mixed_perm_40hcp_1e6/sig_pheno_pairs.txt) <(zcat HDBR_nominal_20hcp/all.chunks.prep.txt.gz | sort -t $'\t' -k1,1) > HDBR_nominal_20hcp/mega_permsig_in_hdbr1.txt