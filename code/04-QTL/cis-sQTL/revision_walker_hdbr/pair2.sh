#!/bin/bash -l
#$ -wd /u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/
#$ -l h_data=8G,h_rt=24:00:00
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.pair2
#$ -m a

join -t $' ' -1 1 -2 1 -a 1 -o '1.1 2.2' <(sort -t $'\t' -k1,1 Walker_perm_10hcp/sig_pheno_pairs.txt) <(zcat mixed_nominal_40hcp_1e6/all.chunks.prep.txt.gz | sort -t $'\t' -k1,1) > mixed_nominal_40hcp_1e6/walker_permsig_in_mega.txt

join -t $' ' -1 1 -2 1 -a 1 -o '1.1 2.2' <(sort -t $'\t' -k1,1 HDBR_perm_20hcp/sig_pheno_pairs.txt) <(zcat mixed_nominal_40hcp_1e6/all.chunks.prep.txt.gz | sort -t $'\t' -k1,1) > mixed_nominal_40hcp_1e6/hdbr_permsig_in_mega.txt
