#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load R/3.6.0

logfile=$1
chr=$2
cutoff=$(head -n1 ./log_residQN_h2g | cut -d" " -f2)
# for i in {1..22}; do Rscript scripts/summarize_sig_genes_fromP.R $i $cutoff; done
# for i in {1..22}; do Rscript scripts/summarize_gene_pair.R $i $cutoff; done
Rscript scripts/summarize_gene_pair.R $chr $cutoff