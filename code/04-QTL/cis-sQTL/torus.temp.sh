#!/bin/bash -l
#$ -wd /u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/
#$ -l h_data=4G,h_rt=30:00:00
#$ -pe shared 12
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/sqtl_new/log/job.out.torus.CW
#$ -m a

annot=/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/annot/mixed_variant_annot_ensembl_vep.txt.gz
fastqtl=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/gtex.allpairs.torus.txt.gz
output=/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/mixed_nominal_40hcp_1e6/torus.out

. /u/local/Modules/default/init/modules.sh
module load python/3.7.2
/u/project/gandalm/shared/apps/torus-master/src/torus \
    -d ${fastqtl} \
    --fastqtl \
    -est \
    -annot ${annot} \
    > ${output}
