#!/bin/bash -l
#$ -cwd
#$ -l h_data=16G,h_rt=4:00:00
#$ -j y
#$ -o ../log/job.out.ecaviar.bash.new
#$ -m a
#$ -t 1-2256

. /u/local/Modules/default/init/modules.sh
module load gcc/10.2.0

list=gwas_feature_fetal_scz_eqtl_torun.txt
INLINE=`head -n ${SGE_TASK_ID} ${list} | tail -n 1`
IFS=$'\t' 
params=($INLINE)


locus=${params[0]}
gene=${params[1]}
gwas=${params[2]}
annot=${params[3]}


cd ../out_1MB/${gwas}/locus${locus}/

# note for fetal e/iso/sQTL, separate LD files; for fetal trimester QTL and adult (thistle, MB), same LD file
/u/project/gandalm/shared/apps/caviar/CAVIAR-C++/eCAVIAR \
    -o ${gene}_fetal_${annot}_ecaviar \
    -l ${gene}_fetal_${annot}_gwas.ld \
    -z ${gene}_fetal_${annot}_gwas_zscore.txt \
    -l ${gene}_fetal_${annot}.ld \
    -z ${gene}_fetal_${annot}_zscore.txt \
    -f 1 \
    -c 2 \
    -r 0.95
