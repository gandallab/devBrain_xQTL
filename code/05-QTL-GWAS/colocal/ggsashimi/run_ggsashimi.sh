#!/bin/bash -l
#$ -cwd
#$ -l h_data=16G,h_rt=96:00:00,highp
#$ -j y
#$ -o job.out
#$ -m a

. /u/local/Modules/default/init/modules.sh
module load R/4.1.0-BIO
module load python/3.9.6

/u/project/gandalm/shared/apps/ggsashimi/ggsashimi.py \
    -b rs10276352_input_bams.tsv \
    -c chr7:21467661-21554440 \
    -g /u/project/gandalm/shared/refGenomes/hg19/Gencode/v33/gencode.v33lift37.annotation.gtf \
    -C 3 \
    -O 3 \
    --shrink \
    --alpha 0.25 \
    --base-size=20 \
    --ann-height=4 \
    --height=3 \
    --width=18 \
    -P /u/project/gandalm/shared/apps/ggsashimi/examples/palette.txt \
    -o sashimi3


    # -M 10 \
