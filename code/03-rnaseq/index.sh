#!/bin/bash
#$ -l h_data=4G,h_rt=4:00:00,highp
#$ -cwd
#$ -pe shared 4
#$ -j y
#$ -o /u/project/gandalm/cindywen/isoform_twas/salmon/log/job.out.index
#$ -m a

source ~/anaconda3/bin/activate SalmonTools
# Salmon index used has 1 few tx mapped to v29 gene, not v33
# Here, re-run generateDecoyTranscriptome script to generate gentrome.fa to check
# The script uses command realpath which is not on H2
# function from Connor
function realpath()
{
    f=$@
    if [ -d "$f" ]; then
        base=""
        dir="$f"
    else
        base="/$(basename "$f")"
        dir=$(dirname "$f")
    fi
    dir=$(cd "$dir" && /bin/pwd)
    echo "$dir$base"
}

export -f realpath

/u/project/gandalm/shared/refGenomes/hg19/Gencode/v33/generateDecoyTranscriptome.sh \
    -j 4 \
    -a /u/project/gandalm/shared/refGenomes/hg19/Gencode/v33/gencode.v33lift37.annotation.gtf \
    -g /u/project/gandalm/shared/refGenomes/hg19/Gencode/GRCh37.primary_assembly.genome.fa \
    -t /u/project/gandalm/shared/refGenomes/hg19/Gencode/v33/gencode.v33lift37.transcripts.fa \
    -o /u/project/gandalm/cindywen/isoform_twas/salmon/salmon_index_gencodev33lift37_decoys
