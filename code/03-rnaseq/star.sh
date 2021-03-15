#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load samtools

STAR_INDEX=/u/project/gandalm/pampas/refGenomes/GRCh37/STARindex/
INPUT_FILE=$1
GTF_PATH=/u/project/gandalm/pampas/refGenomes/GRCh37/gencode.v29lift37.annotation.gtf

INLINE=`head -n ${SGE_TASK_ID} ${INPUT_FILE} | tail -n 1`
IFS=$'\t' params=($INLINE)
SAMPLE_ID=${params[0]}             
IN_FASTQS_R1=${params[1]}          
IN_FASTQS_R2=${params[2]}          

SAMPLEPATH=./hdbr/${SAMPLE_ID}
mkdir -p ${SAMPLEPATH}
# echo $SAMPLE_ID
# echo $IN_FASTQS_R1
# echo $IN_FASTQS_R2

# --limitSjdbInsertNsj 100000000 --limitBAMsortRAM 40000000000
# --sjdbGTFfile $GTF_PATH: include annotations on the fly at the mapping step, without including them at the genome generation step; not necessary here, as index was built with annotation.gtf
# --quantMode TranscriptomeSAM: default -; if given, in addition to genome alignments (Aligned.*.sam/bam), will output transcriptome coordinates alignments (Aligned.toTranscriptome.out.bam); can be used in quantifications like RSEM, we don't really need it here
./STAR-2.7.3a/bin/Linux_x86_64_static/STAR \
                    --genomeDir $STAR_INDEX \
                    --runThreadN 6 \
                    --readFilesIn $IN_FASTQS_R1 $IN_FASTQS_R2 \
                    --twopassMode Basic \
                    --sjdbGTFfile $GTF_PATH \
                    --readFilesCommand gunzip -c \
                    --outFilterMultimapNmax 20 \
                    --outSAMattributes NH HI AS NM MD nM \
                    --outSAMtype BAM SortedByCoordinate \
                    --alignSJoverhangMin 8 \
                    --alignSJDBoverhangMin 1 \
                    --outFilterMismatchNmax 999 \
                    --alignIntronMin 20 \
                    --alignIntronMax 1000000 \
                    --alignMatesGapMax 1000000 \
                    --outSAMunmapped Within \
                    --outFilterMismatchNoverReadLmax 0.04 \
                    --sjdbScore 1 \
                    --quantMode TranscriptomeSAM \
                    --limitSjdbInsertNsj 100000000 \
                    --outFileNamePrefix $SAMPLEPATH/$SAMPLE_ID.STAR

# Sort transcriptome coordinates alignments
mv $SAMPLEPATH/$SAMPLE_ID.STARAligned.toTranscriptome.out.bam $SAMPLEPATH/Tr.bam
cat <( samtools view -H $SAMPLEPATH/Tr.bam ) <( samtools view -@ 6 $SAMPLEPATH/Tr.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S 40G -T ./ | tr ' ' '\n' ) | samtools view -@ 6 -bS - > $SAMPLEPATH/$SAMPLE_ID.STARAligned.toTranscriptome.out_sorted.bam
rm $SAMPLEPATH/Tr.bam
touch $SAMPLEPATH/$SAMPLE_ID.done
