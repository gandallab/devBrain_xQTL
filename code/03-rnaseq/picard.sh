#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load java
module load R/3.5.1

#JAVA_BIN=/u/project/gandalm/pampas/bin/anaconda3/envs/rnaseq/bin/java
JAVA_FLAGS="-Xmx40g" #30G memory in JVM (Java Virtual Machine)
PICARD_TOOLS=./picard.jar

REF_FILE=/u/project/gandalm/pampas/refGenomes/GRCh37/GRCh37.primary_assembly.genome.fa
TX_FLAT_FILE=/u/project/gandalm/pampas/refGenomes/GRCh37/reorder.gencode.v29lift37.refFlat
RIBOSOMAL_INTERVALS=/u/project/gandalm/pampas/refGenomes/GRCh37/hg19.rRNA.interval_list
RNA_STRAND=SECOND_READ_TRANSCRIPTION_STRAND  # default: 1st read is reverse strand

PICARD_INPUT=./sortedBAM_list_hdbr.txt #STARAligned.sortedByCoord.out.bam
paramline=`head -n ${SGE_TASK_ID} ${PICARD_INPUT} | tail -n 1`
params=($paramline)
SAMPLE_ID=${params[0]}
IN_BAM=${params[1]}

mkdir -p ./hdbr/${SAMPLE_ID}
METRICS_OUT_BASE=./hdbr/${SAMPLE_ID}/${SAMPLE_ID}.picard
PICARD_DONE=./hdbr/status/${SAMPLE_ID}
mkdir -p ${PICARD_DONE}

if [ ! -f $PICARD_DONE ]; then
  ALIGNMENT_METRICS_OUT=$METRICS_OUT_BASE.alignment_metrics.txt
  RNA_SEQ_METRICS_OUT=$METRICS_OUT_BASE.rna_seq_metrics.txt
  GC_METRICS_OUT=$METRICS_OUT_BASE.gc_bias_metrics.txt
  GC_SUMMARY_OUT=$METRICS_OUT_BASE.gc_bias_metrics.summary.txt
  CHART_OUT=$METRICS_OUT_BASE.gc_bias.chart.pdf
  DUP_METRICS_OUT=$METRICS_OUT_BASE.duplication_metrics.txt
  DUP_BAM_OUT=$METRICS_OUT_BASE.temp.dup.bam
  INSERT_SIZE_OUT=$METRICS_OUT_BASE.insert_size_metrics.txt
  INSERT_HIST_OUT=$METRICS_OUT_BASE.insert_size_histogram.pdf
  MULTIPLE_METRICS_OUT=$METRICS_OUT_BASE.multiple_metrics.txt

  java $JAVA_FLAGS -jar $PICARD_TOOLS CollectAlignmentSummaryMetrics \
  R=$REF_FILE \
  O=$ALIGNMENT_METRICS_OUT \
  I=$IN_BAM
  SUCCESS_AM=$?

  java $JAVA_FLAGS -jar $PICARD_TOOLS CollectRnaSeqMetrics \
  I=$IN_BAM \
  O=$RNA_SEQ_METRICS_OUT \
  STRAND=$RNA_STRAND \
  RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS \
  REF_FLAT=$TX_FLAT_FILE
  SUCCESS_RM=$?

  java $JAVA_FLAGS -jar $PICARD_TOOLS CollectGcBiasMetrics \
  R=$REF_FILE \
  I=$IN_BAM \
  O=$GC_METRICS_OUT \
  S=$GC_SUMMARY_OUT \
  CHART=$CHART_OUT
  SUCCESS_GM=$?  # always fails unless Rscript is in the $PATH (usually no)

  java $JAVA_FLAGS -jar $PICARD_TOOLS MarkDuplicates \
  I=$IN_BAM \
  O=$DUP_BAM_OUT \
  M=$DUP_METRICS_OUT
  SUCCESS_MD=$?

  java $JAVA_FLAGS -jar $PICARD_TOOLS CollectInsertSizeMetrics \
  I=$IN_BAM \
  O=$INSERT_SIZE_OUT \
  H=$INSERT_HIST_OUT
  SUCCESS_IS=$?
  if [ $SUCCESS_MD -eq 0 ]; then
    rm $DUP_BAM_OUT*
  fi
  if [ $SUCCESS_AM -eq 0 ] && \
     [ $SUCCESS_RM -eq 0 ] && \
     [ -f $GC_SUMMARY_OUT ] && \
     [ $SUCCESS_IS -eq 0 ] && \
     [ $SUCCESS_MD -eq 0 ]; then
    touch $PICARD_DONE
  fi
fi
