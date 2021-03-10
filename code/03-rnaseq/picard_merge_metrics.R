setwd("~/Desktop/fetalBrainQCreport/picard/all_multiqc/multiqc_data/")
alignment <- read.delim("multiqc_picard_AlignmentSummaryMetrics.txt", header = TRUE)
dups <- read.delim("multiqc_picard_dups.txt", header = TRUE)
gcbias <- read.delim("multiqc_picard_gcbias.txt", header = TRUE)
insert <- read.delim("multiqc_picard_insertSize.txt", header = TRUE)
rnaseq <- read.delim("multiqc_picard_RnaSeqMetrics.txt", header = TRUE)

# drop NA columns
drops <- c("SAMPLE","LIBRARY", "READ_GROUP")
alignment <- alignment[ , !(names(alignment) %in% drops)]
gcbias <- gcbias[ , !(names(gcbias) %in% drops)]
insert <- insert[ , !(names(insert) %in% drops)]
rnaseq <- rnaseq[ , !(names(rnaseq) %in% drops)]

# merge metrics
temp <- merge(alignment, dups, by="Sample", all.y = TRUE)
temp_2 <- merge(temp, gcbias, by="Sample", all.y = TRUE)
temp_3 <- merge(temp_2, insert, by.x = "Sample", by.y = "SAMPLE_NAME", all.y = TRUE)
picard_merged <- merge(temp_3, rnaseq, by="Sample", all.y = TRUE)

picard <- picard_merged[,-1]
rownames(picard) <- picard_merged[,1]
# drop non-numeric/no variance columns
drop <- c('CATEGORY','LIBRARY','ACCUMULATION_LEVEL','READS_USED','Sample.y','PAIR_ORIENTATION')
picard <- picard[ , !(names(picard) %in% drop)]
# remove columns with zero variance
picard <- picard[ , which(apply(picard, 2, var) != 0)]
write.table(picard, file = "~/Desktop/picard_QC_compiled.tsv", row.names = TRUE, col.names = TRUE)

# calculate seqPC, not used in this project
# seq_qc_pca <- prcomp(picard, scale. = TRUE)
# library(factoextra)
# # show percentage of variance explained by each PC
# fviz_eig(seq_qc_pca)
