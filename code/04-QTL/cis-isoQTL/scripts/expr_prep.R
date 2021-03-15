#! /usr/bin/env Rscript
library(data.table)
library(argparser)
# library(edgeR)
library(DESeq2)
library(WGCNA)
library(sva)
library(dplyr)

# Input: scaled counts from salmon/tximport
# Prepare expression file for QTL mapping
# 0. match with GTF file, remove chrX, chrY, chrM
# 1. round expression for DESeq2
# 2. remove lowly expressed based on TPM
# 3. VST normalization and transformation
# 4. remove outlier subjects by connectivity
# 5. ComBat
# 6. pheno.bed for FastQTL; note strand; note 0/1-based

p <- arg_parser("QTL expression preparation")
p <- add_argument(p, "--gene_count", help="Salmon, tximport output")
p <- add_argument(p, "--gene_gtf", help="GENCODE annot")
p <- add_argument(p, "--gene_tpm", help="Salmon, tximport output")
p <- add_argument(p, "--geno_subj_dir", help="Subject lists")
p <- add_argument(p, "--outdir", help="")
args <- parse_args(p)


##### 1
datExpr <- fread(args$gene_count, data.table=F)
gtf <- fread(args$gene_gtf, data.table=F)
# remove chrM, chrX, chrY
# walker and hdbr duplicate 1707
# make hdbr 1707.1
names(datExpr)[647] <- "1707.1"
gtf <- gtf %>% filter(!(V1 %in% c("chrX","chrY","chrM")))
datExpr <- datExpr %>% filter(V1 %in% gtf$V5)
gtf <- gtf %>% filter(V5 %in% datExpr$V1)

rownames(datExpr) <- datExpr$V1
datExpr <- datExpr[,-1]
datExpr <- round(datExpr)

#####
# View the distribution of expression for each sample.
# box plot, looking for big differences in read depth (raw counts), symmetry in distribution across samples
# par(mfrow=c(1,2))
# boxplot(datExpr, range=0, main='Scaled Transcript Counts', xlab='Samples',xaxt="n")
# boxplot(log2(.1+datExpr), range=0, main = "log2(counts+.1)", xlab="Samples",xaxt="n")

# # Histogram/density plot
# # Look for: how well do the distributions line up, outlier samples, zero counts
# par(mfrow=c(1,1))
# i <- 1
# plot(density(log2(.1+datExpr[,i])), main = 'Scaled Transcript Counts', xlab="log2(counts+.1)", xlim=c(-10,30), ylim = c(0,1))
# for(i in 2:ncol(datExpr)){
#   lines(density(log2(.1+datExpr[,i])))
# }


##### 2
# Remove lowly expressed
datExpr.tpm <- fread(args$gene_tpm, data.table=F)
names(datExpr.tpm)[647] <- "1707.1"
datExpr.tpm <- datExpr.tpm %>% filter(V1 %in% gtf$V5)
rownames(datExpr.tpm) <- datExpr.tpm$V1
datExpr.tpm <- datExpr.tpm[,-1]

# # Select cutoff
# # cutoff = 0
# cut0_df <- data.frame("subjProportion" = seq(0.1, 1, 0.1),
#                       "numTx" = NA)
# for(i in seq(0.1, 1, 0.1)) {
#   keep <- (rowSums(datExpr.tpm >0)) > i*dim(datExpr.tpm)[2]
#   cut0_df[which(cut0_df$subjProportion==i),"numTx"] <- sum(keep)
# }
# cut0_df$cutoff <- rep(0, 10)

# # cutoff = 0.1
# cut1_df <- data.frame("subjProportion" = seq(0.1, 1, 0.1),
#                       "numTx" = NA)
# for(i in seq(0.1, 1, 0.1)) {
#   keep <- (rowSums(datExpr.tpm >.1)) > i*dim(datExpr.tpm)[2]
#   cut1_df[which(cut1_df$subjProportion==i),"numTx"] <- sum(keep)
# }
# cut1_df$cutoff <- rep(.1, 10)

# # cutoff = 0.5
# cut5_df <- data.frame("subjProportion" = seq(0.1, 1, 0.1),
#                       "numTx" = NA)
# for(i in seq(0.1, 1, 0.1)) {
#   keep <- (rowSums(datExpr.tpm >.5)) > i*dim(datExpr.tpm)[2]
#   cut5_df[which(cut5_df$subjProportion==i),"numTx"] <- sum(keep)
# }
# cut5_df$cutoff <- rep(.5, 10)

# # cutoff = 1
# cut10_df <- data.frame("subjProportion" = seq(0.1, 1, 0.1),
#                       "numTx" = NA)
# for(i in seq(0.1, 1, 0.1)) {
#   keep <- (rowSums(datExpr.tpm > 1)) > i*dim(datExpr.tpm)[2]
#   cut10_df[which(cut10_df$subjProportion==i),"numTx"] <- sum(keep)
# }
# cut10_df$cutoff <- rep(1, 10)


# cut_df <- rbind(cut0_df, cut1_df, cut5_df, cut10_df)
# cut_df$cutoff <- as.factor(cut_df$cutoff)

# library(ggplot2)
# p <- ggplot(cut_df, aes(x=subjProportion, y=numTx, color=cutoff)) +
#   geom_point() +
#   labs(x="Proportion of subjects needed with TPM greater than cutoff",
#        y="Number of transcript passing cutoff",
#        title="Transcript TPM cutoff\n(total 220872 transcripts)",
#        color="Cutoff") +
#   theme_light() +
#   theme(axis.text = element_text(size=12),
#         axis.title = element_text(size=14),
#         plot.title = element_text(size=16, face="bold", hjust = 0.5)) +
#   scale_x_continuous(breaks = seq(0,1,0.1)) +
#   ylim(0,200000)

# p
# ggsave("~/Desktop/final/figures/iso/03-isoCutoff.png", p, width = 5.5, height = 5)


# cutoff=0.1, 25% subjects
keep <- (rowSums(datExpr.tpm > .1)) > 0.25*ncol(datExpr.tpm)
# table(keep)

# # Count filter, use cpm, not raw or scaled from TPM counts!
# keep <- filterByExpr(datExpr, min.count = 10, min.prop=0.25)

datExpr <- datExpr[keep,]

# histogram for filtered tx
# i <- 1
# plot(density(log2(.1+datExpr[,i])), main = "Scaled Counts of\nLow Expression Filtered Transcripts", xlab="log2(counts+.1)", xlim=c(-10,30), ylim=c(0,1))
# for(i in 2:ncol(datExpr)){
#   lines(density(log2(.1+datExpr[,i])))
# }


##### 3
# Normalization
# # 1. TMM, cpm
# counts <-  DGEList(datExpr)
# counts <- calcNormFactors(counts, method = 'TMM')  #Perform Trimmed Mean of M-value (TMM) library size normalization
# datExpr.cpm <- as.data.frame(cpm(counts, log=T))
#
# i <- 1
# plot(density(datExpr.cpm[,i]), main = "TMM Nomalizaed Gene Counts", xlab="log2(CPM)", xlim = c(-10, 30), ylim = c(0,0.3))
# for(i in 2:ncol(datExpr.cpm)){
#   lines(density(datExpr.cpm[,i]))
# }

# 2. variance stabilizing transformation (normalize by depth, and stabilize variance)

# before variance stabilized
# library("vsn")
# library(hexbin)
#
# meanSdPlot(as.matrix(log2(.1+datExpr)))
# meanSdPlot(as.matrix(datExpr.cpm))

# Because there's no design matrix, input is a count matrix not a DESeqDataSet
# blind = T or F is no difference
# vst is using a subset of genes to estimate dispersion, faster
# use varianceStabilizingTransformation
# datExpr.vst <- vst(as.matrix(datExpr), blind = TRUE)
# datExpr.vst2 <- vst(as.matrix(datExpr), blind = FALSE)
datExpr.vst <- varianceStabilizingTransformation(as.matrix(datExpr), blind = TRUE)
# meanSdPlot(as.matrix(datExpr.vst))
datExpr.vst <- as.data.frame(datExpr.vst)

# Density after VST normalization
# i <- 1
# plot(density(datExpr.vst[,i]), main = "VST Normalized Transcript Counts", xlab="Counts (in log2)", xlim=c(-10,30), ylim=c(0,1))
# for(i in 2:ncol(datExpr.vst)){
#   lines(density(datExpr.vst[,i]))
# }


# check depth normalization
# par(mfrow=c(2, 2))
# boxplot(log2(.1+datExpr), range=0, main = "log2(gene_counts+.1)", xlab="Samples",xaxt="n")
# boxplot(datExpr.tpm, range=0, main = "Gene TPM", xlab="Samples",xaxt="n")
# boxplot(datExpr.cpm, range=0, main = "TMM normalized, logCPM", xlab="Samples",xaxt="n")
# boxplot(datExpr.vst, range=0, main = "VST normalized, log2(counts)", xlab="Samples",xaxt="n")

# # plot(colSums(datExpr), main = "Gene counts", xlab= "Samples")
# plot(colSums(log2(.1+datExpr)), main = "log2(tx_counts+.1)", xlab="Samples")
# plot(colSums(datExpr.tpm), main = "Tx TPM", xlab="Samples")
# # plot(colSums(2**datExpr.cpm), main = "TMM normalized, CPM", xlab = "Samples")
# plot(colSums(2**datExpr.vst), main = "VST normalized counts", xlab = "Samples")


##### 4
# Remove outliers
normadj <- adjacency(datExpr.vst,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity   #Extract connectivity of each sample
Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score

# par(mfrow=c(1,1))
# plot(1:length(Z.C),Z.C,main="Outlier Plot",xlab = "Samples",ylab="Connectivity Z Score")
# abline(h=-3, col="red")
outliers <- (Z.C < -3)

datExpr.final <- datExpr.vst
datExpr.final <- datExpr.final[,!outliers]
write.table(datExpr.final, paste0(args$outdir,"tx.counts.processed.noComBat.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")


##### 5
# ComBat
walker <- read.table(paste0(args$geno_subj_dir, "geno_subj_walker.txt"))
obrien <-read.table(paste0(args$geno_subj_dir, "geno_subj_obrien.txt"))
werling <-read.table(paste0(args$geno_subj_dir, "geno_subj_werling.txt"))
hdbr <- read.table(paste0(args$geno_subj_dir, "geno_subj_hdbr.txt"))
libd <- read.table(paste0(args$geno_subj_dir, "geno_subj_libd.txt"))

data.batch <- c()

for (i in 1:ncol(datExpr.final)) {
  sample <- colnames(datExpr.final)[i]
  if (sample %in% walker[,1]) {
    data.batch[i] <- 1
  }
  if (sample %in% obrien[,1]) {
    data.batch[i] <- 2
  }
  if (sample %in% werling[,1]) {
    data.batch[i] <- 3
  }
  if (sample %in% hdbr[,1]) {
    data.batch[i] <- 4
  }
  if (sample %in% libd[,1]) {
    data.batch[i] <- 5
  }
}

exprMat <- as.matrix(datExpr.final)
combat_expr <- ComBat(dat = exprMat, batch = data.batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
combat_expr <- as.data.frame(combat_expr)

write.table(combat_expr, paste0(args$outdir, "tx.counts.processed.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")

outlier.df <- data.frame(c(which(outliers)))
outlier.df$c.which.outliers.. <- rownames(outlier.df)
write.table(outlier.df,paste0(args$outdir,"tx.outlier.txt"), sep="\t", quote=F, col.names = F, row.names = F)


##### 6
# phenotype.bed
gtf <- gtf %>% filter(V5 %in% rownames(datExpr.final))

# prepare BED file
setDT(combat_expr, keep.rownames = TRUE)
bed <- merge(combat_expr, gtf, by.x = "rn", by.y = "V5", all = TRUE)
bed <- as.data.frame(bed)

bed <- bed[, c(ncol(bed)-3, ncol(bed)-2, ncol(bed)-1, ncol(bed), 1:(ncol(bed)-4))]

# special note! strand; 0/1-based conversion
bed$V1 <- gsub("chr","",bed$V1)

for (i in 1:dim(bed)[1]){
	if (bed[i,4] == "+") {
		gtf_start=bed[i,2]
		bed[i,2]=gtf_start-1
		bed[i,3]=gtf_start
	}
	else if (bed[i,4] == "-") {
		gtf_end=bed[i,3]
		bed[i,2]=gtf_end-1
		bed[i,3]=gtf_end
	}
}

bed <- bed[order(bed$V1, bed$V2),]
bed <- bed %>% select(-V4)
colnames(bed)[1:4] <- c("#Chr", "start", "end", "ID")

write.table(bed, paste0(args$outdir, "tx.counts.scaled.normalized.bed"), quote=F, sep="\t", row.names=F, col.names=T)
