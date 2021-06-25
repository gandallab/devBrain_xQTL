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
p <- add_argument(p, "--level", help="gene or tx")
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


##### 2
# Remove lowly expressed
datExpr.tpm <- fread(args$gene_tpm, data.table=F)
names(datExpr.tpm)[647] <- "1707.1"
datExpr.tpm <- datExpr.tpm %>% filter(V1 %in% gtf$V5)
rownames(datExpr.tpm) <- datExpr.tpm$V1
datExpr.tpm <- datExpr.tpm[,-1]

# cutoff=0.1, 25% subjects
keep <- (rowSums(datExpr.tpm > .1)) > 0.25*ncol(datExpr.tpm)
# table(keep)

# # Count filter, use cpm, not raw or scaled from TPM counts!
# keep <- filterByExpr(datExpr, min.count = 10, min.prop=0.25)

datExpr <- datExpr[keep,]


##### 3
# Normalize
datExpr.vst <- varianceStabilizingTransformation(as.matrix(datExpr), blind = TRUE)
datExpr.vst <- as.data.frame(datExpr.vst)


#### 4
# Remove outliers
normadj <- adjacency(datExpr.vst,type = 'signed',corFnc = 'bicor')   #Calculate network adjacency
netsummary <- fundamentalNetworkConcepts(normadj)
C <- netsummary$Connectivity   #Extract connectivity of each sample
Z.C <- (C-mean(C))/sqrt(var(C))   #Covert to Z-score
outliers <- (Z.C < -3)

datExpr.final <- datExpr.vst
datExpr.final <- datExpr.final[,!outliers]
write.table(datExpr.final, paste0(args$outdir,args$level,".counts.processed.noComBat.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")


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

write.table(combat_expr, paste0(args$outdir, args$level, ".counts.processed.tsv"), col.names = T, row.names = T, quote = F, sep = "\t")

outlier.df <- data.frame(c(which(outliers)))
outlier.df$c.which.outliers.. <- rownames(outlier.df)
write.table(outlier.df,paste0(args$outdir, args$level, ".outlier.txt"), sep="\t", quote=F, col.names = F, row.names = F)


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

write.table(bed, paste0(args$outdir, args$level, ".counts.scaled.normalized.bed"), quote=F, sep="\t", row.names=F, col.names=T)