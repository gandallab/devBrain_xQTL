#! /usr/bin/env Rscript
library(qvalue)

args <- commandArgs(trailingOnly=T)
trans_dir <- args[1]

allcoef <- NULL
allp <- NULL
for(chr in c(1:22)){
    dat <- read.table(paste(trans_dir,"chr",chr,"_allp_h2g_filtered_interchrom.txt", sep = ""))
    allp <- c(allp, as.numeric(dat[,1]))
}
test <- qvalue(allp, fdr.levels = 0.1)
pvalues <- test$pvalues
print(max(pvalues[test$qvalues <= 0.1]))
print(max(pvalues[test$qvalues <= 0.05]))