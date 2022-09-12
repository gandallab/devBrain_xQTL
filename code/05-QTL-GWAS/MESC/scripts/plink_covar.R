#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(dplyr)
p <- arg_parser("Make plink covar, number of subjects is final for QTL mapping")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

dat <- read.table(args$input, stringsAsFactors = F, header = T, check.names = F)
new <- data.frame(t(dat[,-1]))
colnames(new) <- dat$id
new <- cbind(data.frame(rownames(new)), data.frame(rownames(new)), new)
colnames(new)[1:2] <- c("FID", "IID")
# 8/30/22: it seems this is not supported by R/4.1.0; for features other than gene, iso, intron, manually generated MESC covar files
new$sex <- as.numeric(new$sex) # F=1, M=2
write.table(new, args$out, row.names = F, col.names = T, quote = F, sep = "\t")
