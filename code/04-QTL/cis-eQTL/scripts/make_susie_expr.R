#! /usr/bin/env Rscript
library(argparser)
library(data.table)
library(dplyr)

p <- arg_parser("Make expr input for susie")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--outname", help="")
p <- add_argument(p, "--exclude_rel", help="List of relatives to remove")
args <- parse_args(p)

dat <- fread(args$input, data.table = F)
rel <- read.table(args$exclude_rel, header = F, stringsAsFactors = F)
colnames(dat)[1] <- "phenotype_id"
dat <- dat[,!colnames(dat) %in% rel$V1]
write.table(dat, args$outname, col.names=T, row.names=F, quote=F, sep="\t")

