#! /usr/bin/env Rscript
library(argparser)
library(data.table)
library(dplyr)

p <- arg_parser("Make expr input for susie")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--outname", help="")
args <- parse_args(p)

dat <- fread(args$input, data.table = F)
colnames(dat)[1] <- "phenotype_id"
write.table(dat, args$outname, col.names=T, row.names=F, quote=F, sep="\t")

