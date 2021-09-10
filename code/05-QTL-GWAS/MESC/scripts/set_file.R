#! /usr/bin/env Rscript
library(argparser)
library(data.table)
library(dplyr)
p <- arg_parser("Make gene set file for MESC expression score estimation")
p <- add_argument(p, "--input", help="FastQTL significant features")
p <- add_argument(p, "--out", help="")
p <- add_argument(p, "--name", help="Set name")
args <- parse_args(p)

dat <- fread(args$input, data.table = F)
dat <- data.frame(unique(dat$pid))
dat <- t(dat)
rownames(dat) <- args$name
write.table(dat, args$out, col.names = F, row.names = T, quote = F, sep = "\t")