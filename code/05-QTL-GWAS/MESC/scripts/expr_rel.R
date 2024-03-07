#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(dplyr)
p <- arg_parser("Remove relatives in expression files. Used --exclude-samples in FastQTL")
p <- add_argument(p, "--rel", help="")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

dat <- fread(args$input, data.table = F)
rel <- read.table(args$rel, header = F, stringsAsFactors = F)
rel <- rel %>% filter(V1 %in% colnames(dat))
dat <- dat %>% select(-rel$V1)
write.table(dat, args$out, col.names = T, row.names = F, quote = F, sep = "\t")