#! /usr/bin/env Rscript

library(data.table)
library(argparser)
p <- arg_parser("Generate variant coord file")
p <- add_argument(p, "--info", help="variant info")
p <- add_argument(p, "--out", help="output file")
args <- parse_args(p)

dat <- fread(args$info, data.table = F)
colnames(dat) <- c("CHR", "START", "GENE")
dat$END <- dat$START+1
dat <- dat[dat$CHR %in% c(1:22),]
dat <- dat[,c("GENE","CHR","START","END")]
write.table(dat, file = args$out, col.names = T, row.names = F, quote = F, sep = "\t")
