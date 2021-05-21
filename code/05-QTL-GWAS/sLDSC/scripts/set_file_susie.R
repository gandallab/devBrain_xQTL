#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(dplyr)
p <- arg_parser("Make set file for susie output")
p <- add_argument(p, "--susie", help="susie finemapping output, all variants that belong to CS")
p <- add_argument(p, "--info", help="variant info")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

dat <- fread(args$susie, data.table = F)
info <- fread(args$info, data.table = F)
dat <- dat%>%left_join(info, by=c("chr"="V1","pos"="V2"))
df <- as.data.frame(unique(dat$V3))
write.table(df, args$out, col.names=F, row.names=F, sep="\t", quote=F)
