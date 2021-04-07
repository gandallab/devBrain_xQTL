#! /usr/bin/env Rscript

library(argparser)
p <- arg_parser("Make set file")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

dat <- read.table(args$input, header=T, stringsAsFactors=F)
df <- as.data.frame(unique(dat$sid))
write.table(df, args$out, col.names=F, row.names=F, sep="\t", quote=F)
