#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(tidyr)

argp <- arg_parser("")
argp <- add_argument(argp, "--input", help="")
argp <- add_argument(argp, "--output", help="")
argv <- parse_args(argp)

dat <- fread(argv$input, data.table=F)
dat <- dat %>% separate(chrom, c("chr","start","end","clu"), sep="[:]")
df <- data.frame(unique(dat$chr))
colnames(df) <- "chr"
remove <- data.frame("chr"=df[grepl(pattern="GL|X|Y|M", x=df$chr),"chr"])
write.table(remove, argv$output, col.names=F, row.names=F, quote=F, sep="\t")
