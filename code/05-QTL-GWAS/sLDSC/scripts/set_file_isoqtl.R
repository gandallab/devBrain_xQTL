#! /usr/bin/env Rscript

library(argparser)
library(dplyr)
p <- arg_parser("Make set file, top isoQTL per isoGene")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--out", help="")
p <- add_argument(p, "--tx2gene", help="")
p <- add_argument(p, "--out2", help="")
args <- parse_args(p)

dat <- read.table(args$input, header=T, stringsAsFactors=F)
tx2gene <- read.table(args$tx2gene, header=T, stringsAsFactors=F)
dat2 <- dat %>% left_join(tx2gene, by=c("pid"="Tx"))
write.table(dat2, args$out2, col.names=T, row.names=F, sep="\t", quote=F)

# top isoQTL: per isoGene, QTL with the best bpval
dat3 <- dat2 %>% arrange(Gene, bpval)
top <- dat3[!duplicated(dat3$Gene),]
df <- data.frame(unique(top$sid))
write.table(df, args$out, col.names=F, row.names=F, sep="\t", quote=F)
