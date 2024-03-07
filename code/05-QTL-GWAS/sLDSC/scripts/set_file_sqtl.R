#! /usr/bin/env Rscript

library(argparser)
library(dplyr)
library(data.table)
p <- arg_parser("Make set file, top sQTL per sGene")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

dat <- fread(args$input, data.table = F)


# top sQTL: per sGene, QTL with the best bpval
dat1 <- dat %>% arrange(ensg, bpval)
# temporary workaround: remove introns that are not mapped to genes (with GTEx's pipeline, only run FastQTL for introns that are mapped to genes)
# note there are introns mapped to genes but has NA transcripts, criptic unanchored
dat1 <- dat1 %>% filter(substring(dat1$ensg, 1, 4) == "ENSG")
top <- dat1[!duplicated(dat1$ensg),]
df <- data.frame(unique(top$sid))
write.table(df, args$out, col.names = F, row.names = F, sep = "\t", quote = F)
