#! /usr/bin/env Rscript

library(argparser)
library(dplyr)
library(data.table)
p <- arg_parser("Make set file, top iso/sQTL per iso/sGene, GTEx grouped permutation")
p <- add_argument(p, "--input", help="")
p <- add_argument(p, "--out", help="")
args <- parse_args(p)

# dat has top association per phenotype group (gene)
# a few duplicated feature (tx, intron) and group feature (gene) exist
dat <- fread(args$input, data.table = F)
dat <- dat %>% filter(qval < 0.05)
df <- data.frame(unique(dat$variant_id))
write.table(df, args$out, col.names = F, row.names = F, sep = "\t", quote = F)