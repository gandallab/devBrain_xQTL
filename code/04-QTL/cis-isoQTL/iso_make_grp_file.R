#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(tidyr)

p <- arg_parser("Write isoQTL group file for FastQTL mapping with grouping")
p <- add_argument(p, "--tx2gene", help="")
p <- add_argument(p, "--out", help="")

argv <- parse_args(p)

tx2gene <- fread(argv$tx2gene, data.table = F)
tx2gene <- tx2gene %>% unite("Tx", Tx, Gene, sep = ":", remove = FALSE)
write.table(tx2gene, argv$out, col.names = F, row.names = F, quote = F, sep = "\t")