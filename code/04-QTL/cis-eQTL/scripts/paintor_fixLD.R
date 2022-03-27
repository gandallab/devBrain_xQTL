#! /usr/bin/env Rscript
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

p <- arg_parser("Exclude NA in LD files")
p <- add_argument(p, "--gene", help = "")
args <- parse_args(p)

setwd(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/",
    args$gene, "_dir"))

# Only AFR has NA in LD file
eur <- fread(paste0(args$gene, ".LD1"), data.table = F)
amr <- fread(paste0(args$gene, ".LD2"), data.table = F)
afr <- fread(paste0(args$gene, ".LD3"), data.table = F)
zscore <- read.table(args$gene, header = T)
annot <- read.table(paste0(args$gene, ".annotations"), header = T)
variants <- read.table("shared_variants.txt")

vec <- which(zscore$ZSCORE.P3 == 0)
# note if vec is integer(0), the following code will make LD matrices empty
eur <- eur[-vec, -vec]
amr <- amr[-vec, -vec]
afr <- afr[-vec, -vec]
zscore <- zscore[-vec, ]
annot <- annot[-vec, ]
variants <- variants[-vec, ]

write.table(eur, paste0(args$gene, ".LD1"),
    col.names = F, row.names = F, quote = F, sep = " ")
write.table(amr, paste0(args$gene, ".LD2"),
    col.names = F, row.names = F, quote = F, sep = " ")
write.table(afr, paste0(args$gene, ".LD3"),
    col.names = F, row.names = F, quote = F, sep = " ")
write.table(zscore, args$gene,
    col.names = T, row.names = F, quote = F, sep = " ")
write.table(annot, paste0(args$gene, ".annotations"),
    col.names = T, row.names = F, quote = F, sep = " ")
write.table(variants, "shared_variants.txt",
    col.names = F, row.names = F, quote = F, sep = "\t")
