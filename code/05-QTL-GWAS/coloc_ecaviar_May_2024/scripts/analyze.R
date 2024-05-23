#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Analyze eCAVIAR results")
p <- add_argument(p, "--gwas", help = "")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--annot", help = "")
p <- add_argument(p, "--gene", help = "")

args <- parse_args(p)

result <- fread(paste0(
    "../out_1MB/", args$gwas, "/locus", args$locus, "/", args$gene, "_", args$annot, "_ecaviar_col"
), data.table = F)

if (max(result$CLPP > 0.01)) {
    sig <- result %>% filter(CLPP > 0.01)
    sig$GWAS <- args$gwas
    sig$locus <- args$locus
    sig$annot <- args$annot
    sig$gene <- args$gene
    write.table(sig, paste0(
        "../out_1MB/", args$gwas, "/locus", args$locus, "/", args$gene, "_", args$annot, "_ecaviar_col_sig.txt"
    ),
    col.names = T, row.names = F, quote = F, sep = "\t"
    )
}
