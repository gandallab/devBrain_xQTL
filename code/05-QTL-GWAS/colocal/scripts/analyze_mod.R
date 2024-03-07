#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Analyze eCAVIAR results")
p <- add_argument(p, "--level", help = "genemod, isomod")
p <- add_argument(p, "--trait", help = "")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--module_gene", help = "modxx_genexx")

args <- parse_args(p)

# Read CLPP
result <- fread(paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_",
    args$level, "_ieqtl/",
    args$trait, "/locus_", args$locus, "/", 
    args$module_gene, "_ecaviar_col"
), data.table = F)

if (max(result$CLPP > 0.01)) {
    sig <- result %>% filter(CLPP > 0.01)
    sig$locus <- args$locus
    sig$module_gene <- args$module_gene
    write.table(sig, paste0(
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_",
        args$level, "_ieqtl/",
        args$trait, "/locus_", args$locus, "/", 
        args$module_gene, "_ecaviar_col_sig.txt"
    ),
    col.names = T, row.names = F, quote = F, sep = "\t"
    )
}

