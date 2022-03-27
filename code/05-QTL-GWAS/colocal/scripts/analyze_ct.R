#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Analyze eCAVIAR results")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--gene", help = "")
p <- add_argument(p, "--trait", help = "")
p <- add_argument(p, "--type", help = "")
p <- add_argument(p, "--hcp", help = "")

args <- parse_args(p)

# Read CLPP
result <- fread(paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_",
    args$type, "_", args$hcp, "hcp_eqtl/",
    args$trait, "/locus_", args$locus, "/", args$gene, "_ecaviar_col"
), data.table = F)

if (max(result$CLPP > 0.01)) {
    sig <- result %>% filter(CLPP > 0.01)
    sig$locus <- args$locus
    sig$gene <- args$gene
    sig$type <- args$type
    write.table(sig, paste0(
        "/u/project/gandalm/cindywen/isoform_twas/colocal/results_",
        args$type, "_", args$hcp, "hcp_eqtl/",
        args$trait, "/locus_", args$locus, "/", args$gene,
        "_ecaviar_col_sig.txt"
    ),
    col.names = T, row.names = F, quote = F, sep = "\t"
    )
}
