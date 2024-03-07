#! /usr/bin/env Rscript
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

p <- arg_parser("Make paintor annotation file")
p <- add_argument(p, "--gene", help="")
args <- parse_args(p)

setwd(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/", args$gene, "_dir"))
variant <- read.table("shared_variants.txt", header = F, stringsAsFactors = F)
colnames(variant)[1] <- "RSID"
annot <- fread("/u/project/gandalm/cindywen/isoform_twas/genotype/all_data/isec_R2_greater_than_3/ancestry/annot/eur_variant_annot.txt.gz", data.table = F)

variant <- variant %>% left_join(annot, by=c("RSID"="SNP"))
locus <- variant %>% select(-1)
write.table(locus, paste0(args$gene, ".annotations"), col.names = T, row.names = F, quote = F, sep = " ")
