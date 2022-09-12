#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Get eGene in locus")
p <- add_argument(p, "--locus", help = "")
p <- add_argument(p, "--trait", help = "")
p <- add_argument(p, "--group", help = "male, female, tri1, tri2")

args <- parse_args(p)

locus_file <- read.table(paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_eqtl/",
    args$trait, "/locus_", args$locus, "/gwas_in_bim.txt"
),
header = F, stringsAsFactors = F
)

if(args$group == "male") {
    egene <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/ALL_m_nominal_50hcp/significant_assoc.txt", data.table = F)
}

if(args$group == "female") {
    egene <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/ALL_f_nominal_50hcp/significant_assoc.txt", data.table = F)
}

if(args$group == "tri1") {
    egene <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T1-nominal_significant_assoc.txt", data.table = F)
}

if(args$group == "tri2") {
    egene <- fread("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/results/eur_trimester/T2-nominal_significant_assoc.txt", data.table = F)
}

egene_in_locus <- egene %>% filter(sid %in% locus_file$V1)
egene_list <- as.data.frame(unique(egene_in_locus$pid))
write.table(egene_list, paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_", args$group, "_eqtl/", args$trait, "/locus_", args$locus, "/locus_egene.txt"
), col.names = F, row.names = F, quote = F, sep = "\t")
