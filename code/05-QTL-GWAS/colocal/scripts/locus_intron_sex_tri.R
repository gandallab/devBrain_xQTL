#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Get intron in locus")
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

# if(args$group == "male") {
#     egene <- fread("~/project-gandalm/isoform_twas/sqtl_new/results/outputs_fetal_ALL_sex_specific_sqtl_male/", data.table = F)
# }

# if(args$group == "female") {
#     egene <- fread("~/project-gandalm/isoform_twas/sqtl_new/results/outputs_fetal_ALL_sex_specific_sqtl_female/", data.table = F)
# }

if(args$group == "tri1") {
    intron <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri1_nominal_15hcp/significant_assoc.txt", data.table = F)
}

if(args$group == "tri2") {
    intron <- fread("/u/project/gandalm/cindywen/isoform_twas/sqtl_new/results/tri2_nominal_10hcp/significant_assoc.txt", data.table = F)
}

intron_in_locus <- intron %>% filter(sid %in% locus_file$V1)
intron_list <- as.data.frame(unique(intron_in_locus$pid))
write.table(intron_list, paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_", args$group, "_sqtl/", args$trait, "/locus_", args$locus, "/locus_intron.txt"
), col.names = F, row.names = F, quote = F, sep = "\t")

