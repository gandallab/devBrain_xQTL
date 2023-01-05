#! /usr/bin/env Rscript

library(data.table)
library(dplyr)
library(argparser)

p <- arg_parser("Get isoform in locus")
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
    isoform <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri1_nominal_35hcp/significant_assoc.txt", data.table = F)
}

if(args$group == "tri2") {
    isoform <- fread("/u/project/gandalm/cindywen/isoform_twas/isoqtl_new/results/tri2_nominal_20hcp/significant_assoc.txt", data.table = F)
}

isoform_in_locus <- isoform %>% filter(sid %in% locus_file$V1)
isoform_list <- as.data.frame(unique(isoform_in_locus$pid))
write.table(isoform_list, paste0(
    "/u/project/gandalm/cindywen/isoform_twas/colocal/results_", args$group, "_isoqtl/", args$trait, "/locus_", args$locus, "/locus_isoform.txt"
), col.names = F, row.names = F, quote = F, sep = "\t")

