#! /usr/bin/env Rscript
# NOTE: paintor wants signed pearson correlation r. Don't use this
suppressMessages(library(argparser))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidyverse))

p <- arg_parser("Make paintor LD file")
p <- add_argument(p, "--gene", help="")
args <- parse_args(p)

setwd(paste0("/u/project/gandalm/cindywen/isoform_twas/eqtl_new/paintor/", args$gene, "_dir"))
num_var <- system("wc -l shared_variants.txt", intern = TRUE) %>% strsplit(split = " ")
num_var <- as.numeric(num_var[[1]][1])

# EUR
myFile <- file("eur.ld.bin", "rb")
eur_ld <- readBin(myFile, what = "numeric", n = num_var*num_var, size = 8)
close(myFile)
eur_matrix <- matrix(eur_ld, nrow = num_var)
write.table(eur_matrix, paste0(args$gene, ".LD1"), col.names = F, row.names = F, quote = F, sep = " ")

# AMR
myFile <- file("amr.ld.bin", "rb")
amr_ld <- readBin(myFile, what = "numeric", n = num_var*num_var, size = 8)
close(myFile)
amr_matrix <- matrix(amr_ld, nrow = num_var)
write.table(amr_matrix, paste0(args$gene, ".LD2"), col.names = F, row.names = F, quote = F, sep = " ")

# AFR
myFile <- file("afr.ld.bin", "rb")
afr_ld <- readBin(myFile, what = "numeric", n = num_var*num_var, size = 8)
close(myFile)
afr_matrix <- matrix(afr_ld, nrow = num_var)
write.table(afr_matrix, paste0(args$gene, ".LD3"), col.names = F, row.names = F, quote = F, sep = " ")
