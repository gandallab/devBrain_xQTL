#! /usr/bin/env Rscript

library(argparser)
library(data.table)
library(dplyr)

argp <- arg_parser("Fix APEX factors")
argp <- add_argument(argp, "--known", help="known covariates")
argp <- add_argument(argp, "--factor", help="APEX factor")
argp <- add_argument(argp, "--outfile", help="")
argv <- parse_args(argp)

known <- fread(argv$known, data.table = F)
factor <- fread(argv$factor, data.table = F)
# Get inferred factors
factor <- factor[grep("factor",factor$'#id'),]
# Reorder factor columns
factor <- factor[,c('#id',colnames(known)[2:ncol(known)])]
# Cat to known covariates
colnames(known)[1] <- '#id'
cov <- rbind(known, factor)

write.table(cov, argv$outfile, col.names = T, row.names = F, quote = F, sep = "\t")