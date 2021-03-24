#! /usr/bin/env Rscript

library(argparser)
library(Matrix)

argp <- arg_parser("Make GRM for apex from GENESIS pcrelate")
argp <- add_argument(argp, "--pcrel", help="pcrelate RData, GRM")
argp <- add_argument(argp, "--out", help="Output file, in apex format")
argv <- parse_args(argp)

# Matrix.RData is made using GENESIS::pcrelateToMatrix with scaleKin = 2

load(argv$pcrel)
df <- as.data.frame(summary(km))
for(row in 1:nrow(df)){
df[row,"i"] <- colnames(km)[as.numeric(df[row,"i"])]
df[row,"j"] <- colnames(km)[as.numeric(df[row,"j"])]}
colnames(df) <- c("id1","id2","kinship")

write.table(df, argv$out, col.names=T, row.names=F, quote=F, sep=" ")
