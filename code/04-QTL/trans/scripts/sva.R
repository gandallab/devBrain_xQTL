#! /usr/bin/env Rscript
library(sva)

args <- commandArgs(trailingOnly=TRUE)
input <- args[1]
output <- args[2]
num_sv <- as.numeric(args[3])
tool_dir <- args[4]

source(paste(tool_dir, 'make_sva.R', sep=""))
load(input)
sv <- make_sva(ex, num_sv)
write.table(sv, paste(output, sep=""), quote = F, row.names = F)