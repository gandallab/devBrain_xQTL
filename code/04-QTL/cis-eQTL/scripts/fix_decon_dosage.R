#! /usr/bin/env Rscript
library(data.table)
library(argparser)

p <- arg_parser("Fix dosage file for Decon-QTL")
p <- add_argument(p, "--dosage", help="Dosage file to fix")
p <- add_argument(p, "--sample", help="Samples, or header of the dosage file to fix")
p <- add_argument(p, "--cellcount", help="Dosage file header order needs to match that of cell count and expression file")
p <- add_argument(p, "--out", help="fixed dosage file")
args <- parse_args(p)

# Read dosage file
dosage <- fread(args$dosage, data.table = F)
rownames(dosage) <- dosage$V1
dosage <- dosage[,-1]

# Column names (subjects) are in the order of genotype files
sample <- read.table(args$sample, header = F, stringsAsFactors = F )
colnames(dosage) <- sample$V1

# Column names (i.e. subject IDs) of dosage, cell count, and expression files should be in the same order
cg <- fread(args$cellcount, data.table = F)
dosage <- dosage[,cg$V1]

# Write to output, note SNP IDs should be rownames
write.table(dosage, args$out, col.names = T, row.names = T, quote = F, sep = "\t")
