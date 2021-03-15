#! /usr/bin/env Rscript

library(argparser)
library(TopmedPipeline)
library(SeqArray)
library(SNPRelate)
sessionInfo()

argp <- arg_parser("Convert GDS to BED")
argp <- add_argument(argp, "--gds_file", help="GDS file, all chr")
argp <- add_argument(argp, "--bed_file", help="BED file, all chr")
argv <- parse_args(argp)

gds <- seqOpen(argv$gds_file)
snpfile <- tempfile()
seqGDS2SNP(gds, snpfile)
seqClose(gds)

gds <- snpgdsOpen(snpfile)
snpgdsGDS2BED(gds, argv$bed_file)
snpgdsClose(gds)

unlink(snpfile)
